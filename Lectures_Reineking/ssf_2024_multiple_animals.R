#' # Step selection functions: multiple animals
#' Björn Reineking, 2024-06-20
#'
#' We will use the buffalo data set.
#' The code in this example is based on (slightly modified) from Muff, Signer & Fieberg (2020) J Anim Ecol 89: 80-92.
#'
#' # Loading packages
#+ results='hide', message=FALSE, warning=FALSE
library(animove)
library(dplyr)
library(amt)
library(glmmTMB)
library(sf)

#' Estimating mixed effects models for step selection functions is not trivial.
#' The computations can take a long time, and the model may not converge.
#'
#' Also, we do not always need mixed effects models when analysing movement data for multiple animals.
#' Often we have enough data for each animal to estimate the movement parameters.
#' In a mixed effects model, we make the
#' assumption that the variation in parameter values across individuals follows a normal distribution.
#' This may not be appropriate.
#' The advantage of the mixed effects model is that we get estimates of the average population-level effect and
#' a measure of the variation between individuals.
#'
#' If the data are sufficient to fit a model to each individual, I would recommend fitting individual models as well and looking
#' at the distribution of the parameter estimates, to see if they are reasonably close to a normal distribution.
#'

#' Load data
data(buffalo_utm_mv2)
#' Convert MoveStack to amt::track object.
#'
#' Currently, there is no conversion function from move::moveStack to amt::track implemented, so we do it by hand (thank you, Anne!)
move2_TO_track.xyt <- function(mv2){
  if(mt_is_move2(mv2)){
    warning("!!INFO!!: only coordinates, timestamps and track IDs are retained")
    track(
      x=sf::st_coordinates(mv2)[,1],
      y=sf::st_coordinates(mv2)[,2],
      t=mt_time(mv2),
      id=mt_track_id(mv2),
      crs = sf::st_crs(mv2)
    )
  }
}

buffalo_tracks <- move2_TO_track.xyt(buffalo_utm_mv2)

#' Environmental data: topography, waterways, and NDVI
data(buffalo_env)
#' We exclude the full first day, since at least for buffalo Cilla there are some weird locations
#' We perform the data preparation separately for each individual, but we use the same proposal distribution for step length and turning angles.
#' ## Thin movement data and split to bursts
#'
#' - We reduce the data set to observations that are within a certain time step range. The SSF assumes Brownian motion, so we should thin sufficiently, so that the velocities of successive steps are uncorrelated. Here we go for 3 hours.
#' - There is some tolerance around the target time interval of 3 hours. When two observations are separated by less than the threshold, the second observation is removed
#' - When two observations are separated by more than the upper threshold, the observations are assigned to different bursts.
#'
step_duration <- 3

buffalos <- buffalo_tracks |> nest(track = -"id")
buffalos <- buffalos |> 
  mutate(steps = map(track, function(x) { 
    x |> dplyr::filter(t_ > min(t_) + days(1)) |> 
      track_resample(hours(step_duration), tolerance = minutes(15)) |>
      filter_min_n_burst(3) |>
      steps_by_burst()})) |> 
  dplyr::select(id, steps)


#' The estimated parameters related to the movement kernel (e.g. steplength,
#' log(steplength), cos(turnangle) are relative to the proposal distribution.
#' 
#' We have to ensure that we use the same proposal distribution for all animals,
#' because otherwise the movement kernel parameters would mean (slightly) different
#' things for different animals and it would not make sense to assume that
#' the estimates come from a common distribution
#' 
steps_all <- buffalos |> unnest(cols = steps)
combined_sl_distr <- fit_distr(steps_all$sl_, "gamma")
combined_ta_distr <- fit_distr(steps_all$ta_, "vonmises")

#' Create random steps, using the same proposal distribution for all individuals.
#' We use only 100 control steps to speed up calculations.
#' 
set.seed(23)
ssf_buffalos <- buffalos |> 
  mutate(steps = map(steps, random_steps, n_control = 100, 
                     sl_distr = combined_sl_distr,
                     ta_distr = combined_ta_distr)) |> 
  unnest(cols = steps)

#' We need to modify the values for burst_ and step_id_ to ensure that not 
#' two individuals are assigned the same values for step_id_
#' 
ssf_buffalos <- ssf_buffalos |> 
  mutate(burst_ = paste(id, burst_), step_id_ = paste(id, step_id_))
table(table(ssf_buffalos$step_id_)) # Now OK. step_id_ is used in the statistical modelling, so we need to make sure that the right observations are combined.

#' Extract environmental covariates
ssf_buffalos <- extract_covariates(ssf_buffalos, buffalo_env)
#' Remove NA's
ssf_buffalos <- ssf_buffalos[complete.cases(ssf_buffalos),]

#' Convert to standard data frame
ssf_df <- as.data.frame(ssf_buffalos)

#' Remove steps where there is no observed step (i.e. there is no case_ == TRUE) 
ssf_df <- ssf_df |> dplyr::group_by(step_id_) |> dplyr::filter(any(case_))

#' ## Scaling predictors
#' We scale the predictors (mean of zero and standard deviation of one) 
#' because this helps with model convergence.
#' The values of the parameter estimates depend on the scaling of the predictors. 
#' So when we want to make predictions with the model,
#' we have to transform the maps of the predictors in the same way. 
#' The parameters for this transformation are stored in
#' the output of the scale function.
#'
scaled_env <- scale(ssf_df[,c("elev", "slope", "water_dist", "mean_NDVI", 
                              "var_NDVI")])
colnames(scaled_env) <- paste0(colnames(scaled_env), "_scaled")
ssf_df <- cbind(ssf_df, scaled_env)

#' ## Model without taking individual variation into account
m_3 <- fit_clogit(ssf_df, case_ ~ cos(ta_) + sl_ + log(sl_) +
                    elev_scaled + mean_NDVI_scaled + strata(step_id_))
summary(m_3)

#' ## Fit separate parameter values for each individual
m_4 <- fit_clogit(ssf_df, 
                  case_ ~ (cos(ta_) + sl_ + log(sl_) + elev_scaled + 
                             mean_NDVI_scaled):id + strata(step_id_))
summary(m_4)

#' Alternatively, we can fit a model for each animal and then e.g. extract the 
#' parameter estimates
buffalo_models <- ssf_df |> nest(data = -"id") |> 
  mutate(model = map(data, fit_clogit, case_ ~ cos(ta_) + sl_ + log(sl_) +
                       elev_scaled + mean_NDVI_scaled + strata(step_id_)))
buffalo_coef <- buffalo_models |> 
  mutate(coef = map(model, function(x) as.data.frame(t(coef(x))))) |>
  dplyr::select(id, coef) |> 
  unnest(col = coef)

buffalo_coef

#' ## Fixed-effects model glmmTMB
#' To see if we get (more or less) the same result with the glmmTMB library
#' for the m_3 equivalent model.
#'
#' I generally had better convergence when using the optim optimizer 
#' and method BFGS
TMBStruc.fix = glmmTMB(case_ ~ -1 + cos(ta_) + sl_ + log(sl_) +
                           elev_scaled + mean_NDVI_scaled  + (1|step_id_),
                         family = poisson, data = ssf_df, doFit = FALSE,
                         control = glmmTMBControl(optimizer = optim, 
                                                  optArgs = list(method="BFGS")))
TMBStruc.fix$parameters$theta[1] = log(1e3)
TMBStruc.fix$mapArg = list(theta = factor(c(NA)))
TMBStruc.fix <- glmmTMB::fitTMB(TMBStruc.fix)
# glmmTMB::diagnose(TMBStruc.fix)
#' diagnose does not currently work
#' Error in names(pp) <- nn : 
#' 'names' attribute [6] must be the same length as the vector [5]
#' 

summary(TMBStruc.fix)

#' Estimates are similar.
fixef(TMBStruc.fix)
coef(m_3)

#' Confidence intervals as well.
confint(TMBStruc.fix)
confint(m_3$model)

#' ## Random effects model
#'
#' The variable representing different animals, ANIMAL_ID, needed to be numeric in earlier versions of glmmTMB;
ssf_df$ANIMAL_ID <- factor(ssf_df$id)

TMBRandomStruc = glmmTMB(case_ ~ -1 + cos(ta_) + sl_ + log(sl_) + elev_scaled + mean_NDVI_scaled +
                           (1|step_id_) +
                           (0 + cos(ta_) | ANIMAL_ID) +
                           (0 + sl_ | ANIMAL_ID) +
                           (0 + log(sl_) | ANIMAL_ID) +
                           (0 + elev_scaled | ANIMAL_ID) +
                           (0 + mean_NDVI_scaled | ANIMAL_ID), family = poisson,
                         data=ssf_df, doFit=FALSE,
                         control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

TMBRandomStruc$parameters$theta[1] = log(1e3)

#' Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
TMBRandomStruc$mapArg = list(theta=factor(c(NA,1:5)))

#' Fit the model and look at the summary:
glmm.TMB.random <- glmmTMB::fitTMB(TMBRandomStruc)

# diagnose(glmm.TMB.random)
# diagnose does not currently work:
# Error in names(pp) <- nn : 
#   'names' attribute [11] must be the same length as the vector [10]
summary(glmm.TMB.random)
#' So both elevation and mean_NDVI are significant at the population level (animals avoid higher elevation sites and prefer higher NDVI).
#'

#' 95\% CIs for fixed and random effects (standard deviations) are obtained via the confint() function:
confint(glmm.TMB.random)
#' The standard deviation of the random effect is large relative to the fixed effect size,
#' so there can be animals that react positively to elevation or negatively to NDVI
#'


