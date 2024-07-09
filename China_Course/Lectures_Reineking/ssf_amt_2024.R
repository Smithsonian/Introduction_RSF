#' # Step selection functions
#'  Björn Reineking, 2024-06-19
#'
#' ## Step selection functions: Literature
#' - Forester et al. (2009) Ecology 90: 3554–3565
#' - Avgar et al. (2016) Methods in Ecology and Evolution 7: 619–630
#'
#' ## Step selection functions: A receipe
#'
#'- Problem formulation
#'- Data collection
#'- Data pre-processing
#'- Modelling
#'- Interpretation
#'- Iterate process
#'
#' We will use the buffalo data set.
#'
#' ## Loading packages
#+ results='hide', message=FALSE, warning=FALSE
library(ggplot2)
library(animove)
library(survival)
library(MASS)
library(lubridate)
library(dplyr)
library(nlme)
library(pbs)
library(circular)
library(CircStats)
library(amt)
library(ctmm)
library(terra)
library(sf)

#' ## Load buffalo data
#' 
#' 
data(buffalo_utm_mv2) # Load buffalo data from animove package

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

#' ## Environmental data: topography, waterways, and NDVI
data(buffalo_env)
buffalo_env <- rast(buffalo_env) # conversion from raster to terra objects

#' ## Always inspect your data: summary statistics
summary(buffalo_tracks)

#' We can plot the buffalo tracks, but they do not have a plot method, so
#' we need to give a bit more information to the points() function.
plot(buffalo_env[[1]])
points(y_ ~ x_, data = buffalo_tracks, col = factor(buffalo_tracks$id))

#' To speed up analyses, we will only work with one individual, Cilla
cilla <- filter(buffalo_tracks, id == "Cilla")

#' ## Inspect data
hist(step_lengths(cilla))
which(step_lengths(cilla) > 5000)

#' The very first step is unusually long; let us plot the first day in red on top of the full trajectory.
plot(cilla, type = "l")
lines(filter(cilla, t_ < min(t_) + days(1)), col = "red")

#' Let us exclude the full first day
cilla <- filter(cilla, t_ > min(t_) + days(1))

#' ## Thin movement data and split to bursts
#' - We reduce the data set to observations that are within a certain time step range. The SSF assumes that velocities of successive steps are independent, so we should thin sufficiently. Here we go for 3 hours.
#' - There is some tolerance around the target time interval of 3 hours. When two observations are separated by less than the threshold, the second observation is removed.
#' - When two observations are separated by more than the upper threshold, the observations are assigned to different bursts.
#'
#' It is a good idea to perform the analysis at several temporal scales, i.e. different step durations.
#'
#' The initial sampling rate of Cilla is about 1 hour:
#'
summarize_sampling_rate(cilla)

#' ## Prior analysis of spatio-temporal autocorrelation
#' SSF assumes that the velocities are not autocorrelated. So we cannot simply do the analysis
#' at the highest temporal resolution of the data.
#' We could try to use the downweighting trick that we use for the RSF, but for the SSF,
#' the time between successive steps will also affect the point estimates of the parameters,
#' so the problem is a bit more tricky.
#' As a rule-of-thumb, I would downsample the data such that the autocorrelation in velocities
#' has decayed to something between 2% and 5%.
#' Say you had a sampling of 1 hour, and the SSF should be done at 3 hours given this
#' rule of thumb, then you could still use all data, by constructing 3 step data sets, each starting one hour apart, and
#' give each a weight of 1/3 in the likelihood.
#'
#'
cilla_telemetry <- as_telemetry(cilla)
plot(variogram(cilla_telemetry), xlim = c(0, 10 * 3600))

GUESS <- ctmm.guess(cilla_telemetry, interactive=FALSE)

if (file.exists("cilla_fit.rds")) {
  FIT <- readRDS("cilla_fit.rds")
} else {
  FIT <- ctmm.select(cilla_telemetry, GUESS) 
  saveRDS(FIT, "cilla_fit.rds")
}

plot(variogram(cilla_telemetry), xlim = c(0, 10 * 3600))
abline(v = FIT$tau["velocity"]  * -log(0.01) / 3600, col = "blue")
abline(v = FIT$tau["velocity"]  * -log(0.02) / 3600, col = "red")
abline(v = FIT$tau["velocity"]  * -log(0.05) / 3600, col = "yellow")
legend("bottomright", lty = 1, col = c("blue", "red", "yellow"), legend = c("1%", "2%", "5%"), title = "Velocity\nautocorrelation",
       bty = "n")

#' Now we resample to 3 hour intervals, with a tolerance of 15 minutes
step_duration <- 3
cilla <- track_resample(cilla, hours(step_duration), tolerance = minutes(15))

#' Look at the new sampling rate
summarize_sampling_rate(cilla)

#' So there is at least two observations that are more than 3 hours 15 minutes apart, so there should be at least two bursts:
table(cilla$burst_)

#' If there are bursts, we may want to filter bursts with very few locations. For example, to calculate a turning angle, we need at least three locations. So we often will want to filter out bursts with at least 3 observations:
cilla <- filter_min_n_burst(cilla, 3)

#' Convert locations to steps. We will have fewer rows in the step data frame than in the track data frame because the final position is not a complete step.
ssf_cilla <- steps_by_burst(cilla)

#' We still have steps without a turning angle (the first step in a burst)
which(is.na(ssf_cilla$ta_))
ssf_cilla <- filter(ssf_cilla, !is.na(ta_))

#' ## Empirical distances and turning angles
par(mfrow = c(1, 2))
hist(ssf_cilla$sl_, breaks = 20, main = "",
  xlab = "Distance (m)")
hist(ssf_cilla$ta_,  main="",breaks = seq(-pi, pi, len=11),
      xlab="Relative angle (radians)")

#' ## Fit gamma distribution to distances
fexp <- fitdistr(ssf_cilla$sl_, "exponential")
fgamma <- amt::fit_distr(ssf_cilla$sl_, "gamma")
par(mfrow = c(1, 1))
hist(ssf_cilla$sl_, breaks = 50, prob = TRUE,
     xlim = c(0, 8000), ylim = c(0, 2e-3),
     xlab = "Step length (m)", main = "")
plot(function(x) dexp(x, rate = fexp$estimate), add = TRUE, from = 0.1, to = 8000, col = "red")
plot(function(x) dgamma(x, shape = fgamma$params$shape,
                        scale = fgamma$params$scale), add = TRUE, from = 0.1, to = 8000, col = "blue")
legend("topright", col = c("red", "blue"), lty = 1,
       legend = c("exponential", "gamma"), bty = "n")

#' ## Fit von Mises distribution to angles
fvmises <- fit_distr(ssf_cilla$ta_, "vonmises")
par(mfrow = c(1, 1))
hist(ssf_cilla$ta_, breaks = 50, prob = TRUE,
     xlim = c(-pi, pi),
     xlab = "Turning angles (rad)", main = "")
plot(function(x) dvonmises(x, mu = 0, kappa = fvmises$params$kappa), add = TRUE, from = -pi, to = pi, col = "red")


#' Create random steps. We typically get a warning that "Step-lengths or turning angles contained NA, which were removed", because of the missing turning angles at the start of a burst.
set.seed(2)
ssf_cilla <- steps_by_burst(cilla)
ssf_cilla <- random_steps(ssf_cilla, n_control = 200)

#' ## Sanity check: plot the choice set for a given step
my_step_id <- 4
ggplot(data = filter(ssf_cilla, step_id_ == my_step_id | (step_id_ %in% c(my_step_id - 1, my_step_id - 2) & case_ == 1)),
       aes(x = x2_, y = y2_)) + geom_point(aes(color = factor(step_id_))) + geom_point(data = filter(ssf_cilla, step_id_ %in% c(my_step_id, my_step_id - 1, my_step_id - 2) & case_ == 1), aes(x = x2_, y = y2_, color = factor(step_id_), size = 2))

#' ## Extract environmental covariates
#' I recommend to always use the option "both", which provides the environmental conditions at the start and the end of the step.
#' The condition at the end are what we use for selection, and the conditions at the start can be used
#' to modify e.g. turning angles and step length.
ssf_cilla <- amt::extract_covariates(ssf_cilla, buffalo_env, where = "both")

#' ## Add variable hour
#' Adding hour modelling diurnal variation in step lengths, turning angles, and preference for environmental conditions
ssf_cilla <- mutate(ssf_cilla, "hour" = hour(t1_) + minute(t1_) / 60)

#' Remove NA's
ssf_cilla <- ssf_cilla[complete.cases(ssf_cilla),]

#' ## A first model
m_1 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + slope_end + elev_end +
                    mean_NDVI_end + var_NDVI_end + strata(step_id_))
summary(m_1)

#' ## Collinearity
#' In statistics, multicollinearity (also collinearity) is a phenomenon in which one predictor variables in a multiple regression model can be linearly predicted from the others with a substantial degree of accuracy. In this situation the coefficient estimates of the multiple regression may change erratically in response to small changes in the model or the data." [Wikipedia, accessed 29.08.2017](https://en.wikipedia.org/wiki/Multicollinearity)
#'
#' One way of dealing with collinearity is to select a subset of variables that is sufficiently uncorrelated [Dormann et al. 2013](http://onlinelibrary.wiley.com/doi/10.1111/j.1600-0587.2012.07348.x/abstract). Here we simply look at pairwise correlation between predictors.
#'
round(cor(ssf_cilla[, c("slope_end", "elev_end",
  "water_dist_end", "mean_NDVI_end", "var_NDVI_end")]), 2)

#' elev and water_dist are positively correlated > 0.7
#'
#' ## Which collinear variable to pick?
#' - The one that is more relevant
#' - The one that is by itself a better predictor
#'

m1_water <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + water_dist_end + strata(step_id_))
m1_elev <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + strata(step_id_))
AIC(m1_water$model)
AIC(m1_elev$model)

#' So we pick elev, because it by itself explains the movement better.
#'
#' Fit step selection function
m_1 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + slope_end + elev_end +
                    mean_NDVI_end + var_NDVI_end + strata(step_id_))
summary(m_1)

#' slope and var_NDVI do not contribute significantly to the fit.
#'
#' ## Model selection
#' Model selection is a vast topic. I recommend using only few models with ecological justification, rather than
#' searching for the "best" model in a huge model space.
#' Here we just use stepwise backward selection based on AIC

m_2 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + slope_end + elev_end + mean_NDVI_end + strata(step_id_))
AIC(m_1$model)
AIC(m_2$model)
summary(m_2)

m_3 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + mean_NDVI_end + strata(step_id_))
AIC(m_3$model)
summary(m_3)

#' ## Model checking: serial autocorrelation
#' Forester et al. 2009 Ecology 90:3554–3565.
#' Calculate  deviance residuals for each stratum (i.e., the sum of the residuals for the case and all associated controls).
#'
ssf_residuals <- function(m, data) {
  df <- tibble::as_tibble(data.frame("time" = data$t1_,
                                       "residuals" = residuals(m$model, type = "deviance")))
  df <- df |> dplyr::group_by(time) |> dplyr::summarise(residuals = sum(residuals))
  df$group <- 1
  df
}

resid_df <- ssf_residuals(m_3, ssf_cilla)

#' Fit an intercept-only mixed-effects model using lme() from the nlme package.
#'
rm1 <- lme(residuals ~ 1, random = ~ 1 | group,
           data = resid_df)
plot(ACF(rm1, maxLag = 40), alpha = 0.05)

#' So we see that there is some significant autocorrelation at lags up to 8.
#'
#' One effect of residual temporal autocorrelation is too extreme p-values, but it may also cause bias in parameter estimates.
#' Forester et al. 2009 suggest a way to estimate more appropriate confidence intervals and p-values.
#'
#' ## Model evaluation
#' - R2 is low. Always is for this type of model.
#' - Not yet clear to me what a good performance index would be. Perhaps https://github.com/aaarchmiller/uhcplots will help.
#' - Cross-validation
#'     - split steps in e.g. 5 folds (long stretches better - should be long enough so that autocorrelation in residuals has tapered off)
#'     - leave each fold out, refit and predict to left-out fold
#'
#' ## Interpretation
#'  - Map preference
#' - Response functions
#'
#' ## Map habitat preference ("Habitat map")
#' Caveat: This habitat preference in general will not match the range distribution of the animal (i.e. how much time it spends where).
#' For the UD, a generic approach is to do simulations from the fitted model (see below).
#' See for more information [Signer et al. (2024)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14263) Simulating animal space use from fitted integrated Step-Selection Functions (iSSF)
#' 
#' The habitat_preference function assumes that all environmental layers are
#' represented in one raster stack. The habitat_preference function does not work for the "barrier" model.
#'

habitat_preference <- function(model, object, include_avail = TRUE, data_crs = crs(object)) {
  remove_end <- function(x) {
    names(x) <- gsub("_end", "", names(x))
    x
  }
  beta <- coef(model)
  exclude_terms <- c("strata","ta_", "sl_")
  for(i in exclude_terms) {
    if (length(grep(i, names(beta)))) {
      beta <- beta[-grep(i, names(beta))]
    }
  }
  beta <- remove_end(beta)
  if(include_avail) {
    xy <- as.data.frame(xyFromCell(object, 1:ncell(object)))
    xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = crs(object))
    xy <- st_transform(xy, data_crs)
    xy <- st_coordinates(xy)
    colnames(xy) <- c("x2_", "y2_")
  } else {
    xy <- cbind("x2_" = rep(0, ncell(object)), "y2_" = rep(0, ncell(object)))
  }
  newdata <- as.data.frame(cbind(xy, values(object)))
  f <- as.formula(paste("~", paste(names(beta), collapse = " + ")))
  attr(newdata, "na.action") <- "na.pass"
  X <- stats::model.matrix.default(f, data = newdata, na.action = na.pass)
  lambda <- X[,-1] %*% beta
  r <- object[[1]]
  r[] <- exp(lambda - max(lambda, na.rm = TRUE))
  r <- r / sum(values(r), na.rm = TRUE)
  r
}

habitat_map <- habitat_preference(m_1, buffalo_env)

plot(habitat_map)
lines(ssf_cilla[, c("x1_", "y1_")])

#' Now zoom in
plot(crop(habitat_map, ext(as_sf(cilla)) + 5000))
lines(ssf_cilla[, c("x1_", "y1_")])

#' The model is strongly driven by elevation
plot(crop(buffalo_env[["elev"]], ext(as_sf(cilla)) + 5000))
lines(ssf_cilla[, c("x1_", "y1_")])

#' ## Iterate
#' Here: a model with time-varying preference for mean_NDVI
#' We group observation in 3 hour bins to smooth the picture
#'
boxplot(mean_NDVI_end ~ I(floor(hour/3)*3), data =
  filter(ssf_cilla, case_ == 1),xlab = "Time of day",
  ylab = "mean NDVI")

#' We can do the same for other variables, including those of the "movement kernel", e.g. distance
boxplot(sl_ ~ floor(hour), data =
          filter(ssf_cilla, case_ == 1), xlab = "Time of day",
        ylab = "dist")

#' What behavioural rhythm do these figures suggest?
m_time_ndvi <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + mean_NDVI_end +
                   mean_NDVI_end:pbs(hour, df = 5, Boundary.knots = c(0,24)) + strata(step_id_))
m_time_dist <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + mean_NDVI_end +
                   sl_:pbs(hour, df = 5, Boundary.knots = c(0,24)) + strata(step_id_))

#' ## Predictions with the model: response function
#'
extract_stratum <- function(object) {
  attr(object$terms, 'special')$strata[1]
}
stratum <- extract_stratum(m_time_ndvi$model)
pred_data_ndvi <- data.frame("step_id_" = stratum, ta_ = 0, sl_ = 1, elev_end = 0, mean_NDVI_end = 1, hour = seq(0, 24, len = 101))
m_time_ndvi_clogit <- clogit(case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + mean_NDVI_end +
                            mean_NDVI_end:pbs(hour, df = 5, Boundary.knots = c(0,24)) + strata(step_id_), data = ssf_cilla)
pred_time <- survival:::predict.coxph(m_time_ndvi_clogit, newdata = pred_data_ndvi, se.fit = TRUE)
upper <- pred_time$fit + 1.96 * pred_time$se.fit
lower <- pred_time$fit - 1.96 * pred_time$se.fit

par(mfrow = c(1, 1))
plot(pred_data_ndvi$hour, pred_time$fit, type = "l",
  ylim = range(c(upper, lower)), xlab = "Time of day",
  ylab = "Preference mean_NDVI")
lines(pred_data_ndvi$hour, upper, lty = 2)
lines(pred_data_ndvi$hour, lower, lty = 2)
abline(h = 0, lty = 3)


#' ## Parameter estimates for the movement model
#' When using importance sampling in SSF (the default in amt), 
#' the parameter estimates for the movement kernel (e.g. step length, 
#' log(step length), cos(turn angle)) represent differences to the parameter 
#' values of the proposal distribution. 
#' Parameter values for the step length and turn angle distributions have to be 
#' calculated from the fitted ssf model and the proposal distribution. 
#' See [Avgar et al. (2016), Appendix 3](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12528) for more information.

# Proposal distribution - step length
# sl_distr(m_time_ndvi_fit)
initial_sl_distr <- sl_distr(ssf_cilla)
updated_sl_distr <- update_sl_distr(m_time_ndvi, beta_log_sl = "log(sl_)")

plot_sl_distr <- function(x, n = 1000, upper_quantile = 0.99, plot = TRUE) {
  # based on amt:::plot_sl_base
  stopifnot(x$name == "gamma")
  xx <- x$params
  to <- qgamma(upper_quantile, shape = xx$shape, scale = xx$scale)
  xs <- seq(0, to, length.out = n)
  ys <- dgamma(xs, shape = xx$shape, scale = xx$scale)
  if (plot) {
    plot(xs, ys, type = "l", ylab = "Probability", xlab = "Distance")
  }
  invisible(data.frame(sl = xs, d = ys))
}

initial_df <- plot_sl_distr(initial_sl_distr)
updated_df <- plot_sl_distr(updated_sl_distr)

plot(initial_df, type = "l", xlab = "Distance", ylab = "Density")
lines(updated_df, col = "red")
#' Proposal distribution - turning angle
initial_ta_distr <- ta_distr(ssf_cilla)
updated_ta_distr <- update_ta_distr(m_time_ndvi, beta_cos_ta = "cos(ta_)")
initial_ta_distr$params$kappa
updated_ta_distr$params$kappa

ta_values <- seq(-pi, pi, len = 101)
plot(ta_values, dvonmises(ta_values, 0, initial_ta_distr$params$kappa), type = "l", xlab = "Turning angle (radians)",
     ylab = "Density")
lines(ta_values, dvonmises(ta_values, 0, updated_ta_distr$params$kappa), col = "red")

#' The calculation of movement parameters becomes more interesting (challenging) 
#' when we have interactions between 
#' the movement model parameters (e.g. with environmental conditions).
#' See ?update_gamma or ?update_vonmises
#' 

#' ## Simulating with the model
set.seed(23)
k1 <- amt::redistribution_kernel(m_3, map = buffalo_env, 
                                 start = amt:::make_start.steps_xyt(ssf_cilla[1, ]),
                                  landscape = "continuous", tolerance.outside = 0.2, 
                                 as.rast = FALSE,
                                  n.control = 1e3)
s1 <- amt::simulate_path(k1, n.steps = 500)

extent_tracks <- function(x, y) {
  df <- data.frame(na.omit(rbind(x[,c("x_", "y_")], y[,c("x_", "y_")])))
  terra::ext(c(range(df$x_), range(df$y_)))
}

elev_crop <- crop(buffalo_env[["elev"]], extent_tracks(s1, cilla) + 2000)
plot(elev_crop)
lines(cilla)
lines(s1$x_, s1$y_, col = "red")

#' # Home ranging behaviour (Thanks, Chris!)
m_hr <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) +
                    elev_end + water_dist_end + x2_ + y2_ + I(x2_^2 + y2_^2) + strata(step_id_))
summary(m_hr)

k2 <- amt::redistribution_kernel(m_hr, map = buffalo_env, start = amt:::make_start.steps_xyt(ssf_cilla[1, ]),
                                 landscape = "continuous", tolerance.outside = 0.01, as.rast = FALSE,
                            n.control = 1e3)
set.seed(2)
s2 <- amt::simulate_path(k2, n.steps = 1000)

plot(crop(buffalo_env[["elev"]], extent_tracks(s2, cilla) + 2000))
lines(cilla)
lines(s2$x_, s2$y_, col = "red")

#' Plot simulated track on the habitat preference map
habitat_map_hr <- habitat_preference(m_hr, buffalo_env)
terra::plot(crop(habitat_map_hr, extent_tracks(s2, cilla) + 2000))
lines(cilla)
lines(s2$x_, s2$y_, col = "red")


#' ## Barriers
water <- buffalo_env[["water_dist"]] > 100
water <- crop(water, ext(water) - 5000)
plot(water)
ww <- patches(water, zeroAsNA = TRUE)
plot(ww)
ww[is.na(ww)] <- 2
names(ww) <- "water_crossed"
buffalo_env_2 <- c(ww, crop(buffalo_env, ww))

set.seed(2)
ssf_cilla <- cilla |> steps_by_burst() |> random_steps(n_control = 1000) |>
  extract_covariates(buffalo_env_2, where = "both") |>
  mutate(hour = hour(t1_) + minute(t1_) / 60) |>
  na.omit()


m_crossing <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) +
                    elev_end + water_dist_end + x2_ + y2_ + I(x2_^2 + y2_^2) +
                    I(water_crossed_end != water_crossed_start) + strata(step_id_))

summary(m_crossing)

set.seed(2)
k3 <-  amt::redistribution_kernel(m_crossing, map = buffalo_env_2,
                                  start = amt:::make_start.steps_xyt(ssf_cilla[1, ]),
                                  landscape = "continuous", tolerance.outside = 0.01,
                                  as.rast = FALSE, n.control = 1e3)

s3 <- amt::simulate_path(k3, n.steps = 1000)

plot(crop(buffalo_env[["elev"]], extent_tracks(s3, cilla) + 2000))
lines(cilla)
lines(s3$x_, s3$y_, col = "red")


#' ## From preference to range distribution maps
#' There is a fast method if we have a symmetric and temporally stable jump kernel, e.g. exponential, and no effect of step angles:
#' Barnett, A. & Moorcroft, P. (2008) Analytic steady-state space use patterns and rapid computations in mechanistic home range analysis. [Journal of Mathematical Biology, 57, 139–159](https://link.springer.com/article/10.1007/s00285-007-0149-8).
#' The generic but computationally expensive method is to do simulations: Signer et al. (2023) Simulating animal space use from fitted integrated Step-Selection Functions (iSSF). [Methods Ecol Evol 15: 43-50](https://doi.org/10.1111/2041-210X.14263)
#'
do_steady_state_UD <- TRUE
if(do_steady_state_UD) {
  set.seed(23)
  s_looong <- amt::simulate_path(k2, n.steps = 1e4) # We should use many more steps here.
  xx <- amt::make_track(s_looong, x_, y_, t_, crs = crs(buffalo_env))
  xx <- amt::as_telemetry(xx)
  ctmm::projection(xx) <- crs(buffalo_env, proj = TRUE)
  xx_guess <- ctmm.guess(xx, ctmm(isotropic = TRUE, range = TRUE), interactive = FALSE)
  xx_fit <- ctmm.fit(xx, CTMM = xx_guess)
  xx_akde <- akde(xx, xx_fit, grid = buffalo_env[[1]])
  xx_r <- ctmm::raster(xx_akde, DF = "PMF")
  plot(xx_r, col = rev(gray.colors(51, start = 0, end = 1, alpha = NULL)))
}

#' # Dependence of results on step interval
#'
step_durations <- 3:12
do_run <- TRUE
buffalo_ids <- levels(factor(buffalo_tracks$id))
if(do_run) {
step_interval_simulation_amt <- lapply(buffalo_ids, function(animal) {
  lapply(step_durations, function(step_duration) {
    ssf_animal <- filter(buffalo_tracks, id == animal) |>
      track_resample(hours(step_duration), tolerance = minutes(15)) |>
      filter_min_n_burst(3) |>
      steps_by_burst() |>
      random_steps(n = 200) |> 
      extract_covariates(buffalo_env) |>
      na.omit() |>
      filter(sl_ > 0)
    m_1 <- clogit(case_ ~ cos(ta_) + sl_ + log(sl_) + elev +
                    mean_NDVI + strata(step_id_), data = ssf_animal)
    list("coef" = coef(m_1), "confint" = confint(m_1))
  })
})
}

#' ## Parameter estimates
do_run <- TRUE
if (do_run)  {
for (j in seq(step_interval_simulation_amt)) {
  model_list <- step_interval_simulation_amt[[j]]
  name <- names(split(buffalo_utm))[j]
  coefs <- sapply(model_list, function(x) x$coef)
  ci_lower <- sapply(model_list, function(x) x$confint[, 1])
  ci_upper <- sapply(model_list, function(x) x$confint[, 2])
  par(mfrow = c(3, 2), mar=c(4,5,1,1), oma = c(1,1,5,1))
  for (i in rownames(coefs)) {
    plot(c(0,0), xlim = range(step_durations),
         ylim = range(c(0, ci_lower[i, ],
                        coefs[i,],
                        ci_upper[i, ])), type = "n", xlab = "Step duration (hr)",
         ylab = i)
    abline(h = 0, lty = 2, col = "red")
    lines(step_durations, ci_lower[i, ], lty = 3)
    lines(step_durations, ci_upper[i, ], lty = 3)
    lines(step_durations, coefs[i, ])
  }
  mtext(name, outer = TRUE)
}
}

#' ## Some of the stuff to expand on
#' - Mixed effects models -> see ssf_multiple_animals
#' - Interaction with changes of internal states (coming from another model/other observations)
#'

