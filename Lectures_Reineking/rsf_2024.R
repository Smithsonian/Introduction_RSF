#' # Resource selection functions
#' Björn Reineking, 2024-06-20
#'
#' Model the relative density of animals (also called range distribution or utilisation distribution) as a function of environmental predictors.
#'
#' We will use the buffalo data set.
#' Install animove metapackage
#' https://github.com/AniMoveCourse/animove_R_package
#' install.packages("remotes")
#' remotes::install_github("AniMoveCourse/animove_R_package")
#'
#' This code needs ctmm version 1.1.1 (or newer)
#' remotes::install_github("ctmm-initiative/ctmm")

#'Loading packages
#+  results='hide', message=FALSE, warning=FALSE
library(animove)
library(ctmm)
library(sf)
library(mvtnorm)
library(terra)

#' ## Load environmental data: topography, waterways, and NDVI
data(buffalo_env)
# Convert from raster to terra object. Terra is the recommended package.
buffalo_env <- rast(buffalo_env) 
plot(buffalo_env)


#' ## Load buffalo data from ctmm
data(buffalo)
projection(buffalo) <- crs(buffalo_env) # project buffalo tracks to projection of raster layers

#' ## Animals and elevation
plot(buffalo_env[[1]])
buffalo_mv <- move2::mt_as_move2(buffalo) # Convert to move2 object for plotting
plot(buffalo_mv, add = TRUE)

#' To speed up analyses, we will only work with one individual, Cilla
cilla_mv <- subset(buffalo_mv, track == "Cilla") 
plot(buffalo_env[[1]], ext = extent(cilla_mv) * 2)
plot(st_geometry(cilla_mv), add = TRUE, type = "p", 
     col = "red", pch = 19, cex = 0.5)

#' # Minimal example of rsf.fit
#' Select one individual, Cilla
cilla <- buffalo$Cilla

#' Fit ctmm model
#' We start with a guesstimate of parameter values; we demand a
#' isotropic ctmm model, because this is a requirement for rsf.fit
cilla_guess <- ctmm.guess(cilla, CTMM = ctmm(isotropic = TRUE), 
                          interactive = FALSE)

# Fit ctmm model and perform model selection, using all but one core
if (file.exists("cilla_fit_rsf.rds")) {
  cilla_fit <- readRDS("cilla_fit_rsf.rds")
} else {
  cilla_fit <- ctmm.select(cilla, cilla_guess, cores = -1)
  saveRDS(cilla_fit, "cilla_fit_rsf.rds")
}

#' Create named list of rasters
#' rsf.fit can deal with raster layers that differ in resolution or even projection
buffalo_env_list <- list("elev" = raster::raster(buffalo_env[["elev"]]),
           "slope" = raster::raster(buffalo_env[["slope"]]),
           "var_NDVI" = raster::raster(buffalo_env[["var_NDVI"]]))

#' Definition of reference grid the akde predictions should get aligned to
reference_grid <- crop(buffalo_env_list[["elev"]], extent(cilla_mv) * 2)

#' Fit akde
cilla_akde <- akde(cilla, cilla_fit, grid = reference_grid)
plot(cilla_akde)
# Create a raster representation of the akde, for the reference grid
# PMF is "probability mass function". The values sum up to 1 for the whole grid,
# and represent for a given pixel the probability of finding the animal in that pixel
plot(ctmm::raster(cilla_akde, DF = "PMF"))

#' The integrator = "Riemann" option is much faster
cilla_rsf_riemann <- rsf.fit(cilla, cilla_akde, 
                             R = buffalo_env_list, integrator = "Riemann")

# Estimates and confidence intervals for all model parameters
# - Selection: The three environmental layers
# - Movement: home range area,time scales of position and velocity autocorrelation;
# derived quantities speed and diffusion
summary(cilla_rsf_riemann)

#' A suitability map - with 95% confidence interval
suitability_riemann <- ctmm::suitability(cilla_rsf_riemann, R = buffalo_env_list, 
                                         grid = reference_grid)
terra::plot(suitability_riemann)

#' Range distribution (includes the ranging behaviour)
agde_cilla <- agde(CTMM = cilla_rsf_riemann, R = buffalo_env_list, 
                   grid = reference_grid)
plot(agde_cilla)
# Plot raster of range distribution
agde_raster <- ctmm::raster(agde_cilla, DF = "PMF")
plot(agde_raster)

#' Selection-informed akde
akde_rsf <- akde(cilla, CTMM = cilla_rsf_riemann, R = buffalo_env_list, 
                 grid = reference_grid)
plot(ctmm::raster(akde_rsf, DF = "PMF"))

# Plot all three probability mass functions with a common color scale
pmf <- raster::stack(ctmm::raster(cilla_akde, DF = "PMF"),
              ctmm::raster(agde_cilla, DF = "PMF"),
              ctmm::raster(akde_rsf, DF = "PMF"))
names(pmf) <- c("akde", "rsf.fit", "akde_rsf.fit")
plot(pmf, zlim=c(0, max(raster::getValues(pmf)))) # use zlim to force same color scale

#' # Traditional RSF with downweighted Poisson regression
#' Functions to generate quadrature points and predict with the model
rsf_points <- function(x, UD, R = NULL, n = 1e5, k = 1e6, type = "Riemann",
                       rmax = 6*sqrt(UD@CTMM$sigma[1,1]),
                       interpolation = FALSE) {
  # Samples background points from a 2D normal distribution fitted to the relocation data, and extracts environmental
  # information from a raster object "R"
  # x: telemetry object
  # UD: UD object
  # R: terra object
  # n: number of background points to sample
  # k: weight of presence points
  # rmax: maximum distance for Riemann-type integration
  # interpolation: do interpolation when sampling the grid
  # When type 0´= "MonteCarlo", importance sampling is done
  stopifnot(UD@CTMM$isotropic)
  stopifnot(type %in% c("Riemann", "MonteCarlo"))
  if (type == "Riemann") {
    quadrature_pts <- values(R) #raster::getValues(R)
    xy <- as.data.frame(xyFromCell(R, 1:ncell(R)))
    xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = crs(R))
    xy <- st_transform(xy, UD@info$projection)
    xy <- st_coordinates(xy)
    r <- sqrt(((xy[,1] - UD@CTMM$mu[1]))^2 + ((xy[,2] - UD@CTMM$mu[2]))^2)
    bg <- data.frame(case_ = 0,
                     x_ = xy[r<rmax,1], y_ = xy[r<rmax,2],  w_ = prod(res(R)), k_ = k)
    bg <- cbind(bg, quadrature_pts[r<rmax,])
    bg <- sf::st_as_sf(bg, coords = c("x_", "y_"), crs = UD@info$projection)
    xx <- data.frame(case_ = 1, x_ = x$x, y_ = x$y,
                     w_ = 1/k * UD$weights * mean(UD$DOF.area), k_ = k)
    xx <- sf::st_as_sf(xx, coords = c("x_", "y_"), crs = UD@info$projection)
    xx[names(R)] <- as.data.frame(extract(R, st_transform(xx, crs(R)), 
                                          ID = FALSE,
                                          method = ifelse(interpolation, "bilinear", "simple")))
    xx <- rbind(bg, xx)
  } else {
    quadrature_pts <- MASS::mvrnorm(n, mu = UD@CTMM$mu, Sigma = UD@CTMM$sigma)
    xx <- data.frame(case_ = 0, x_ = quadrature_pts[, 1], y_ = quadrature_pts[, 2], w_ = UD@CTMM$sigma[1,1]/n, k_ = k)
    xx <- rbind(xx, data.frame(case_ = 1, x_ = x$x, y_ = x$y,
                               w_ = 1/k * UD$weights * mean(UD$DOF.area), k_ = k
    ))
    xx <- sf::st_as_sf(xx, coords = c("x_", "y_"), crs = UD@info$projection)
    xx[names(R)] <- as.data.frame(extract(R, st_transform(xx, crs(R)), 
                                          ID = FALSE,
                                          method = ifelse(interpolation, "bilinear", "simple")))
  }
  xy <- st_coordinates(xx)
  colnames(xy) <- c("x_", "y_")
  sd <- sqrt(UD@CTMM$sigma[1,1])
  xy[,1] <- (xy[,1] - UD@CTMM$mu[1])/sd
  xy[,2] <- (xy[,2] - UD@CTMM$mu[2])/sd
  xx <- cbind(xy, xx)
  xx
}

predict_rsf <- function(model, UD, object, include_avail = TRUE, data_crs = UD@info$projection) {
  if(include_avail) {
    xy <- as.data.frame(xyFromCell(object, 1:ncell(object)))
    xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = crs(object))
    xy <- st_transform(xy, data_crs)
    xy <- st_coordinates(xy)
    colnames(xy) <- c("x_", "y_")

    sd <- sqrt(UD@CTMM$sigma[1,1])
    xy[,1] <- (xy[,1] - UD@CTMM$mu[1])/sd
    xy[,2] <- (xy[,2] - UD@CTMM$mu[2])/sd

  } else {
    xy <- cbind("x_" = rep(0, ncell(object)), "y_" = rep(0, ncell(object)))
  }
  newdata <- as.data.frame(cbind(xy, values(object)))
  lambda <- as.numeric(predict(model, newdata, type = "link"))
  r <- rast(object[[1]])
  r[] <- exp(lambda - max(lambda, na.rm = TRUE))
  r <- r / sum(values(r), na.rm = TRUE)
  r
}

#' ## A minimal "classic" example
#' Generate quadrature points ("background points")
set.seed(2)
rsf_cilla_df <- rsf_points(cilla, cilla_akde, buffalo_env, interpolation = TRUE)
rsf_cilla_df <- rsf_cilla_df[!is.na(rsf_cilla_df$slope),] # Remove lines with NA values
#' Fit a downweighted Poisson regression
#' the homeranging behaviour is represented by x_ + y_ + I(-(x_^2 + y_^2)/2)
m_rsf_cilla <- glm(case_*k_ ~ x_ + y_ + I(-(x_^2 + y_^2)/2) + elev + slope + var_NDVI,
                  family = poisson(), data= rsf_cilla_df, weights = w_)
#' Summary of model and confidence intervals for parameter estimates
#+ message=FALSE, warning=FALSE
summary(m_rsf_cilla)
confint(m_rsf_cilla)

#' Map of suitability, including the home ranging behaviour
suitability_glm <- predict_rsf(m_rsf_cilla, cilla_akde, 
                               crop(buffalo_env, extent(reference_grid)))
plot(suitability_glm)

#' Map of suitability, without the home ranging behaviour
suitability_no_avail_glm <- predict_rsf(m_rsf_cilla, cilla_akde, 
                                        crop(buffalo_env, extent(reference_grid)),
                                        include_avail = FALSE)
plot(suitability_no_avail_glm)

#' Comparison of estimates from the two approaches (ctmm and "classic" via glm)
vars <- c("elev", "slope", "var_NDVI")

# Estimates are not identical but similar
cilla_rsf_riemann$beta[vars]
coef(m_rsf_cilla)[vars]

# Confidence intervals are of similar width, but location is different
#+ message=FALSE, warning=FALSE
ci_glm <- confint(m_rsf_cilla)
ci_glm[vars, ]
summary(cilla_rsf_riemann)$CI[sapply(vars, function(x) grep(x, rownames(summary(cilla_rsf_riemann)$CI))), ]

#' Create maps of suitability based on parameter estimates
M <- terra::values(buffalo_env)
eta_glm <- M[,vars] %*% coef(m_rsf_cilla)[vars]
eta_rsf <- M[,vars] %*% cilla_rsf_riemann$beta[vars]

suit_r_glm <- buffalo_env[[1]]
suit_r_glm[] <- exp(eta_glm)
suit_r_rsf <- buffalo_env[[1]]
suit_r_rsf[] <- exp(eta_rsf)

suit_raster <- c(suit_r_rsf, suit_r_glm)
names(suit_raster) <- c("ctmm", "classic")
plot(suit_raster)
plot(suit_raster, ext = extent(reference_grid))


#' ## Use vs availability
#' h_avail: Histogram of elevation values within the reference grid (here taken as available elevations)
#' h_use: Histogram of elevation values at locations used by Cilla
h_avail <- hist(values(reference_grid), prob = TRUE, col=rgb(0,0,1,1/4), main = "",
                xlab = "Elevation (m a.s.l.)", ylim = c(0,0.015))
h_use <- hist(extract(reference_grid, cilla_mv), prob = TRUE,
              breaks = h_avail$breaks, add = TRUE, col=rgb(1,0,0,1/4))
legend("topright",
       fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)),
       legend = c("available","use"),
       bty = "n"
)

# A simple rsf model containing only selection for elevation, but with
# a quadratic relationship. No homeranging is included
m_rsf_cilla_elev <- glm(case_*k_ ~ poly(elev, 2), family = poisson(), 
                        data= rsf_cilla_df, weights = w_)

# Estimate of density of available elevation. 
# We need this to calculate the integration constant to convert the 
# habitat selection values (i.e. the prediction of the habitat model)
# to the use/availability ratio.

d_avail <- density(values(reference_grid)) 
elev_pred <- data.frame("elev" = d_avail$x)
elev_pred$w <- predict(m_rsf_cilla_elev, 
                       newdata = elev_pred,
                       type = "response")

#' Calculate the integration constant
#' $$u(x) = \frac{w(x) a(x)}{\int_{E}w(X) a(X)dX} $$
#' Where u(x): use of resource level x ; a(x): availablility of resource level x
#' w(x): selection of resource level x, E: the range of resource levels.
#' K is the integration constant $\int_{E}w(X) a(X)dX$.
#' 
K <- sum(elev_pred$w * d_avail$y * diff(d_avail$x[1:2]))

plot(h_avail$mids, h_use$density/h_avail$density,
     xlab = "Elevation (m a.s.l.)",
     ylab = "use/availability ratio")
abline(h = 1, lty = 3)
lines(elev_pred$elev, elev_pred$w / K, col = "green", lwd = 2)
legend("topright", 
       col = c("black", "green"), 
       lty = c(NA,1), lwd = 2, pch = c(1, NA),
       legend = c("use/availability","rsf: poly(elev, 2)/K"),
       bty = "n"
)

#' ## Multiple animals
#'
#' - with rsf.fit: you can use the mean function on a list of rsf.fit objects
#' - "classic" approach: use glmmTMB and a mixed-effects Poisson regression: Muff, Signer & Fieberg (2020) J Anim Ecol 89: 80-92.
#'


