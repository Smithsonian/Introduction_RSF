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
#' 
#' Loading packages
#+  results='hide', message=FALSE, warning=FALSE
library(animove)
library(ctmm)
library(future.apply) # for parallel processing

plan("multisession", workers = min(6, parallelly::availableCores()))
#' ## Load buffalo data
data("buffalo")

#' ## Environmental data: topography, waterways, and NDVI
data(buffalo_env)

#' Project buffalo data to projection of rasters (needed later when
# doing akde with a provided grid
ctmm::projection(buffalo) <- raster::crs(buffalo_env)

#' Create named list of rasters
buffalo_env_list <- list("elev" = raster::raster(buffalo_env, "elev"),
           "slope" = raster::raster(buffalo_env, "slope"),
           "var_NDVI" = raster::raster(buffalo_env, "var_NDVI"))


#' Create a list of ctmm models
buffalos_guess <- lapply(buffalo, ctmm.guess, 
                         CTMM = ctmm(isotropic = TRUE), 
                         interactive = FALSE)
buffalos_ctmm <- future_mapply(ctmm.select, buffalo, buffalos_guess, 
                               SIMPLIFY = FALSE, future.seed = TRUE)

#' Calling akde on the list of buffalos and their ctmms ensures 
#' that a common grid is used; here we explicitly provide a common grid
buffalos_akde <- akde(buffalo, buffalos_ctmm, grid = buffalo_env_list[[1]])
#' See ctmm_3_meta.R from Inês and Chris for more information

#' Calculate rsf.fit for each individual
buffalos_rsf <- future_mapply(rsf.fit, buffalo, buffalos_akde, 
                       MoreArgs = list(R = buffalo_env_list, 
                                       integrator = "Riemann"), 
                       SIMPLIFY = FALSE, future.seed = TRUE)

#' Have a look at the individual rsf selection estimates
selection_coef <- sapply(buffalos_rsf, "[[", 3)
colnames(selection_coef) <- names(buffalo)
selection_coef
#' Except for Queen, buffalos have negative estimates for elevation 
#' (not necessarily statistically significant)

#' Calculation of population-level rsf
mean_rsf <- mean(buffalos_rsf)
summary(mean_rsf)

#' Map the suitability for the "average" animal
suitability_mean <- ctmm::suitability(mean_rsf, 
                                      R = buffalo_env_list, 
                                      grid = buffalo_env_list[[1]])
raster::plot(suitability_mean)
