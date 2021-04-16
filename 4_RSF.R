## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----RSF Libraries, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------
# Remove everything from memory
rm(list=ls())

set.seed(533)

# Make sure required libraries are loaded
library(amt) 
library(dplyr)
library(ggplot2)
library(jtools)
library(lme4)
library(lubridate)
library(raster)
library(sp)
library(sjPlot)
library(usdm)
library(visreg)

# TimeZone
Timezone1 <- "UTC"
Timezone2 <- "Africa/Nairobi" 

# Spatial reference information
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326
#https://spatialreference.org/ref/esri/102022/
AEA.Africa.proj <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"   #ESRI:102022


## ----RSF Load--------------------------------------------------------------------------------------------------------------------------------
# Read GPS collar data saved as an object in previous exercise
# ******************************
# ******************************
load("./Data/wild.Rdata")

# Read spatial data - Rasters for Athi-Kaputiei Plains
# ******************************
# ******************************
load("./Data/wild_SpatialData.Rdata")

# Stack the rasters together
rsf.stack <- stack(anth_risk,fence_dist,prirds_dist,secrds_dist,river_dist,waterpts_dist,woody_dist)

# All files have been project to Albers Equal Area, with a spatial resolution of 250m
# I have included:
# anth_risk - Anthropogenic Risk, simply an index of human footprint made specifically for this ecosystem (expected negative response)
# Fence_dist - Fence Distance, with fences manually map by a Kenyan field team (expected negative response)
# prirds_dist - Primary Road Distance, the distance from primary/paved roads (expected negative response)
# secrds_dist - Secondary Road Distance, the distance from secondary/unpaved roads (expected null response)
# river_dist - River Distance, the distance from permanent rivers (depending on season, but wildebeest must drink - greater attraction in dry season, but a predation risk)
# waterpts_dist - Water Point Distance, the distance to mapped water wells (expected positive response, stronger in dry season)
# woody_dist -Woody Distance, the distance to woody vegetation (based on VCF) (expected negative response - predation risk)

# Project to Albers Equal Area (maintains distances across study areas and gets around the problem of working across UTM zones)
# DEFINE the projection first
wild.Athi <- SpatialPointsDataFrame(coords = wild.Athi[,c("location.long","location.lat")],
                                    data = wild.Athi,
                                    proj4string = CRS(LatLong.proj)) 
# Simple plot, to make sure things look correct
#plot(wild.Athi, xlab="Longitude",ylab="Latitude", main = "Athi-Kaputiei Plains Wildebeest Data", axes=T) # Don't need to specify the coordinates because R recognizes the file as spatial.  If the file was a dataframe, we would need to specify the X and Y fields

# Transform to Albers Equal Area Project - Good practice to see a different projection and to see that the code is essential the same.  This projection is good for maintaining areas across regions that span multiple UTM zones.  The parameters from the projection are specific to Africa.  
wild.Athi <- spTransform(wild.Athi,
                         CRS =  CRS(AEA.Africa.proj))
#plot(wild.Athi, xlab="Easting",ylab="Northing", main = "Athi-Kaputiei Plains Wildebeest Data", axes=T)

# Plot one of the images in the raster stack with one of the tracking datasets
plot(rsf.stack[[1]])
points(wild.Athi[wild.Athi$id == 30077,], pch = '.', cex = 2, col='black')

# Convert to dataframe for statistical analysis
wild.Athi <- as.data.frame(wild.Athi)

# Let's clean up this dataset a bit, keeping only the columns important for this analysis
# wild.Athi <- wild.Athi %>% select(
#   id = id,
#   animal = individual.local.identifier,
#   sex = animal.sex,
#   timestamp = timestamp,
#   easting = location.long.1,
#   northing = location.lat.1
# )

# Select function was giving a hard time...re-arranging/re-selecting as follows
wild.Athi <- wild.Athi[,c(1,11,13,2,15:16)]
colnames(wild.Athi) <- c("id","animal","sex","timestamp","easting","northing")

# In this case, we will not be incorporating time (like you would in a SSF), but could be a vastly important component depending on your research question
# Good to check our formats and make sure the timestamps are in the correct format.
str(wild.Athi)


## ----RSF Fitting-----------------------------------------------------------------------------------------------------------------------------
# Select a single individual (e.g., 30077)
wb <- wild.Athi[wild.Athi$id==30077,]
#plot(wb$easting,wb$northing, xlab="Easting", ylab="Northing", main=paste0("Animal: ",unique(wb$id)," - Ntishya"), cex=2, pch=".",col="blue",frame=F)

# We've already checked complete cases, already checked timestamp as as.POSIXct, and removed duplicates
# We need to create a track object using mk_track function (this is similar to as.ltraj in adehabitatLT)
# Note, for this analysis, we actually will not be incorporated time (that would be a SSF)
wb <- mk_track(wb, .x = easting, .y = northing,
                .t = timestamp, crs = CRS(AEA.Africa.proj), order_by_ts = T,
                id = id, animal = animal, sex = sex)

# Add steplengths to object
# Simply calculate the distance moved between steps and adding a new column, called sl (steplength)
# Direction/turning angles can also be easily added using `direction_abs()` command.  Not adding here because not incorporating movement process
wb <- wb %>% mutate(
  sl = step_lengths(.))
summary(wb$sl)

# Summarize the distribution of time intervals between successive locations to get a general impression for the sampling rate.
summarize_sampling_rate(wb)

# Let's sample to 3 hours, but could do others that are more appropriate (daily?)
# Important to be able to defend 'why' you are doing something
wb <- wb %>% track_resample(rate = hours(3), tolerance = minutes(20))
wb
summary(wb$sl) # Not that the movement rate now changes....because we are simply connecting the points...so more representative of our fix collection schedule than actual movement (need continuous time model for this)

# Ad-hoc thinning of Use locations
# we will randomly sample 20% of the locations to use in modeling.
# Here, I'm not using any spatial statistical tests (acf, Moran's I, Geary's C) to determine if this thinning amount is enough
rcd.amt <- ceiling(nrow(wb)*0.20)
wb <- wb[sample(1:nrow(wb),size = rcd.amt),c("x_","y_","t_")]
#wb # Not that object is a tibble (a dataframe of dataframes)
#plot(wb$easting,wb$northing, xlab="Easting", ylab="Northing", main=paste0("Animal: ",unique(wb$id)," - Ntishya"), cex=2, pch=".",col="blue",frame=F)


## ----Availability----------------------------------------------------------------------------------------------------------------------------
# Lot's of ways to assess availability
# Here we are generating 5 times the number of Use points
wb_RSF<-random_points(wb,n=nrow(wb)*5,typ="random") # Note that this crease a 'case' field (True/False: 1/0)
#unique(wb_RSF$case_)
plot(wb_RSF) # Plot shows Use/Availability

# Could also create points in a "regular" pattern....arguably, could be better
#wb2 <- random_points(wb, n=nrow(wb.3hr)*5, typ="regular") 
#plot(wb2)

# Extract all the raster variables at all Use/Available points.
# Nice that this is done all with a single command....and very quickly.  Amazing!
# NOTE: we are not considering time here - see extract_covariates_var_time() function to do so, a bit more complicated
wb_RSF <- wb_RSF %>% extract_covariates(rsf.stack)

# Compare values between Use/Available
round(aggregate(wb_RSF[,4:10], by=list(as.numeric(wb_RSF$case_)), FUN=mean),digits=2)


## ----Collinearity----------------------------------------------------------------------------------------------------------------------------
# Assess collinearity
test.corr <- as.data.frame(wb_RSF)
vifstep(test.corr[,4:10]) # VIF indicates no collinearity problem, but there is a high level of correlation (max: 0.84)
cor(test.corr[,4:10])

# How to address?  Should they be removed or are we less concerned?
# Let's fit different models and see how the standard error of our included variables changes between fitted models
# Notice I have scaled the variables
MA <- glm(case_ ~ scale(anth_risk) + scale(fence_dist) + scale(waterpts_dist) + scale(prirds_dist) + scale(river_dist) + scale(secrds_dist) + scale(woody_dist), family = binomial(link="logit"), 
     data =wb_RSF)
summary(MA)

# Remove waterpts_dist
MB <- glm(case_ ~ scale(anth_risk) + scale(fence_dist) + scale(prirds_dist) + scale(river_dist) + scale(secrds_dist) + scale(woody_dist), family = binomial(link="logit"), 
     data =wb_RSF)
summary(MB)

# Remove fence_dist
MC <- glm(case_ ~ scale(anth_risk) + scale(waterpts_dist) + scale(prirds_dist) + scale(river_dist) + scale(secrds_dist) + scale(woody_dist), family = binomial(link="logit"), 
     data =wb_RSF)
summary(MC)

# Look at coefficients
tab_model(MA,MB,MC, show.se=T, transform=NULL)

# Out of curiosity, how to models compare 
AIC(MA,MB,MC)


## ----Sensitivity-----------------------------------------------------------------------------------------------------------------------------
# Exploring sensitivity of RSF coefficients to the number of available points

# Setup number of available locations to sample
n.frac <- c(1, 5, 20, 50, 100) 

# Total available locations to be generated, based on the number of Use points
n.pts <- nrow(wb) * n.frac

# Number of replicates of each available location (1,5,20,50,100)
# This is number of times each model will be fit, allows to calculate variability
# For a publication, I would increase the n.rep to 100
n.rep <- 20

# Create a table which can hold the n.reps
# Extract covariates from rsf.stack
# Rename the variables
# Fit glm model for each available location (1,5,20,50,100) and for each replication (re-sampling)

# Commenting out because the process is time consuming
# *******************************************
# *******************************************
# wb.sim <- tibble(
#   n.pts = rep(n.pts, n.rep), 
#   frac = rep(n.frac, n.rep), 
#   result = map(
#     n.pts, ~
#       wb %>% random_points(n = .x) %>%
#       extract_covariates(rsf.stack) %>%
#       mutate(anth_risk= scale(anth_risk),
#              woody_dist = scale(woody_dist),
#              fence_dist = scale(fence_dist),
#              prirds_dist = scale(prirds_dist),
#              river_dist = scale(river_dist),
#              secrds_dist = scale(secrds_dist),
#              waterpts_dist = scale(waterpts_dist)) %>%
#       glm(case_ ~ anth_risk + woody_dist + fence_dist + prirds_dist + river_dist + secrds_dist + waterpts_dist,
#           data = ., family = binomial(link = "logit")) %>%
#       tidy()))

wb.sim # We see the structure has 100 rows (n.frac (5) x n.rep (20) = 100) # File loaded at the beginning of script
#str(wb.sim)
wb.sim$result[[1]]

# We now "unnest" the results to plot the coefficient estimates from individual model fits with different sets of available data (1,5,20,100) 
wb.sim %>% unnest(cols = result) %>% 
  mutate(wb.sim = recode(term, "(Intercept)"="Intercept",
                         anth_risk = "Anthropogenic disturbance", 
                         woody_dist = "Woody Distance",
                         fence_dist = "Fence Distance", 
                         prirds_dist = "Pri Rds Distance",
                         river_dist = "River Distance",
                         secrds_dist = "Sec Rds Distance",
                         waterpts_dist = "Water Point Distance")) %>% 
  ggplot(aes(factor(frac), y = estimate)) +
  geom_boxplot() + facet_wrap(~ wb.sim, scale  ="free") +
  geom_jitter(alpha = 0.2) + 
  labs(x = "Number of available points (multiplier of no. of used locations)", 
       y = "Estimate") +
  theme_light()


## ----RSF 50----------------------------------------------------------------------------------------------------------------------------------
# Generating 50 available locations per each used location within animal homerange
wb_RSF50 <- random_points(wb,n=nrow(wb)*50,typ="random")

# Extracting the covariates
wb_RSF50 <- wb_RSF50 %>% extract_covariates(rsf.stack)

# Fitting model
M2 <- glm(case_ ~ scale(anth_risk) + scale(fence_dist) + scale(waterpts_dist) + scale(prirds_dist) + scale(river_dist) + scale(secrds_dist) + scale(woody_dist), 
         family = binomial(link="logit"), 
         data =wb_RSF50)
# Get Summary & confidence intervals
summary(M2)
confint(M2)

# Plot coefficients (an easy alternative to coefplot - THANKS Pawel!!)
plot_model(M2,transform = NULL) 

# There are multiple ways to graph the response curves.  In addition to visreg, we could also plot using jtools or sjplot (many others)
visreg(M2,"anth_risk",
       scale="response", 
       ylab="Relative Selective Strength",
       xlab="Anthropogenic risk",
       partial=F,
       rug=F,
       line=list(col="black"), 
       fill=list(col="light gray"))

# Let's loop over all the variables
par(mfrow=c(2,4))
var1 <- c("anth_risk","fence_dist","waterpts_dist","prirds_dist","river_dist","secrds_dist","woody_dist")
var2 <- c("Anthropogenic Risk","Fence Distance","Water Point Distance","Primary Roads Distance","River Distance","Secondary Roads Distance","Woody Distance")

par(mfrow=c(2,4))
for (i in 1:length(var1)){
visreg(M2,var1[i],
       scale="response", 
       ylab="Probability of Selection",
       xlab=var2[i],
       partial=F,
       rug=F,
       line=list(col="black"), 
       fill=list(col="light gray"))
}
par(mfrow=c(1,1))


## ----Prediction------------------------------------------------------------------------------------------------------------------------------
# Extract the coefficient of each predictor from the model summary
coeff <- M2$coefficients

# Then, use the logistic equation to generate predictions (no B_0)
# w*(x)=exp(β1x1 + β2x2 +.... + βnxn)/exp(1 + β1x1 + β2x2 +.... + βnxn)

# where w*(x) is the relative probability of selection, dependent upon covariates X1 through Xn, and their estimated regression coefficients β1 to βn, respectively.

# Remember that we must scale the raster layers too:
anth.scale <- (anth_risk - mean(wb_RSF50$anth_risk)) / sd(wb_RSF50$anth_risk)
fence.scale <- (fence_dist - mean(wb_RSF50$fence_dist)) / sd(wb_RSF50$fence_dist)
water.scale <- (prirds_dist - mean(wb_RSF50$prirds_dist)) / sd(wb_RSF50$prirds_dist)
prirds.scale <- (secrds_dist - mean(wb_RSF50$secrds_dist)) / sd(wb_RSF50$secrds_dist)
river.scale <- (river_dist - mean(wb_RSF50$river_dist)) / sd(wb_RSF50$river_dist)
secrds.scale <- (waterpts_dist - mean(wb_RSF50$waterpts_dist)) / sd(wb_RSF50$waterpts_dist) 
woody.scale <- (woody_dist - mean(wb_RSF50$woody_dist)) / sd(wb_RSF50$woody_dist) 

# Prediction
pred <- exp(anth.scale*coeff[[2]]+ fence.scale*coeff[[3]] + water.scale*coeff[[4]] + prirds.scale*coeff[[5]] + river.scale*coeff[[6]] + secrds.scale*coeff[[7]] + woody.scale*coeff[[8]])/(1+exp(anth.scale*coeff[[2]]+ fence.scale*coeff[[3]] + water.scale*coeff[[4]] + prirds.scale*coeff[[5]] + river.scale*coeff[[6]] + secrds.scale*coeff[[7]] + woody.scale*coeff[[8]]))

# Provide Spatial Prediction - Based off of the coefficients from this single animal
plot(pred)

# We could then mask this data layer, and then use tmap an place the layer on top of a map service (like we did with the Addax)
# We loaded a spatial polygon data layer at the start of this script named "Athi_Bound"
# Create new layer and crop it to the analysis extent
pred.30077 <- mask(pred, Athi_Bound)
pred.30077 <- crop(pred.30077,y=extent(Athi_Bound))
plot(pred.30077)

# Save this file in your data directory
save(pred.30077, file="./Data/Prediction30077.Rdata")


## ----GLMER, eval=F, echo=T-------------------------------------------------------------------------------------------------------------------
## # Create object of ids
## wild.An <- unique(wild.Athi$id)
## # Creating list that contain the dataframe for each individual
## wild.Athi<-split(wild.Athi,wild.Athi$id)
## 
## # Creating track for each individual animal in the list
## wild.track <- lapply(wild.Athi, function(x) mk_track(x,
##          .x=easting,.y=northing,
##          .t=timestamp,crs = sp::CRS(AEA.Africa.proj)))
## 
## # Let's reduce the file, selecting only 20% of each use dataset
## wild.track <- lapply(wild.track, function(x) x[sample(1:nrow(x),size = ceiling(nrow(x)*0.20)),c("x_","y_","t_")])
## 
## # Generate 50 random points for each animal
## # Using a list apply here to apply the function
## wb.Athi <- lapply(wild.track, function(x) random_points(x, n=nrow(x) * 50,typ="random"))


## ----GLMER Extract, eval=F, echo=T-----------------------------------------------------------------------------------------------------------
## # Adding anthropogenic risk covariate to each dataframe by using extract_covariates function
## wb.Athi <- lapply(wb.Athi, function(x) extract_covariates(x, rsf.stack))
## 
## # Adding animal ID to each dataframe in the list.  We'll use this as the random effect.
## wb.Athi <- mapply(cbind, wb.Athi, "Animal_ID" = wild.An, SIMPLIFY=F)
## 
## # Binding all dataframes in the list into a single dataframe for analysis
## wb.Athi <- do.call(rbind,wb.Athi)


## ----GLMER Fit, eval=F, echo=T---------------------------------------------------------------------------------------------------------------
## # This model will be slow to execute
## mixed.model <- glmer(case_ ~ scale(anth_risk) + scale(fence_dist) + scale(waterpts_dist) + scale(prirds_dist) + scale(river_dist) + scale(secrds_dist) + scale(woody_dist) + (1|Animal_ID),
##                 family = binomial(link="logit"),
##                 data = wb.Athi)
## 
## # Print Model Summary
## summary(mixed.model)


## ----GLMER Response, eval=F, echo=T----------------------------------------------------------------------------------------------------------
## # Plot coefficients and response curves for inference at the population level
## plot_model(mixed.model,transform = NULL)
## 
## # Exponentiate the coefficients and plot the odds ratios
## tab_model(mixed.model)
## plot_model(mixed.model)
## 
## # Graph response curves.  Use effect plot from jtools.  Executes much quicker
## effect_plot(mixed.model,data=wb.Athi,
##             pred=anth_risk,interval = TRUE,
##             x.label = "Anthropogenic risk",
##             y.label = "Relative Selection Strength") +
##   theme_classic()


## ----GlMER Prediction, eval=F, echo=T--------------------------------------------------------------------------------------------------------
## # Extract the coefficient of each predictor from the model summary
## coef <- summary(mixed.model)
## coeff <- coef$coefficients
## 
## # Scale rasters:
## anth.scale <- (anth_risk - mean(wb.Athi$anth_risk)) / sd(wb.Athi$anth_risk)
## fence.scale <- (fence_dist - mean(wb.Athi$fence_dist)) / sd(wb.Athi$fence_dist)
## water.scale <- (prirds_dist - mean(wb.Athi$prirds_dist)) / sd(wb.Athi$prirds_dist)
## prirds.scale <- (secrds_dist - mean(wb.Athi$secrds_dist)) / sd(wb.Athi$secrds_dist)
## river.scale <- (river_dist - mean(wb.Athi$river_dist)) / sd(wb.Athi$river_dist)
## secrds.scale <- (waterpts_dist - mean(wb.Athi$waterpts_dist)) / sd(wb.Athi$waterpts_dist)
## woody.scale <- (woody_dist - mean(wb.Athi$woody_dist)) / sd(wb.Athi$woody_dist)
## 
## # Prediction
## pred.mixed <- exp(anth.scale*coeff[[2]]+ fence.scale*coeff[[3]] + water.scale*coeff[[4]] + prirds.scale*coeff[[5]] + river.scale*coeff[[6]] + secrds.scale*coeff[[7]] + woody.scale*coeff[[8]])/(1+exp(anth.scale*coeff[[2]]+ fence.scale*coeff[[3]] + water.scale*coeff[[4]] + prirds.scale*coeff[[5]] + river.scale*coeff[[6]] + secrds.scale*coeff[[7]] + woody.scale*coeff[[8]]))
## 
## # Provide Spatial Prediction - Based off of the coefficients from this single animal
## plot(pred.mixed)
## 
## # We could then mask this data layer, and then use tmap an place the layer on top of a map service (like we did with the Addax)
## # We loaded a spatial polygon data layer at the start of this script named "Athi_Bound"
## # Create new layer and crop it to the analysis extent
## pred.mixed.mask <- mask(pred.mixed, Athi_Bound)
## pred.mixed.mask <- crop(pred.mixed.mask,y=extent(Athi_Bound))
## plot(pred.mixed.mask)

