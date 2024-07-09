## ----setup, include=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----Libraries, eval=T, echo=T, message=F, warning=F-----------------------------------
library(amt) 
library(dplyr)
library(ggplot2)
library(jtools)
library(lme4)
library(lubridate)
library(raster)
library(sp)
library(sjPlot)
library(tmap)
library(usdm)
library(visreg)


## ----Load, eval=T, echo=T--------------------------------------------------------------
load(file = "./Data/wild.Rdata")
load(file = "./Data/Wild_SpatialData.Rdata")


## ----Examen, eval=T, echo=T------------------------------------------------------------
# Visually look at all the files
# wild Athi file is a dataframe and most be projected
# All other files are fine
head(wild.Athi)

Athi_Bound
anth_risk
fence_dist
prirds_dist
secrds_dist
river_dist
waterpts_dist
woody_dist

# Spatial reference information
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326
AEA.Africa.proj <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"   #ESRI:102022

# File has Lat/Long coordinates.  Coordinates must be defined first and then projected to Albers
wild.Athi <- SpatialPointsDataFrame(coords = wild.Athi[,c("location.long","location.lat")],
                                    data = wild.Athi,
                                    proj4string = CRS(LatLong.proj)) 

# Transform to Albers Equal Area Project
wild.Athi <- spTransform(wild.Athi,
                         CRS =  CRS(AEA.Africa.proj))


## ----Subset, eval=T, echo=T------------------------------------------------------------
# Set Seed
set.seed(99)

# Convert to dataframe for statistical analysis
wild.Athi <- as.data.frame(wild.Athi)

# Remove columns and rename
wild.Athi <- wild.Athi[,c(1,11,13,2,15:16)] # Also good to assure yourself that these columns are correct.
colnames(wild.Athi) <- c("id","animal","sex","timestamp","easting","northing")

wb2 <- wild.Athi[wild.Athi$id==2840,]


## ----Steps, eval=T, echo=T-------------------------------------------------------------
# How many observations?
nrow(wb2)

# 4777 observations

# Calculate track
wb2 <- mk_track(wb2, .x = easting, .y = northing,
                .t = timestamp, crs = CRS(AEA.Africa.proj),
                order_by_ts = T)

# Add steplengths to object
wb2 <- wb2 %>% mutate(
  sl = step_lengths(.))

# Average step length?
summary(wb2$sl)

# 381.29 meters

# Summarize the distribution of time intervals between successive locations to get a general impression for the sampling rate.
summarize_sampling_rate(wb2)

# 1.67 hours


## ----Resample, eval=T, echo=T----------------------------------------------------------
wb2 <- wb2 %>% track_resample(rate = hours(3), tolerance = minutes(20))
wb2

# Ad-hoc thinning of Use locations - 20% removed
rcd.amt <- ceiling(nrow(wb2)*0.20)
wb2 <- wb2[sample(1:nrow(wb2),size = rcd.amt),c("x_","y_","t_")]

# Generating 50 available locations per each used location within animal homerange
wb2_RSF50 <- random_points(wb2,n=nrow(wb2)*50,typ="random")

# Stack rasters
rsf.stack <- stack(anth_risk,fence_dist,prirds_dist,secrds_dist,river_dist,waterpts_dist,woody_dist)

# Extracting the covariates
wb2_RSF50 <- wb2_RSF50 %>% extract_covariates(rsf.stack)

# Assess collinearity/Correlation
test.corr <- as.data.frame(wb2_RSF50)
vifstep(test.corr[,4:10]) # VIF indicates no collinearity problem, but there is a high level of correlation (max: 0.84)
cor(test.corr[,4:10])

# Once again, no collinearity issues.

# Create Model
M.2840 <- glm(case_ ~ scale(anth_risk) + scale(fence_dist) + scale(waterpts_dist) + scale(prirds_dist) + scale(river_dist) + scale(secrds_dist) +  + scale(woody_dist), family = binomial(link="logit"), 
     data = wb2_RSF50)
summary(M.2840)
confint(M.2840)

# Plot the response
visreg(M.2840,"anth_risk",
       scale="response", 
       ylab="Probability of Selection",
       xlab="Anthropogenic risk",
       partial=F,
       rug=F,
       line=list(col="black"), 
       fill=list(col="light gray"))

# The trend of anthropogenic disturbance is the same with animal 2840 as it is with 30077.  As anthropogenic risk increases, the probability of selection decreases.

# Create prediction - Extract the coefficient of each predictor from the model summary
coeff <- M.2840$coefficients

# Remember that we must scale the raster layers too:
anth.scale <- (anth_risk - mean(wb2_RSF50$anth_risk)) / sd(wb2_RSF50$anth_risk)
fence.scale <- (fence_dist - mean(wb2_RSF50$fence_dist)) / sd(wb2_RSF50$fence_dist)
water.scale <- (prirds_dist - mean(wb2_RSF50$prirds_dist)) / sd(wb2_RSF50$prirds_dist)
prirds.scale <- (secrds_dist - mean(wb2_RSF50$secrds_dist)) / sd(wb2_RSF50$secrds_dist)
river.scale <- (river_dist - mean(wb2_RSF50$river_dist)) / sd(wb2_RSF50$river_dist)
secrds.scale <- (waterpts_dist - mean(wb2_RSF50$waterpts_dist)) / sd(wb2_RSF50$waterpts_dist) 
woody.scale <- (woody_dist - mean(wb2_RSF50$woody_dist)) / sd(wb2_RSF50$woody_dist) 

# Prediction
pred <- exp(anth.scale*coeff[[2]]+ fence.scale*coeff[[3]] + water.scale*coeff[[4]] + prirds.scale*coeff[[5]] + river.scale*coeff[[6]] + secrds.scale*coeff[[7]] + woody.scale*coeff[[8]])/(1+exp(anth.scale*coeff[[2]]+ fence.scale*coeff[[3]] + water.scale*coeff[[4]] + prirds.scale*coeff[[5]] + river.scale*coeff[[6]] + secrds.scale*coeff[[7]] + woody.scale*coeff[[8]]))

# Provide Spatial Prediction - Based off of the coefficients from this single animal
plot(pred)

# Create new layer and crop it to the analysis extent
pred.2840 <- mask(pred, Athi_Bound)
pred.2840 <- crop(pred.2840,y=extent(Athi_Bound))
plot(pred.2840)


## ----tmap, eval=T, echo=T--------------------------------------------------------------
# Load prediction
load(file = "./Data/Prediction30077.Rdata")

# Create mean prediction
Mn.Predict <- mean(pred.2840, pred.30077)
plot(Mn.Predict)

# Create tmap_model
tmap_mode("view")
tm_basemap("OpenStreetMap") +
  tm_shape(Mn.Predict, name = "Habitat  suitability") +
  tm_raster(palette="-inferno", n=8, alpha=0.6, 
            title = "Wildebeest Habitat Suitability")

# Yes, the additional information from two animals, highlights differences in habitat suitability across the region.

