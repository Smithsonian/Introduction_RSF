# Project: Wildebeest
# Description: Introduction to Animal Movement Analyses using wildebeest dataset as a template  
# Author: Jared Stabach and Majaliwa Masolele
# Date: 11-Mar-2021

# ********************************************************************************
# ********************************************************************************

# Goals of the Analysis

#1. Import data and clean
#2. Visualize data, provide summary plots of movement - Interested to see how you do these two steps.  Happy to merge my code with yours
#3. Calculate Home range (MCP and AKDE) - I have code for this
#4. RSF analysis and potentially SSF (iSSF) analysis - This is intro, don't want to get to complicated
#     Start with a single animals
#     Integrate/extract a few data layers (Distance to woody cover, anthropogenic disturbance)
#     Build to multi-animal model
#     Start with simple RSF (logistic regression model).  I'd like to build a model with 2 or 3 varables and then ask students to construct a new model with an additional 1 or 2 variables.  I'd then like them to tell me which model is better.
#     Perhaps move to SSF (conditional logistic regression) or even iSSF, but I view this as much more advanced
#5. Animate Movements - I have code for this

# I am going to put all of this in a Markdown document as it develops and will credit you for your assistance

# *****************************************************************************
# *****************************************************************************

# Remove everything from memory
rm(list=ls())

# Load Libraries
library(ctmm)
library(lubridate)
library(dplyr)
require(ggplot2)
library(raster)
library(proj4)
library(adehabitatLT)
library(sp)
library(amt)

# Load and Clean Data
# *****************************************************************************
# *****************************************************************************

# In this example, we will use a dataset that I collected during my PhD research.  Data are freely available on Movebank:
# Stabach JA, Hughey LF, Reid RS, Worden JS, Leimgruber P, Boone RB (2020) Data from: Comparison of movement strategies of three populations of white-bearded wildebeest. Movebank Data Repository. doi:10.5441/001/1.h0t27719

# Do NOT open the file in Excel before loading it into R.  Excel will mess up the time stamp column and set it to the time on your computer, not what you want for most studies 
load("./Data/wildebeest.rda")
head(wild)
str(wild)

# Set TimeZone
# https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
TimeZone <- "Africa/Nairobi"

# Note that I am overwriting the timestamp field here (UTC)
wild$timestamp <- as.POSIXct(wild$study.local.timestamp, format = "%Y-%m-%d %H:%M:%S", tz=TimeZone)

# Check your attributes
# It is essential that you are certain of your time zone
attr(wild$timestamp, "tzone")

# Note, you could have also taken the UTC time, created a POSIXct object and then converted to local time
# wild$timestamp <- as.POSIXct(wild$timestamp, format = "%Y-%m-%d %H:%M:%S", tz="UTC")
# wild$timestamp <- with_tz(wild$timestamp, "Africa/Nairobi")

# The lubridate has a number of useful functions for time, especially when reading data from a file when your system will assume a certain timezone.
# force_tz changes the time zone without changing the clock time.

# See also some of the manipulations you can do to your time stamp, especially when it is a character field
# I would recommend always specifying your TimeZone
#wild$test <- ymd_hms(wild$study.local.timestamp, tz = TimeZone)
#wild  <- wild %>% mutate(test = ymd_hms(study.local.timestamp, tz = TimeZone)) # Could use mutate to get the same result

# Movebank helps in resolving these issues because it standardizes the columns and formats.

# I had to first open the file in Excel to  Must fix the time, first by correcting the time in the Excel file
# I then needed to set the format of the time to POSIXct using lubridate
# Find the function that matches the format of your data


# There are extra columns in the data frame that we can immediately remove.  Not absolutely necessary to do, but helpful to simplify the dataset.


# Create ID field and Re-arrange Columns
wild <- wild[,c(3:9,11,13:14,16:19)]
wild <- wild %>% 
  mutate(id = tag.local.identifier) %>% 
  relocate(id)

# Let's restrict to a few animals in the Maasai Mara.  Animals that I know that fit the home range assumption.
Mara.An <- c(2829, 2830, 2831, 2832, 2833, 2834, 2835, 2836, 2838, 2839, 2841, 2843, 2844, 2845, 2846)
#Ambo.An <- c(2840, 2842, 30068, 30070, 30071, 30072, 30074, 30077, 30079, 30082, 30084, 30086)
#Athi.An <- c(2837, 30069, 30073, 30075, 30076, 30078, 30081, 30083, 30085)

# Subset to Mara animals
wild2 <- wild[(wild$id %in% Mara.An),]

# Visualize data, provide summary plots of movement 
# ***************************
# ***************************
# check that there are no outliers in x-y directions
plot(wild2$utm.easting,wild2$utm.northing)

# remove any NA's in timestamp column
wild2<-wild2[!is.na(wild2$timestamp),]

# ************************ 4.	Fit a simple RSF model **************************
# *****************************************************************************
# *****************************************************************************
# Here we are going to use amt package which has functions to accomplish this task

# load all raster layers to be used as the predictor variables
anth_risk<-raster("./Data/SpatialData/ASCII/anth_risk.asc") # Anthropogenic disturbance
woody_dist<-raster("./Data/SpatialData/ASCII/woody_dist.asc") # Woody cover

#Select only a single individual-2829
wb_2829<-wild2[wild2$id==2829,]

#Retain only essential columns
wb_2829<-wb_2829[,c(1,2,12,13)]

# First we need to create track object using mk_track function
wb_2829<-mk_track(wb_2829,.x=utm.easting,.y=utm.northing,
                .t=timestamp,crs = sp::CRS("+init=epsg:21036"))

# Summarize the distribution of time intervals between successive locations
# to get a general impression for the sampling rate.
summarize_sampling_rate(wb_2829)

# Generating available points (here we are generating 5 available points using the home range level of the 
# minimum convex polygon. This function conveniently calculates a minimum convex polygon around Lupe’s locations. 
# It then samples 5 random locations for every used location within this polygon
wb_RSF<-random_points(wb_2829,n=nrow(wb_2829)*5,typ="random")

# Adding the covariates
wb_RSF<-wb_RSF %>% extract_covariates(anth_risk)%>%
  extract_covariates(woody_dist)

# Fitting model
M1<- glm(case_ ~ scale(anth_risk) + scale(woody_dist), 
     family = binomial(), data =wb_RSF)

# Exploring sensitivity of RSF coefficients to the number of available points
# *****************************************************************************
# *****************************************************************************
# Here we evaluate how coefficient estimates changes as we increase the number of available locations from 1 available location per used location 
# to 100 available locations per used location.

# Available location
n.frac <- c(1, 5, 20, 50, 100) 

# Total available locations to be generated
n.pts <- ceiling(nrow(wb_2829) * n.frac)

# Replication of each available location (1,5,20,50,100) 20 times
n.rep <- 20

# Here we create a table which contains 20 replications of generated available locations,
# followed by extracting covariates and finally fit
# glm model for each available location (1,5,20,50,100) and for each replication
wb <- tibble(
  n.pts = rep(n.pts, n.rep), 
  frac = rep(n.frac, n.rep), 
  res = map(
    n.pts, ~
      wb_2829%>% random_points(n = .x) %>% 
      extract_covariates(anth_risk)%>%
      extract_covariates(woody_dist) %>% 
      mutate(anth_risk= scale(anth_risk), 
             woody_dist = scale(woody_dist)) %>% 
      glm(case_ ~ anth_risk + woody_dist, 
          data = ., family = binomial()) %>% 
      tidy()))

# Plotting the coefficient estimates from the model 
# for each available location(1,5,20,50,100) and for each predictor
wb %>% unnest(cols = res) %>% 
  mutate(wb2829 = recode(term, "(Intercept)"="Intercept",
                         anth_risk = "Anthropogenic disturbance", 
                         woody_dist = "Woody distance",)) %>% 
  ggplot(aes(factor(frac), y = estimate)) +
  geom_boxplot() + facet_wrap(~ wb2829, scale  ="free") +
  geom_jitter(alpha = 0.2) + 
  labs(x = "Number of available points (X observed locations)", y = "Estimate") +
  theme_light()

# Here the intercept decreases as the number of available points increases,
# but the slope parameter estimates, on average, do not change much once we included at least 20 available points per used point. 
# Thus, we conclude that, in this particular case, having 20 available points per used point is sufficient for interpreting the slope coefficients.

# ******** Here we decided to proceed with 20 available locations **************
# ******************************************************************************
# Generating 20 available locations per each used location
wb_RSF20<-random_points(wb_2829,n=nrow(wb_2829)*20,typ="random")

# Adding the covariates
wb_RSF20<-wb_RSF20 %>% extract_covariates(anth_risk)%>%
  extract_covariates(woody_dist)

# Fitting model
M2<- glm(case_ ~ scale(anth_risk) + scale(woody_dist), 
         family = binomial(), data =wb_RSF20)

# Get summary of the model
summary(M2)

# Getting confidence interval of the model
confint(M2)

# ******************** Multiple individual RSF model ***************************
# *****************************************************************************

#Subsetting only important column we need (timestamp, animal id, latitude and longitude)
wild2<-wild2[,c(1,2,12,13)]

# Creating list that contain that dataframe for 15 individuals we subsetted above for Mara
wild2<-split(wild2,wild2$id)

# Creating track for each individual animal in a list
wild2k<-lapply(wild2, function(x) mk_track(x,
         .x=utm.easting,.y=utm.northing,
         .t=timestamp,crs = sp::CRS("+init=epsg:21036")))

#Generating random points 20 x used units for each animals 
wb15<-lapply(wild2k, function(x) random_points(x,n=nrow(x)*20,typ="random"))

# ************ Extract covariates*******

# Adding anthropogenic risk
wb15<-lapply(wb15, function(x) extract_covariates(x,anth_risk))

# Adding woody distance variable
wb15<-lapply(wb15, function(x) extract_covariates(x,woody_dist))

# Adding animal ID column to each dataframe in a list using the ID we subseted above
wb15<- mapply(cbind, wb15, "AnimalID"=Mara.An, SIMPLIFY=F)

# Binding all dataframe in a list into a single dataframe for all animals
wb15<-do.call(rbind,wb15)

##************ Fitting model with generalised linear mixed model ***************
require(lme4)

mixed<-glmer(case_ ~ scale(anth_risk) + scale(woody_dist) + (1|AnimalID), 
             family = binomial(), data =wb15)

summary(mixed)

# ****************************** References ************************************
# Signer, J., Fieberg, J. and Avgar, T., 2019. Animal movement tools (amt):
# R package for managing tracking data and conducting habitat selection analyses.
# Ecology and evolution, 9(2), pp.880-890.

# Fieberg, J., Signer, J., Smith, B. and Avgar, T.,2021. A ‘How-to’Guide for 
#Interpreting Parameters in Habitat-Selection.Journal of Animal Ecology.
# ******************************************************************************