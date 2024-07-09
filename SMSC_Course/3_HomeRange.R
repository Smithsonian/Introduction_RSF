## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----Libraries, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------
# Remove everything from memory
rm(list=ls())

# Load required libraries
library(adehabitatHR)
library(ctmm)

# Timezone
Timezone2 <- "Africa/Nairobi" 

# Spatial reference information
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326
#https://spatialreference.org/ref/esri/102022/
AEA.Africa.proj <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"   #ESRI:102022


## ----Load----------------------------------------------------------------------------------------------------------------------------------------------
# Load .Rda object
load(file="./Data/wild.Rdata")

# Look at the data
head(wild.Athi)


## ----MCP, message = F, warning = F---------------------------------------------------------------------------------------------------------------------
# Define projection and reproject animal 30084
Athi.sp <- SpatialPointsDataFrame(coords = wild.Athi[c("location.long","location.lat")],
                                    data = wild.Athi,
                                    proj4string = CRS(LatLong.proj))

# Reproject to Albers
Athi.albers <- spTransform(Athi.sp, 
                         CRS =  CRS(AEA.Africa.proj))

# Subset datafrmae
Athi.30084 <- subset(Athi.albers, id == 30084)

# Compute MCP - Note that this is NOT a statistical estimate
Athi.mcp <- mcp(Athi.30084[,"id"], percent = 95, unin = "m", unout = "km2")
Athi.mcp

# Plot Result
sp::plot(Athi.30084, pch = '.', cex=2, main=paste0("Animal ",unique(Athi.30084$id)), xlab="Easting",ylab="Northing", axes = T, asp=1)
sp::plot(Athi.mcp, col = NA, border = "red", lwd = 2, add=T)


## ----Telemetry, message=F------------------------------------------------------------------------------------------------------------------------------
#  Import file into CTMM object, converting the spatial points object to a dataframe and using Albers Equal Area
Athi.sp <- as.data.frame(Athi.sp)

# Ctmm will recognize the location.long and location.lat files
# We will project to Albers Equal Area
wild <- as.telemetry(Athi.sp, timezone = TimeZone2, projection = AEA.Africa.proj, UERE=TRUE)

# How many animals?
length(wild)

# plots all in rainbow color
plot(wild, col=rainbow(length(wild))) # Okay, looks good

# Plot a particular animal...just need to specify the number in the list
# 1 is 30084, Karbolo
An <- 1
plot(wild[[An]], col="red")


## ----Variogram-----------------------------------------------------------------------------------------------------------------------------------------
# Set the dt, which in our case is hourly during the day (6 am - 6 pm) and every 3 hours at night (6 pm - 6 am), and confidence levels to display
dt <- c(1,3) %#% "hour" # The time interval of our data
level <- c(0.5, 0.95) # Confidence intervals placed on the contour estimates
xlim <- c(0,2 %#% "day") # 0-2 day window

# Variogram investigation
vg.karbolo <- variogram(wild[[An]], dt=dt)

# Plot full variogram
plot(vg.karbolo, level = level)

# Investigate different time lags
par(mfrow=c(1,2))
plot(vg.karbolo, xlim=xlim, level=level)
title("zoomed in")

plot(vg.karbolo, fraction = 0.65, level=level)
title("zoomed out")

# Variogram Model Selection, let the program take an initial guess based on the data
#GUESS <- ctmm.guess(wild[[An]], interactive=TRUE)
GUESS <- ctmm.guess(wild[[An]], interactive=FALSE)

# We can then use this GUESS to fit and evaluate a number of different movement models
FITS <- ctmm.select(wild[[An]], GUESS, verbose = TRUE, cores = 2) # Be patient, there's a lot happening in the background.

# Summarize the results in a human-readable format
summary(FITS) # OUF is the best model

# Summarize the model of interest
# This is super cool, but ignore the initial area reported here - that's the Gaussian estimate.  We'll refine that below
summary(FITS[[1]])

# Outputs from the model are Tau - position, Tau - velocity, and Speed (Km/day)
# Tau position is the estimate home range crossing time
# Tau Velocity is the Autocorrelation timescale (how long is the animal doing an activity before doing something different)

# Plot the variogram and the best fitting model
plot(vg.karbolo,FITS[[1]], xlim=xlim, level=level, col.CTMM = "blue")
title("zoomed in")

plot(vg.karbolo, FITS[[1]], fraction = 0.65, level=level, col.CTMM = "blue")
title("zoomed out")
par(mfrow=c(1,1))


## ----AKDE----------------------------------------------------------------------------------------------------------------------------------------------
# Calculate AKDE movement parameters and home range estimate
akde.est <- akde(data = wild[[An]], CTMM = FITS[[1]])
plot(wild[[An]], akde.est) # Plot the estimate

# Summarize AKDE output, which includes a maximum likelihood estimate and confidence intervals..
summary(akde.est)


## ----Multiple Animals, echo=F, eval=F------------------------------------------------------------------------------------------------------------------
## # For all animals......this is likely to take a long time for these animals...output to a list
## B.FITS <- list()
## B.AKDE <- list()
## 
## for(i in 1:length(wild))
##   {
##     print(i)
##     GUESS <- ctmm.guess(wild[[i]], interactive=FALSE)
##     B.FITS[[i]] <- ctmm.select(wild[[i]], GUESS, trace=2)
##     B.AKDE[[i]] <- akde(wild[[i]], B.FITS[[i]], trace=1)
## }
## 
## # Summarize the fits
## summary(B.FITS[[1]], units = FALSE)
## summary(B.AKDE[[1]], units = FALSE)
## 
## # Or all
## lapply(B.FITS, FUN = summary, units = FALSE) # This will make the units the same across individuals
## lapply(B.AKDE, FUN = summary, units = FALSE)
## 
## # Plot results
## plot(B.AKDE[[1]])
## 
## # Use metaanalysis to compare results
## meta(B.AKDE)
## plot(B.AKDE[[2]])

