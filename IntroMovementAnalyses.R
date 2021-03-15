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

# ************************************************************
# ************************************************************

# Remove everything from memory
rm(list=ls())

# Load Libraries
library(ctmm)
library(lubridate)
library(dplyr)
library(raster)
library(proj4)
library(adehabitatLT)
library(sp)
library(amt)

# Load and Clean Data
# ************************************************************
# ************************************************************

# In this example, we will use a dataset that I collected during my PhD research.  Data are freely available on Movebank:
# Stabach JA, Hughey LF, Reid RS, Worden JS, Leimgruber P, Boone RB (2020) Data from: Comparison of movement strategies of three populations of white-bearded wildebeest. Movebank Data Repository. doi:10.5441/001/1.h0t27719

# Do NOT open the file in Excel before loading it into R.  Excel will mess up the time stamp column and set it to the time on your computer, not what you want for most studies.  This could be particularly important (erroneus) if you are interested in diurnal or nocturnal activity.
wild <- read.csv("./Data/White-bearded wildebeest in Kenya.csv", header=TRUE)
head(wild)
str(wild)

# Set TimeZone
# https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
TimeZone <- "Africa/Nairobi"

# Note that I am overwriting the timestamp field here (UTC)
# See the R documentation on time formats, from the strptime library
# https://www.rdocumentation.org/packages/base/versions/3.3/topics/strptime
wild$timestamp <- as.POSIXct(wild$study.local.timestamp, format = "%Y-%m-%d %H:%M:%S", tz=TimeZone)

# Check your attributes
# It is essential that you are certain of your time zone
attr(wild$timestamp, "tzone")

# Note, you could have also taken the UTC time, created a POSIXct object and then converted to local time
# wild$timestamp <- as.POSIXct(wild$timestamp, format = "%Y-%m-%d %H:%M:%S", tz="UTC")
# attr(wild$timestamp, "tzone")
# wild$timestamp <- with_tz(wild$timestamp, "Africa/Nairobi")
# attr(wild$timestamp, "tzone")

# The lubridate package has a number of useful functions for time, especially when reading data from a file when your system will assume a certain timezone.
# force_tz changes the time zone without changing the clock time.

# See also some of the manipulations you can do to your time stamp, especially when it is a character field
#wild$test <- ymd_hms(wild$study.local.timestamp, tz = TimeZone)
#wild  <- wild %>% mutate(test = ymd_hms(study.local.timestamp, tz = TimeZone)) # Could use mutate to get the same result

# Your date and time may also be in two different fields.  In those, cases, you need to append them together.
# You can use the paste0 command or use lubridate
#day <- c("2021-03-15", "2021-03-16", "2021-03-17", "2021-03-18", "2021-03-19")
#time <- c("06:02:01", "09:01:00", "12:00:00", "15:03:00", "18:01:01")
#test.data <- cbind.data.frame(day,time)
#test.data$timestamp <- as.POSIXct(paste0(test.data$day," ",
#                                         test.data$time), 
#                                  format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
#test.data$timestamp2 <- as.POSIXct(lubridate::ymd(test.data$day) +
#                      lubridate::hms(test.data$time), tz = "UTC")
#See lubridate for additional options
#https://cran.r-project.org/web/packages/lubridate/lubridate.pdf

# ******************
# ******************
# Remove extra columsn
# Not absolutely necessary to do, but helpful to simplify the dataset.


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

# ***************************
# ***************************
# Load rasters...projected to UTM 36S, WGS84

# List the Raster files for import
raster.list <- list.files(path ="./Data/SpatialData/ASCII/", pattern = ".asc", full.names = TRUE)
raster.list

# Import as a raster stack
All.rasters <- stack(raster.list)
plot(All.rasters)
projection(All.rasters)

# ***************************
# ***************************
# CTMM Home range analysis
# Data could be directly loaded into CTMM for homerange analysis
wild.ctmm <- as.telemetry(wild)
