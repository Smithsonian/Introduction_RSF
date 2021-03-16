# Project: Introduction to Spatial Analyses
# Description: Animal movement analysis using wildebeest dataset as a template  
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

# In this example, we will use a dataset that I collected during my PhD research.  Data are freely available on Movebank:
# Stabach JA, Hughey LF, Reid RS, Worden JS, Leimgruber P, Boone RB (2020) Data from: Comparison of movement strategies of three populations of white-bearded wildebeest. Movebank Data Repository. doi:10.5441/001/1.h0t27719

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

# NOTE: Do NOT open the file in Excel before loading it into R.  Excel will mess up the time stamp column and set it to the time on your computer, not what you want for most studies.  This could be particularly important (erroneus) if you are interested in diurnal or nocturnal activity.

wild <- read.csv("./Data/White-bearded wildebeest in Kenya.csv", header=TRUE)

# This file was downloaded from Movebank, following procedures detailed during our class period.  Various datasets can be downloaded from Movebank, providing standardized field names that we will follow here.

# Look at the data and examine the data structure
head(wild)
str(wild)

# Dates
# ************************************************************
# ************************************************************

# It's essential that you properly format your dates and include the time zone
# Here, we will change the timesamp from UTC to EAT
Timezone1 <- "UTC" #Default
Timezone2 <- "Africa/Nairobi" 

# Look at the timestamp
wild$timestamp[1:50]

wild$timestamp <- as.POSIXct(wild$timestamp, format = "%Y-%m-%d %H:%M:%S", tz=Timezone1)
attr(wild$timestamp, "tzone")
wild$timestamp <- with_tz(wild$timestamp, tz=Timezone2)

# Look at the timestamp
wild$timestamp[1:50]

# Note that I have overwritten the timestamp field here, recorded in UTC
# Other timezones can be found at https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
# See the R documentation on time formats, from the strptime library
# https://www.rdocumentation.org/packages/base/versions/3.3/topics/strptime

# Check your attributes
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

# ************************************************************
# ************************************************************

# All of this could be completed using piping

wild <- read.csv("./Data/White-bearded wildebeest in Kenya.csv", header=TRUE) %>%
  filter(id == 2829,
         complete.cases(wild, c("location.lat","location.long","timestamp")]))) %>%
  
  mutate(month = lubridate::month(timestamp)) %>%
  filter(!duplicated(timestamp)) %>%
  arrange(desc(timestamp)) %>% #Added the order
  make_track(utm.easting, utm.northing, timestamp, id = id, name = individual.local.identifier, month = month, crs = sp::CRS("+init=epsg:32736"))
  
wild1 <- wild %>% filter(id == 2829)

wild1 <- make_track()

# Summarize the distribution of time intervals
summarize_sampling_rate(wild1)

stps <- track_resample(wild1, rate = minutes(10), tolerance = seconds(60)) %>%
                          steps_by_burst() %>%
                          time_of_day(include.crepuscule = FALSE)
str(stps)  

# ******************
# ******************
# Remove extra columsn
# Not absolutely necessary to do, but helpful to simplify the dataset.


# Create ID field and Re-arrange Columns
wild <- wild[,c(3:9,11,13:14,16:19)]

# Create an id field (some function use the tag.local.identified, so we'll keep that)
wild <- wild %>% 
  mutate(id = tag.local.identifier) %>% 
  relocate(id)

head(wild)

# Let's bring in another table with some ancillary data
# ******************************************************
wild.Anc <- read.csv("Data/Wildebeest_Ancillary.csv",header=T)

wild.Anc$Start.Date <- dmy(wild.Anc$Start.Date)
wild.Anc$End.Date <- dmy(wild.Anc$End.Date)
wild2 <- merge(wild,wild.Anc, by.x = "id", by.y = "ID")

# Subset or filter Mara animals
wild2 <- subset(wild, Study.Area == "Mara")

# Also see Amboseli Basin or Athi-Kaputiei Plains sub-regions. Potentially useful for student projects
unique(wild2$Study.Area)

# Let's choose 1 animals
wild2b <- subset(wild2, id == 2829)

# Same as
#wild2c <- wild2 %>%
#  filter(id == 2829)

# ********************
# ********************

# Determine if any missing data
if(all(complete.cases(wild2b)) == T){print("ALL Looking Good")} else {print("Dallas, we have a problem")}

# My main concern, however, is the Lat/Long field.  I will accept other fields being blank, as long as these fields are not blank
nrow(wild2b)
wild2b <- wild2b[complete.cases(wild2b[,c("location.lat","location.long","timestamp")]),]
# Or just query the entire dataset wild2b <- wild2b[complete.cases(wild2b),]
nrow(wild2b) # Not expecting any loss in rows, since above statement was true

# With GPS data, we have to be especially careful about duplicate timestamps.
any(duplicated(wild2b$timestamp))
wild2b <- wild2b[!duplicated(wild2b$timestamp), ]

# Create a month field
wild2b$month <- month(wild2b$timestamp)

# This dataset is already projected, so we can calculate tracks
tr1 <- make_track(wild2b, utm.easting, utm.northing, timestamp, id = id, name = individual.local.identifier, month = month, crs = sp::CRS("+init=epsg:32736"))

# ***************************
# ***************************
# This entire process could be condensed by using piping
wild2c <- wild %>% filter(id == 2829,
                          complete.cases(wild[,c("location.lat","location.long","timestamp")])) %>%
  mutate(month = lubridate::month(timestamp)) %>%
  filter(!duplicated(timestamp)) %>%
  arrange(desc(timestamp)) %>% #Added the order
  make_track(utm.easting, utm.northing, timestamp, id = id, name = individual.local.identifier, month = month, crs = sp::CRS("+init=epsg:32736"))

#make_track(location.long, location.lat, timestamp, id = id, name = individual.local.identifier, month = month, crs = sp::CRS("+init=epsg:4326"))
wild2c <- transform_coords(wild2c, sp::CRS("+init=epsg:32736"))

wild2c <- wild2c %>% mutate(sl_ = step_lengths(.))
summary(wild2c$sl_)

# How many records per individual?  
nrow(wild2c)
table(wild2c$name, dnn = "Records per individual")

# A note on projections.....I am not projecting the data here.  Instead, I am simply defining the projection.  For this projection (UTM36S, WGS84), this could also be specified as:
# sp::CRS("+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
# Make sure you understand the difference between "Defining your coordinate system" and "Projecting to a New coordinate System"

wild_geo <- transform_coords(tr1, sp::CRS("+init=epsg:4326"))
wild_geo <- wild_geo %>% mutate(sl_ = step_lengths(.))
summary(wild_geo$sl_)

# You could alternatively, take your Geographic Coordinates (your known/true data points) and project them:


# Quick graph of data
# Does it look correct?  If getting your data from Movebank, you can check the display to make sure things match.
plot(tr1, xlab = "Easting", ylab = "Northing", main = unique(tr1$id), bty="n")

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
