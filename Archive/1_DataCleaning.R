## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----Clean Libraries, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------------------------
# Remove from memory
rm(list=ls())

# Load required libraries
library(adehabitatLT)
library(lubridate)
library(move)
library(proj4)
library(plyr)
library(tidyverse)


## ----Clean Timezone------------------------------------------------------------------------------------------------------------------------------------------------------------
# TimeZone and Projection
# Other timezones can be found at: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
Timezone1 <- "UTC"
Timezone2 <- "Africa/Nairobi" 
UTM36s.proj <- "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # EPSG:32736
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326


## ----Clean Read, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------------------------------------
# Read files in Data directory (Data File and Accessory File) - FILE IS NOT SPATIAL
wild <- read_csv("./Data/White-bearded wildebeest in Kenya.csv")
wild.ref <- read_csv("./Data/White-bearded wildebeest in Kenya-reference-data.csv")

# In this case, the column headers have "-" or ":" in them.  These need to be removed:
names(wild) <- gsub("-", "", names(wild))
names(wild) <- gsub(":", "", names(wild))
names(wild.ref) <- gsub("-", "", names(wild.ref))

# NOTE: When printing to screen, it appears that significant digits have been truncated.  This is not the case.  It is just how R prints the data
# wild
# wild$locationlong[1]

# Look at the data and examine the data structure
#head(wild)
#str(wild)

#head(wild.ref)
#str(wild.ref)

# If you have multiple files, rather than a single file with all your animals, you can do something similar to:
# Change the wildcard (pattern) to what makes most sense to your dataset
#wild.all <- list.files(path='./data', pattern='Kenya.csv', all.files=FALSE, full.names=TRUE)

# Bind files together
#wild <- lapply(wild.all, FUN=read_csv)
#wild <- do.call(rbind, wild)

# In most cases, you won't need to do this, but this file need a few things to be fixed
#names(wild) <- gsub("-", "", names(wild))
#names(wild) <- gsub(":", "", names(wild))



## ----Clean Merge---------------------------------------------------------------------------------------------------------------------------------------------------------------
# Clean the reference file, selecting only the columns that you want to include
# Important is that the file contains the deployment on and off date and times. The dates will allow us to filter/remove datapoints outside of the deployment window
wild.ref <- 
  wild.ref %>% 
  transmute(id = tagid,
            start = deployondate,
            end = deployoffdate,
            sex = animalsex,
            site = studysite)

# Now do the same for the tracking dataset
# Rename/Reorganize column headers, keeping only the fields of interest
# Change timezone to local time
wild <- 
  wild %>% 
  transmute(id = taglocalidentifier,
            name = individuallocalidentifier,
            timeUTC = timestamp, # Keep if you want for comparison
            timeLocal = with_tz(timestamp, tz = Timezone2), # Update timezone. Field must be a date-time object (POSIXct or other)
            longitude = locationlong,
            latitude = locationlat,
            temp = externaltemperature,
            dop = gpsdop,
            fixtype = gpsfixtyperaw,
            height = heightaboveellipsoid) %>% 
  
  # Make sure no duplicate id and timestamp exist.
  distinct(id, timeLocal, .keep_all = TRUE) %>% 

  # Remove any records that don't have a timestamp or a Lat/Long location
  filter(!is.na(timeLocal),
         !is.na(latitude)) %>% 

  # Join the Reference dataset to the movement dataset
  left_join(
    wild.ref %>% 
      dplyr::select(id, sex, site, start, end),
    by = c('id' = 'id')) %>% 

  # Now use the start and end timestamps to filter the dataset
  # Do this specifically for each individual (group_by)
  group_by(id) %>% 
  filter(
    timeLocal > start & timeLocal < end) %>% 
  
  # Remove fields you don't need
  dplyr::select(-c(start, end)) %>% 
  
  # Arrange the dataset by id and timestamp
  arrange(id, timeLocal)
  


## ----Summarize, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------
# Look at the data and summarize results
wild %>% 
    
  # Check on which animals included and timezone
  group_by(id) %>% 
  
  # Summarize the dataset to highlight the number of records and duration of tracking period
  summarize(
    Records = n(),
    #MinLat = min(latitude), MaxLat = max(latitude), MinLong = min(longitude), MaxLong = max(longitude),
    StartDate = min(timeLocal),
    EndDate = max(timeLocal), 
    TrackPeriod_yrs = round(as.duration(min(timeLocal) %--% max(timeLocal)) / dyears(1), digits=1)) %>% 
  
  # Arrange by start date
  arrange(StartDate, id)

# Create plot
# Not pretty at this point, but this is just a start
ggplot(
  data = wild,
  mapping = aes(
    x = longitude,
    y = latitude,
    color = as.factor(id))) + 
  geom_point(shape = ".",
             size = 3.5,
             alpha = 0.25)


## ----Summary Load Data, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------------------
# Plot the DOP values for confirmation.
# For this dataset, we will examine whether the position is 2D or 3D and then use a qualitative filter to remove poor positions
# Here, I am being more restrictive on 2D positions (dop < 5.0)

# Example: 
nrow(wild)
ggplot(
  data = wild,
  aes(x = dop)) +
  geom_histogram(color = "black", fill = "white") +
  labs(title="Wildebeest GPS Data", x = "DOP", y = "Frequency") +
  theme_classic()

# Filter/subset the dataset to remove problem records 
# This won't do much here because most of these data have already been filtered
wild <- 
  wild %>% 
  filter(
    fixtype == "3D" & dop < 10.0 | fixtype == "2D" & dop < 5.0)

nrow(wild)

# Could then plot again to show the difference.
# I'm not doing so here because I have already cleaned the data. included

# ggplot(
#   data = wild,
#   aes(x = dop)) +
#   geom_histogram(color = "black", fill = "white") +
#   labs(title="Wildebeest GPS Data", x = "DOP", y = "Frequency")
#   theme_classic()


## ----Create Trajectory, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------------------
# Create matrix
temp <-as.matrix(cbind(wild$longitude,wild$latitude))

# Project to UTM 36S (projection specified above) or other meter projection
xy <- project(temp, UTM36s.proj) # This uses the proj4 package (Spatial Points File)

# Calculate the distance between successive locations, relative and absolute turning angles (in radians), and the time interval between successive locations (in seconds)
# Error message relates to using a tibble instead of a dataframe
traj.raw <-as.ltraj(xy, 
                    date = wild$timeLocal,
                    id = wild$id, 
                    typeII = TRUE, 
                    infolocs = wild[5:8], # columns latitude, longitude, temp, dop
                    slsp = c("remove"))

# Look at all animals
plot(traj.raw)
plot(traj.raw[1]) # Or any other animal
traj.raw 

# Notice that the dataset has no NAs.  This is because the function doesn't recognize the movement interval, something we must set.
# In our case, the data were collected: Every hour from 6 am to 6 pm and Every three hours from 6 pm to 6 am
# For the remainder of our exercises, we will treat these data as a 3 hour dataset

# Create a reference date and use setNA to re-run trajectory
refda <- strptime("00:00", "%H:%M", tz=Timezone2)

# Create NA values and make a regular trajectory based on refda
traj.NA <- setNA(traj.raw, refda, 3, units = "hour") 

# Summarize Trajectory
Summary.traj <- summary(traj.NA)

# Add details to the summary
Summary.traj <- 
  Summary.traj %>% 
  mutate(
    DaysTrack = signif(difftime(date.end,date.begin, units="days"),digits=2),
    Records = nb.reloc-NAs,
    PctComplete = signif((nb.reloc-NAs)/nb.reloc*100,digits=4),
)

# Look at
Summary.traj

# Convert trajectory with movement statistics to a dataframe
wild.df <-ld(traj.NA)

# Calculate basic movement statistics.  This are just some initial statistics
# REMEMBER: this is a 3 hour data.  Therefore, movements are simply based on the linear steps.
Mvmt.Statistics <- ddply(wild.df,"id", summarise,
              AvgMove = round(mean(dist/1000,na.rm=TRUE),digits=2), # Convert to km
              SumMove = round(sum(dist/1000,na.rm=TRUE),digits=2), # Convert to km
              MaxDisp=round(max(sqrt(R2n)/1000,na.rm=TRUE),digits=2)) # Convert to km
Mvmt.Statistics


## ----Visualize-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot trajectory using adehabitat object (list)
plot(traj.NA[1])

# We can also look at the DOP over time or the how the data were collected (every X minutes)
plotltr(traj.NA[1], "dop") # Graphic of DOP over time......these should all be < 10, since we've cleaned them above
#plotltr(traj.NA[1],"dt/60") 

# My preference is to create a custom plot
# Setup a plotting layout with three panels
layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), widths=1, heights=c(1,1))

# Here, we will "loop" over every individuals (by id)
Id.val <- unique(wild.df$id)

#for (i in 1:length(Id.val)){  # Uncomment out
i <- 1 ###### Remove to loop over all animals
# ****************************** 

  # Remove the NAs and subset to id == Id.val[i]
  wild.sub <- subset(wild.df[!is.na(wild.df$x),], id == Id.val[i]) # This is just to account for NAs that have been added from as.ltraj function
  
  # Calculate the total days tracked
  time.diff <- trunc(difftime(format(wild.sub$date[1],tz=Timezone2),format(wild.sub$date[nrow(wild.sub)],tz=Timezone2),units="days"))
  
  # Plot the trajectory
  plot(wild.sub$x,wild.sub$y,typ="l",xlab="Easting",ylab="Northing",main=paste0(wild.sub$id[1]," Movement"),frame=FALSE,axes=FALSE,asp=1)
     mtext(paste0(format(wild.sub$date[1],"%Y-%m-%d")," to ",format(wild.sub$date[nrow(wild.sub)],"%Y-%m-%d")),cex=0.75)
     axis(1, labels=TRUE)
     axis(2, labels=TRUE)
  points(wild.sub$x,wild.sub$y,pch=16,cex=0.5,col="blue")
  points(wild.sub$x[1],wild.sub$y[1],pch=17,cex=1,col="green")
  points(wild.sub$x[nrow(wild.sub)],wild.sub$y[nrow(wild.sub)],pch=15,cex=1,col="red")
  
  # Plot the movements over time (Velocity)
  plot(wild.sub$date, wild.sub$dist/1000, type='l', ylab="Distance moved (km)", xlab="Time", main="Steplengths", frame=FALSE)
    # Calculate the time from release date	
    mtext(paste0(abs(time.diff)," days"),cex=0.75)
  
  # Plot the net displacement per step
  plot(wild.sub$date, sqrt(wild.sub$R2n)/1000, type='l', ylab="Distance (km)", xlab="Days Since Release", main="Net Displacement",frame=FALSE)
    mtext(paste0(abs(time.diff)," days"),cex=0.75)
#} # Uncomment


## ----Clean Export--------------------------------------------------------------------------------------------------------------------------------------------------------------
# Filter the data for analysis
wild.Mara <- 
  wild %>% 
  filter(
    site == "Mara")

wild.Athi <- 
  wild %>% 
  filter(
    site == "Athi-Kaputiei Plains")

wild.Ambo <- 
  wild %>% 
  filter(
    site == "Amboseli Basin")

# How many animals are included in the dataframe?  Should be 15, 12, and 9.
length(unique(wild.Mara$id))
length(unique(wild.Athi$id))
length(unique(wild.Ambo$id))

# Save file to your Data directory for subsequent analyses
save(wild.Mara, wild.Athi, wild.Ambo, file = "./Data/wild.Rdata")

