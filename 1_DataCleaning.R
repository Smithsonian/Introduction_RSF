## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----Clean Libraries, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------
# Remove from memory
rm(list=ls())

# Load required libraries
library(adehabitatLT)
library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(proj4)


## ----Clean Timezone------------------------------------------------------------------------------------------------------------------------------------
# TimeZone and Projection
# Other timezones can be found at: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
Timezone1 <- "UTC"
Timezone2 <- "Africa/Nairobi" 
UTM36s.proj <- "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # EPSG:32736
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326


## ----Clean Read----------------------------------------------------------------------------------------------------------------------------------------
# Read files in Data directory (Data File and Accessory File) - FILE IS NOT SPATIAL
wild <- read.csv("./Data/White-bearded wildebeest in Kenya.csv", header=TRUE)
wild.ref <- read.csv("./Data/White-bearded wildebeest in Kenya-reference-data.csv", header=TRUE)

# Look at the data and examine the data structure
#head(wild)
#str(wild)

#head(wild.ref)
#str(wild.ref)


## ----Clean Merge---------------------------------------------------------------------------------------------------------------------------------------
# Re-organize fields and create new ID field (Not entirely necessary, but useful to remove the number of columns)
wild <- wild[,c(3:9,11,13:14)]

# Create an id field using the tag.local.identifier field
# Could also simply do: wild$id <- wild$tag.local.identifier, but then need to reorganize like we did above 
wild <- wild %>% 
  mutate(id = tag.local.identifier) %>% 
  relocate(id)
#head(wild)

# Look at dataframes
#str(wild)
#str(wild.ref)

# Grab the fields you are interested in merge (easist way)
# Merge the two datasets together (also see join in dplyr...lot's of options), including age, sex, and study study site from db2
wild.ref <- wild.ref[,c(1,4:5,7:8,17)]
wild <- merge(x= wild, y = wild.ref, by.x = "id", by.y = "tag.id", all.x = TRUE)
#head(wild)


## ----Clean TimeZone------------------------------------------------------------------------------------------------------------------------------------
# Look at the timestamp
#str(wild$timestamp)
wild$timestamp[1:10] # Character

# Format the timestamp from Timezone1 (UTC) to Timezone2 (EAT)
# Note that I have overwritten the timestamp field here, recorded initially in GMT/UTC (time zone/time standard)
wild$timestamp <- as.POSIXct(wild$timestamp, format = "%Y-%m-%d %H:%M:%S", tz=Timezone1) # Note the format, which is dependent on YOUR data
attr(wild$timestamp, "tzone")
wild$timestamp <- with_tz(wild$timestamp, tz=Timezone2) # Again, see: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
attr(wild$timestamp, "tzone")

# Look at the timestamp
wild$timestamp[1:10]

# Some additional examples for manipulating the time field, include:
# ****************************************
# ****************************************

# Example 1 (same as above, but using the lubridate package and piping):
# **********
#wild$ts2 <- ymd_hms(wild$study.local.timestamp, tz = Timezone2) # This is the really important step, because you need to properly capture the format
#head(wild$ts2)
# If wanted to change back to UTC
#wild$ts2 <- with_tz(wild$ts2, tz = Timezone1)
#head(wild$ts2)

# Or using piping and the dplyr package, you can do this in 1 step
# wild <- wild %>% mutate(tsEAT = ymd_hms(study.local.timestamp, tz = Timezone2), tsUTC = with_tz(tsEAT, tz = Timezone1))
# head(wild)

# Example 2 (Date and Time fields are in different columns):
# **********
# Fields must first be appended (use paste0 or lubridate).
# day <- c("2021-03-15", "2021-03-16", "2021-03-17", "2021-03-18", "2021-03-19")
# time <- c("06:02:01", "09:01:00", "12:00:00", "15:03:00", "18:01:01")
# test.data <- cbind.data.frame(day,time)
# test.data

# Using paste0
# test.data$timestamp <- as.POSIXct(paste0(test.data$day," ", 
#                                          test.data$time), 
#                                   format = "%Y-%m-%d %H:%M:%S", tz = Timezone1)
# test.data$timestamp

# Using lubridate
# test.data$timestamp2 <- as.POSIXct(lubridate::ymd(test.data$day) +
#                                      lubridate::hms(test.data$time), tz = Timezone1)
# test.data#timestamp2

# One other good hint:
# Use 'force_tz' to "Force" the zone without changing the actual clock time.


## ----Missing Final-------------------------------------------------------------------------------------------------------------------------------------
# Determine if any missing data
# In this particular clase, we're not expecting any NA values, since the data have been clean and uploaded from Movebank
if(all(complete.cases(wild)) == T){print("All Looking Good. No issues found")} else {
  print("!!! Dallas, We Have a Problem !!!")
  # Remove NA rows
  nrow(wild)
  wild <- wild[complete.cases(wild[,c("id","location.lat","location.long","timestamp")]),]
  nrow(wild)
  }


## ----Duplicated Final----------------------------------------------------------------------------------------------------------------------------------
# Duplicate timestamps:
if(anyDuplicated(wild[,c("id","timestamp")]) == F){print("No Duplicates found")} else {
  
  # Remove Duplicates
  print(paste0("Initial Number of Records in Dataset: ",nrow(wild)))
  print("Removing Duplicates")
  wild  <-  wild[!duplicated(wild[c('id', 'timestamp')]),] # or wild <- wild %>% distinct(id, timestamp)
  print(paste0("Number of Records after NA's removed: ",nrow(wild)))
}


## ----Different Start Dates-----------------------------------------------------------------------------------------------------------------------------
# Use the deploy.on/off dates in the dataframe to subset the data.
# Convert character field to a date field
wild$deploy.on.date <- ymd_hms(wild$deploy.on.date)
wild$deploy.off.date <- ymd_hms(wild$deploy.off.date)

# Based on Northrup et al, let's update the deploy.on.date + 1 day. 
# Alternatively, you could load in a dataframe with dates/times that you want to subset your dataframe
unique(wild$deploy.on.date)
wild$deploy.on.date <- wild$deploy.on.date+ days(1)
unique(wild$deploy.on.date)

# We can now use these dates to query the database
# There maybe a simpler way to do this, but I usually use a simple loop here, executing a date query on each animal.

# Create Unique ID
Un.ID <- unique(wild$id)

# Create Null Dataframe to store results
MyAnimals <- NULL

# Loop over each animal
for(i in 1:length(Un.ID)){
  # Subset by animal
  #temp <- subset(wild, id == Un.ID[i])
  temp <- wild[wild$id==Un.ID[i],] # This is simply a different way to subset or filter
  
  # Grab the unique Deploy.on.Date
  New.Start <- unique(temp$deploy.on.date)
  
  # Subset by Un.Start
  temp <- subset(temp, timestamp >= New.Start) # Add & clause if want to also subset by deploy.off.date (e.g., (timestamp >= New.Start & timestemp < New.End))
  
  # Bind the subset together
  MyAnimals <- rbind(MyAnimals, temp)
}

nrow(wild)
wild <- MyAnimals # Overwriting initial dataset
rm(MyAnimals) # Remove unneeded dataframe
nrow(wild)

# Remove unwanted Deploy Date columns
#head(wild)
wild <- wild[,-(12:13)]


## ----Summary Load Data---------------------------------------------------------------------------------------------------------------------------------
# Plot the DOP values for confirmation.
# For this dataset, the included identifying whether the position was 2D or 3D and then using a qualitative measure to remove poor positions
# I was more restrictive on 2D positions (gps.dop < 5.0)

# Example: 
nrow(wild)
hist(wild$gps.dop, xlab="DOP", ylab="Frequency", main="Wildebeest GPS Data")

wild <- subset(wild, gps.fix.type.raw == "3D" & gps.dop < 10.0 | gps.fix.type.raw == "2D" & gps.dop < 5.0) # This doesn't do anything here because I already filtered the data

hist(wild$gps.dop, xlab="DOP", ylab="Frequency", main="Wildebeest GPS Data")
nrow(wild)
# Important to inspect your own dataset to understand the fields included


## ----Create Trajectory---------------------------------------------------------------------------------------------------------------------------------
# Create matrix
temp <-as.matrix(cbind(wild$location.long,wild$location.lat))

# Project to UTM 36S (projection specified above) or other meter projection
xy <- project(temp, UTM36s.proj) # This uses the proj4 package (Spatial Points File)
#plot(xy)

# Use the as.ltraj function to create individual animal trajectories. Use the infolocs command to included attributes of the dataframe.  
# Automatically calculates the distance between succesive locations, relative and absolute turning angles (in radians), and the time interval between successive locations (in seconds)
traj.raw <-as.ltraj(xy, date = wild$timestamp,id = wild$id, typeII = TRUE, infolocs = wild[3:14], slsp = c("remove"))

# This is a list of objects
# To view individual aniamls
head(traj.raw[[1]])

# Notice that the dataset has no NAs.  This is because the function doesn't recognize the movement interval, something we must set.
# In our case, the data were collected: Every hour from 6 am to 6 pm and Every three hours from 6 pm to 6 am
# For the remainder of our exercises, we will treat these data as a 3 hour dataset

# Create a reference date and use setNA to re-run trajectory
refda <- strptime("00:00", "%H:%M", tz=Timezone2)

# Create NA values and make a regular trajectory based on refda
traj.NA <- setNA(traj.raw, refda, 3, units = "hour") 
traj.reg <- sett0(traj.NA, refda, 3, units = "hour")
#is.regular(traj.reg)

# Summarize Trajectory
Summary.traj <- summary(traj.reg)

# Add details to the summary
Summary.traj <- Summary.traj %>% mutate(
  DaysTrack = signif(difftime(date.end,date.begin, units="days"),digits=2),
  Records = nb.reloc-NAs,
  PctComplete = signif((nb.reloc-NAs)/nb.reloc*100,digits=4),
)

# Look at
Summary.traj

# Convert trajectory with movement statistics to a dataframe
wild.df <-ld(traj.reg)

# Calculate basic movement statistics.  This are just some initial statistics
# REMEMBER: this is a 3 hour data.  Therefore, movements are simply based on the linear steps.
Mvmt.Statistics <- ddply(wild.df,"id", summarise,
              AvgMove = round(mean(dist/1000,na.rm=TRUE),digits=2), # Convert to km
              SumMove = round(sum(dist/1000,na.rm=TRUE),digits=2), # Convert to km
              MaxDisp=round(max(sqrt(R2n)/1000,na.rm=TRUE),digits=2)) # Convert to km
Mvmt.Statistics


## ----Visualize-----------------------------------------------------------------------------------------------------------------------------------------
# Plot trajectory using adehabitat object (list)
plot(traj.reg[1]) # Or plot(traj.wild) to view all animals

# We can also look at the DOP over time or the how the data were collected (every X minutes)
plotltr(traj.reg[1], "gps.dop") # Graphic of DOP over time......these should all be < 10, since we've cleaned them above
#plotltr(traj.reg[1],"dt/60") 

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


## ----Clean Export--------------------------------------------------------------------------------------------------------------------------------------
# Filter the data to the Mara region for analysis
wild.Mara <- wild %>% filter(
  study.site == "Mara")

wild.Athi <- wild %>% filter(
  study.site == "Athi-Kaputiei Plains")

# How many animals are included in the dataframe?  Should be 15.
length(unique(wild.Mara$id))
length(unique(wild.Athi$id))

# Save file to your Data directory for subsequent analyses
#save(wild.Mara, wild.Athi, file = "./Data/wild.Rdata")

