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

# We will now explore two different packages for analyzing animal movement data (adehabitat and amt).  We will start with visualizing the movement data, providing summary statistics, and initial graphs that can help to generate research questions.  
# We will then move onto RSF analyses and home range

# ************************************************************
# ************************************************************

# Important for our discussion is projecting the geographic data included in the datafile.  Note the difference between defining your projection and projecting to a new coordinate system.  Review Crego lecture for further details

# The best website to obtain Coordinate Reference System information is:
# https://spatialreference.org/

# Here, you can search for your coordinate system information and find the appropriate text required by R.
# ArcGIS or QGIS can also be used to help you define the parameters you are looking for (EPSG)

# In this case, we have Geographic Coordinate information (Lat/Long).  We will ignore the UTM data included in the file, because many of your datasets won't include this information.  Starting with Geographic coordinates will ensure that you have correctly annotated your dataset.

# Projection information
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326
UTM36s.proj <- "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # EPSG:32736

# In this case, we want to take Geographic data (Lat/Long) and "Project" to UTM 36S
# Create Matrix of positions
# In this case, the dataset already contains projected coordinates, so we could have used those.
wild.dd <- as.matrix(cbind(wild$location.long, wild$location.lat))
xy <- project(wild.dd, UTM36s.proj)

# Calculate animal trajectories
temp.traj <-as.ltraj(xy, wild$timestamp, id = wild$id, infolocs = wild[,1:ncol(wild)], typeII = TRUE, slsp = c("remove"))

summary(temp.traj)
temp.traj

# Note that this does not account for NA values
refda <- min(wild$timestamp)
refda <- strptime("00:00", "%H:%M", tz = Timezone2)
                              # or just "2010-05-30"
interval <- 3 # hour.  This is tricky because the data are 1 hour during the day and 3 hours at night

# These data are tricky because the fix interval is 1 hour during the day and 3 hours at night
# We need to set the missing values to NA
traj2 <- setNA(temp.traj, refda, tol = dt/10, dt = interval, units = "hour")
plotltr(traj2, "dt/3600")
is.regular(traj2)

## dt is nearly regular: round the date:
# Many of the functions in adehabitat require regular trajectories, to do so, use the sett0 function
traj2b <- sett0(traj2, refda, interval, units = "hour") 
plotltr(traj2b, "dt/3600")
is.regular(traj2b)

# Summarize trajectory
(Summary.traj <- summary(traj2))

# Summarize the completeness of the dataset
Summary.traj$DaysTrack <- signif(difftime(Summary.traj$date.end,Summary.traj$date.begin, units="days"),digits=2)
Summary.traj$Records <- Summary.traj$nb.reloc-Summary.traj$NAs
Summary.traj$PctComplete <- signif((Summary.traj$nb.reloc-Summary.traj$NAs)/Summary.traj$nb.reloc*100,digits=4)

# Convert to a dataframe 
Data.Traj <-ld(traj2)

# Plot trajectories
plot(traj2)



# Export date to shapefile
# **************************
# **************************
wild.2829.export <- SpatialPointsDataFrame(coords = xy, data = wild.2829, proj4string = CRS(UTM36s.proj))
writeOGR(obj = wild.2829.export, dsn = "Output", layer = "Wild2829_UTM36S", driver = "ESRI Shapefile", overwrite_layer = T)







# A note on projections.....I am not projecting the data here.  Instead, I am simply defining the projection.  For this projection (UTM36S, WGS84), this could also be specified as:
# sp::CRS("+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
# Make sure you understand the difference between "Defining your coordinate system" and "Projecting to a New coordinate System"


# Create a month field
wild2b$month <- month(wild2b$timestamp)

# This dataset is already projected, so we can calculate tracks
tr1 <- make_track(wild2b, utm.easting, utm.northing, timestamp, id = id, name = individual.local.identifier, month = month, crs = sp::CRS("+init=epsg:32736"))











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
