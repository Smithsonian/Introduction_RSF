# ***************************************************************************************************
# ***************************************************************************************************
# Project: Maasai Mara landscape and wildlife changes - Movement Animation
# Description: Following guidelines/code provided by Jakob Schwalb-Willmann
#       Using Wildebeest as example datasete
# Author: Jared Stabach
# Date: 15 October 2020

# ***************************************************************************************************
# ***************************************************************************************************

# Clear the cache
rm(list=ls())

# Load libraries
library(moveVis)
library(move)
library(raster)
library(MODIS)

# Importing and subsetting trajectories
load("Movement_Animation/Data/Mara_Test.Rda")
head(data.Mara)
unique(data.Mara$ID)

# Drop levels
data.Mara <- droplevels(data.Mara)

# Query dataset
colnames(data.Mara)
range(data.Mara$timestamp)
unique(data.Mara[["individual.local.identifier"]])

# Converting to a move class
m <- move(x = data.Mara[["Long"]], y = data.Mara[["Lat"]],
          time = data.Mara[["timestamp"]], animal = data.Mara[["Name"]],
          proj = "+proj=longlat +datum=WGS84 +no_defs",
          removeDuplicatedTimestamps = TRUE)

# Check time stamps and sampling rates (temporal resolution) of the trajectories:
lag <- unlist(timeLag(m, unit = "mins"))
median(lag)
sd(lag)

# In general, animals were monitored every hour during day (6 am to 6 pm) and every three hours at night (6 pm to 6 am)
# moveVis needs to assign each location of a trajectory to a specific frame, the sampling times of all locations across all trajectories need to be aligned to share a uniform temporal resolution and uniform time stamps that can be assigned to frames. moveVis includes a dedicated function to align trajectories using linear interpolation, named align_move():

# Aligning to a 3 hour interval
m <- align_move(m, res = 180, digit = 0, unit = "mins")
length(unique(timestamps(m)))

# subset by character times
m <- subset_move(m, from = "2010-05-27 12:00:00", to = "2010-08-01 12:00:00")

# check min and max of result
min(timestamps(m))
max(timestamps(m))

# Animating trajectories on a static OpenStreetMap base map
# Create new extent
ext <- extent(m) * 1.15
ext@ymin <- ext@ymin * 1.05
ext@ymax <- ext@ymax + abs(ext@ymax * 0.05)

# The custom extent can be passed to frames_spatial() via the ext argument.
frames <- frames_spatial(m, trace_show = TRUE, equidistant = FALSE,
                         ext = ext, map_service = "osm",
                         map_type = "terrain_bg")

frames[[300]]

# moveVis provides a set of functions to customize the appearance of frames after they have been created. For example, the user can add title, caption and axis texts using add_labels(), time stamps using add_timestamps(), a progress bar using add_progress() or basic map elements such as a north arrow using add_northarrow() and a scale bar using add_scalebar().

frames <- frames %>% add_labels(title = "White-bearded wildebeest (Connochaetes taurinus) Migration", caption = "Trajectory data: Stabach et al. (2015) 
Map: OpenStreetMap/Stamen; Projection: Geographic, WGS84", 
                                x = "Longitude", y = "Latitude") %>%
  add_timestamps(type = "label") %>% 
  add_progress(colour = "white") %>%
  add_northarrow(colour = "white", position = "bottomleft") %>% 
  add_scalebar(colour = "black", position = "bottomright",
               distance = 10)

frames[[30]]

# Frames can be passed to animate_frames() to turn them into an animation, written as a GIF image or a video file. The user can check available file formats by calling suggest_formats(). The function returns a vector of file suffixes that can be created on the running operating system:

#suggest_formats()
#> [1] "gif"  "mov"  "mp4"  "flv"  "avi"  "mpeg" "3gp"  "ogg"

# In this example, a mov video is created:
#animate_frames(frames, width = 800, height = 800, out_file = "wildebeest_osm.mov", end_pause = 1)

# Okay, this part works fine

# Animating trajectories on a dynamic NDVI base map
# *************************************************
# *************************************************

# Load raster data to use as background
temp <- list.files(path="./Movement_Animation/Data/MODIS", pattern=".tif", recursive = T, full.names=TRUE)
ndvi <- lapply(temp,FUN=stack) 
#ndvi <- lapply(ndvi, FUN = function(x) {x*0.0001}) # These are the real values, but makes file large

# Project to Lat/Long to match movement data
newproj <- "+proj=longlat +datum=WGS84 +no_defs"
ndvi <- lapply(ndvi, FUN = function(x) {projectRaster(x, crs = newproj, method="bilinear")})
#save(ndvi, file = "./Movement_Animation/ndvi.rda")

#load("Data/ndvi.rda")
plot(ndvi[[1]])

# Extract the dates
# *************************
# *************************
#temp <- list.files(path="./Movement_Animation/Data/MODIS", pattern=".tif", recursive = T, full.names=TRUE)
#temp <- extractDate(temp, asDate = TRUE) # Or use the position of the date in the file

#ndvi_times <- temp$inputLayerDates

# Subset times to match the movement dataset
#ndvi_times <- ndvi_times[ndvi_times <= "2010-08-05"]
#ndvi_times <- as.POSIXct(ndvi_times)
#attr(ndvi_times, "tzone") <- "UTC"

#save(ndvi_times, file = "./Movement_Animation/ndvi_time.rda")
load("./Movement_Animation/ndvi_time.rda")

# Both the raster list and the respective dates can be passed to frames_spatial() instead of defining a map service and type. Using fade_raster, the user can decide whether moveVis should continuously fade between raster images by interpolating them over time or instead should switch between raster images discretely. In this example, fading is activated to create a dynamically changing base map.

# Reduce the ndvi list to first 7, starting with 23 April 2010 (Day 113)
ndvi <- ndvi[1:7]

# This part doesn't work..............
frames <- frames_spatial(m, r_list = ndvi, r_times = ndvi_times,
                         fade_raster = TRUE, ext = ext,
                         trace_show = T, trace_colour = "white")

#frames[[200]]

#  Build summary graphs
frames.flow <- frames_graph(m, ndvi, ndvi_times, path_legend = FALSE, graph_type = "flow")
frames.hist <- frames_graph(m, ndvi, ndvi_times, path_legend = FALSE, graph_type = "hist")

# check lengths (must be equal)
sapply(list(frames, frames.flow, frames.hist), length)

# Let's join the graph frames vertically
frames.join.gr <- join_frames(list(frames.flow, frames.hist), ncol = 1, nrow = 2)
frames.join.gr[[200]]
