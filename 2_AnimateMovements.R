## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----Animation Libraries, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------
# Remove everything from memory
rm(list=ls())

# Make sure required libraries are loaded
library(raster)
library(rgdal)
library(move)
library(moveVis)
library(ggplot2)
library(dplyr)
library(lubridate)

# TimeZone and Projection for resulting animation
Timezone1 <- "UTC"
Timezone2 <- "Africa/Nairobi" 
UTM36s.proj <- "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # EPSG:32736
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326


## ----Animation Movement Data---------------------------------------------------------------------------------------------------------------------------
# Read GPS collar data, saved as an object in previous excercises
load("./Data/wildMara.Rdata")

# Subset to the animals to animate
unique(wild.Mara$individual.local.identifier)
data.Mara.An <- subset(wild.Mara, individual.local.identifier == "Naboisho" | individual.local.identifier == "Kayioni" | individual.local.identifier == "Nkairowua" | individual.local.identifier == "Ledama" | individual.local.identifier == "Nkoko" | individual.local.identifier == "Fifteen")
# Or, alternatively (same thing):
#data.Mara.An <- data.Mara %>% filter(individual.local.identifier == "Naboisho" | individual.local.identifier == "Kayioni" | individual.local.identifier == "Nkairowua" | individual.local.identifier == "Ledama" | individual.local.identifier == "Nkoko" | individual.local.identifier == "Fifteen")

# Load shapefile and Reproject to Lat/Long
# **************************************
# **************************************
Mara <- readOGR(dsn = "./Data/SpatialData", layer = "MMNR_UTM36S") 
Mara <- spTransform(Mara, CRS = CRS(LatLong.proj))

# Fortify object - helps to simplify the object for ggplot
mara.DF <- fortify(Mara)


## ----Animation Move Class, message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------
# Convert dataframe to move object
m <- df2move(data.Mara.An, x = "location.long", y = "location.lat", 
             time = "timestamp", 
             track_id = "individual.local.identifier", 
             proj = LatLong.proj,
             removeDuplicatedTimestamps = TRUE)

# Check time stamps and sampling rates (temporal resolution) of the trajectories:
#lag <- unlist(timeLag(m, unit = "mins"))

#summary(lag)
#median(lag)
#sd(lag)


## ----Animation Prepare, warning=FALSE------------------------------------------------------------------------------------------------------------------
# Aligning to a 1 day interval (or whatever makes sense for your data, just be aware of the amount of data included in the animation)
m <- align_move(m, res = 1, digit = 0, unit = "days")
print(paste0("How many days: ",length(unique(timestamps(m)))))

# Subset the dates so runs within a manageable time period.
m <- subset_move(m, from = mdy_hms("05-27-2010 03:00:00"), to = mdy_hms("08-31-2010 12:00:00"))

# check min and max of result
print(paste0("Minimum timestamp: ",min(timestamps(m))))
print(paste0("Maximum timsteamp: ",max(timestamps(m))))

# Create extent for animal (This can take some tweaking, based on the actual animal movements included).
# I'm increasing the extent slightly for visualization purposes
ext <- extent(m) * 1.15
ext@ymin <- ext@ymin * 1.05
ext@ymax <- ext@ymax + abs(ext@ymax * 0.05)


## ----Create Animation, warning=F-----------------------------------------------------------------------------------------------------------------------
# Create a sequence spatial movement maps, customizing the spatial extent, and specifying the map_service and map_type.
frames <- frames_spatial(m, trace_show = TRUE, equidistant = FALSE,
                         ext = ext, map_service = "osm",
                         map_type = "terrain_bg")

# See other basemap options available:
#get_maptypes()

# Add specific details to map, including title, timestamps, progress bar, north arrow, scalebar, and text.
frames <- frames %>% 
  add_labels(title = "White-bearded wildebeest (Connochaetes taurinus)", caption = "Map: OpenStreetMap/Terrain; Projection: Geographic, WGS84", 
                                x = "Longitude", y = "Latitude") %>%
  add_timestamps(type = "label") %>% 
  add_progress(colour = "white") %>%
  add_northarrow(colour = "white", position = "bottomleft") %>% 
  add_scalebar(colour = "black", position = "bottomright",
               distance = 10) %>%
  add_gg(gg=expr(geom_polygon(aes(x=long, y=lat), 
                               size=1, color="goldenrod", data=mara.DF, alpha=0)))

# Add text (Piping doesn't seem to work here)
frames <- add_text(frames, "Loita Plains", x = 35.583, y = -1.296,
                     colour = "black", size = 4) 
frames <- add_text(frames, "Mara Plains", x = 35.2, y = -1.21,
                   colour = "black", size = 4)
frames <- add_text(frames, "Maasai Mara National Reserve", x = 35.228, y = -1.58,
                   colour = "black", size = 4)

# Check to make sure everything is looking okay (important step because the animation processing time can be lengthy):
frames[[30]]


## ----Animation Export, warning=FALSE-------------------------------------------------------------------------------------------------------------------
# Additional file formats:
#suggest_formats()

# Create animation
animate_frames(frames, width = 1000, height = 1000, out_file = "Wildebeest_Example.mov", end_pause = 1, overwrite = TRUE)

