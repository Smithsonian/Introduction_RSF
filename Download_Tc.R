library(gfcanalysis)
library(tidyverse)
library(sf)

Athi.Bound <- st_read("China_Course/Data/Athi.shp")
# Put in Geographic
Athi.LL <- Athi.Bound %>% 
  st_transform(crs("EPSG:4326"))

# Define Area
tc.aoi <- calc_gfc_tiles(Athi.LL)

# Pull data
tc <- download_tiles(tc.aoi, output = "China_Course/Output", images = "treecover2000", "GFC-2022-v1.10")

rast <- rast("China_Course/Data/ak_raster_stack.tif"
             )
