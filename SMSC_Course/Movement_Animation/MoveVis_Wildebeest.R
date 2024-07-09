library(moveVis)
library(move)

# Load Data
load("Move_Data_Example.Rda")

# Align movement
m <- align_move(m, res = 180, unit = "mins") # Align to 3 hours

# Create spatial frames
frames.sp <- frames_spatial(m, r_list = ndvi, r_times = ndvi_times, r_type = "gradient",
                            fade_raster = TRUE)
#frames.sp <- add_colourscale(frames.sp, type = "gradient",
#                             colours = c("orange", "white", "darkgreen"), legend_title = "NDVI")
frames.flow <- frames_graph(m, ndvi, ndvi_times, path_legend = FALSE, graph_type = "flow")
frames.hist <- frames_graph(m, ndvi, ndvi_times, path_legend = FALSE, graph_type = "hist")

# check lengths (must be equal)
sapply(list(frames.flow, frames.hist), length)

# Let's join the graph frames vertically
frames.join.gr <- join_frames(list(frames.flow, frames.hist), ncol = 1, nrow = 2)
frames.join.gr[[200]]
