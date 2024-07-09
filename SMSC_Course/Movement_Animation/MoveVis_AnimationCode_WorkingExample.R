library(moveVis)
library(move)

data("move_data", "basemap_data")

# align movement
m <- align_move(move_data, res = 4, unit = "mins")

# create spatial frames and graph frames:
r_list <- basemap_data[[1]]
r_times <- basemap_data[[2]]

frames.sp <- frames_spatial(m, r_list = r_list, r_times = r_times, r_type = "gradient",
                            fade_raster = TRUE)
frames.sp <- add_colourscale(frames.sp, type = "gradient",
                             colours = c("orange", "white", "darkgreen"), legend_title = "NDVI")
frames.flow <- frames_graph(m, r_list, r_times, path_legend = FALSE, graph_type = "flow")
frames.hist <- frames_graph(m, r_list, r_times, path_legend = FALSE, graph_type = "hist")

# check lengths (must be equal)
sapply(list(frames.sp, frames.flow, frames.hist), length)

# Let's join the graph frames vertically
frames.join.gr <- join_frames(list(frames.flow, frames.hist), ncol = 1, nrow = 2)
frames.join.gr[[100]]

# Now, let's join the joined graph frames with the spatial frames horizontally
# in 2:1 ration and align all axis
frames.join <- join_frames(list(frames.sp, frames.join.gr),
                           ncol = 2, nrow = 1, rel_widths = c(2, 1), axis = "tb")
frames.join[[100]]
# in a standard graphics device, this looks a bit unproportional
# however when setting the correct width, height and resolution of a graphic device,
# it will come out well aligend.

# Do so for example with animate_move() with width = 900, dheight = 500 and res = 90
animate_frames(frames.join, out_file = tempfile(fileext = ".gif"), fps = 25, 
               width = 900, height = 500, res = 90, display = TRUE, overwrite = TRUE)

## End(Not run)