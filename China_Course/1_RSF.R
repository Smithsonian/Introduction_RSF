## ----setup, include=FALSE----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----RSF Libraries, message=FALSE, warning=FALSE-----------------------------
# Remove everything from memory
rm(list=ls())

# Set.seed - This function means that any generation of "random" points/numbers will be the same each time you run the script.
set.seed(533)

# You may need to install these packages first
#install.packages('amt', 'tidyverse', 'lme4', 'terra', 'sf', 'sjPlot', 'visreg', 'usdm', 'tmap', 'jtools')

# Load libraries
library(amt)
library(tidyverse)
library(lme4)
library(terra)
library(sf)
library(sjPlot) 
library(visreg) 
library(usdm)
library(tmap)
library(jtools)

# Set all tmaps to plot in view mode
tmap_mode("view")


## ----Load, message=FALSE, warning=FALSE, echo=TRUE---------------------------
# Load a polygon shapefile of the study area.  File is project to Albers Equal Area.  This will be our prediction area.
Athi.Bound <- st_read("Data/Athi.shp")

# Load all raster layers for Athi-Kaputiei Plains study area.  The raster stack consists of 7 raster layers.  To be in a stack they all have to have the exact same resolution (250-m) and spatial extent. 
rsf.stack <- terra::rast("data/ak_raster_stack.tif")
plot(rsf.stack)

# Data Layers included:
# anth_risk - Anthropogenic Risk, simply an index of human footprint made specifically for this ecosystem (expected negative response)
# Fence_dist - Fence Distance, with fences manually mapped by a Kenyan field team (expected negative response)
# prirds_dist - Primary Road Distance, the distance from primary/paved roads (expected negative response)
# secrds_dist - Secondary Road Distance, the distance from secondary/unpaved roads (expected null response)
# river_dist - River Distance, the distance from permanent rivers (depending on season, but wildebeest must drink - greater attraction in dry season, but a predation risk)
# waterpts_dist - Water Point Distance, the distance to mapped water wells (expected positive response, stronger in dry season)
# woody_dist -Woody Distance, the distance to woody vegetation (based on VCF) (expected negative response - predation risk)

# The fence plot is hard to interpret because there are so many low values in the center of the plot. This means there are a lot of fences there. Let's color the distance values and only show the first 100 to better see where the fences are located.
# plot(rsf.stack$fence_dist)
plot(rsf.stack$fence_dist,
            range = c(0, 100),
            col = "black")

# We can review the range of values in each raster like this.
# summary(values(rsf.stack$fence_dist))

# What is the crs of the raster stack?
crs(rsf.stack, proj=TRUE)

# *******************************
# Now let's import the GPS tracking dataset.  
# This is the 3-hour, regular trajectory.
# Let's import and remove the trajectory information
# We also need to remove the na values in the x/y coordinates
WB.data <- read_rds("Data/wildebeest_3hr_adehabitat.rds") %>% 
  select(x,
         y,
         date,
         id,
         sex) %>% 
  filter(!is.na(x),
         !is.na(y),
         !is.na(date))

# Is this a spatial object?
#class(WB.data)
# Is it projected?
#st_crs(WB.data)
# What is the timezone?
#tz(WB.data$date)

# Here's where we need to be careful.
# Our raster layers are projected to Albers Equal Area.  Our point layer is not projected, but the x/y coordinates are derived from when the layer was projected to UTM 37S, WGS84.
# This means we should sent the file to UTM37S, WGS84 and project to Albers Equal Area
WB.data.sf <- WB.data %>%  
  st_as_sf(coords = c('x', 'y'),
           crs = "EPSG:32737") %>% # This is UTM 37S, WGS84 (UtmZone.proj <- "EPSG:32737")
  st_transform(crs(rsf.stack)) # Easiest way to set the projection is grab it from the crs of our already defined raster stack

# The specific project is Albers Equal Area, specific to Africa.  This would be defined as:
#AEA.Africa.proj <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Now what's the CRS?
# st_crs(WB.data.sf)

# Let's plot to make sure everything is aligned.
# Here I'm choosing to plot with a sequential palette with 10 colors on the ramp. You can run tmaptools::palette_explorer() to explore the palette options. 
tm_shape(rsf.stack$anth_risk,
         name = "Human Risk Level") +
  tm_raster(palette = "YlOrRd", n = 10,
            alpha = 0.6,
            title = "Anth. Risk") +
  tm_shape(WB.data.sf %>% 
             filter(id == "Noontare"),
           name = "Noontare Locations") +
  tm_dots(size = 0.01,
          col = "black") + 
  tm_shape(Athi.Bound,
           name = "Athi-Kaputiei Plains") +
  tm_polygons(alpha = 0,
              border.col = "green") + 
  tm_layout("Example Graph")

# This first raster represents anthropogenic risk. And you can see that this wildebeest appears to be responding to this variable. If you zoom into the northern locations, as there are no locations in the areas of elevated human risk.

# Let's plot one more. Let's look at this wildebeest's points with distance to primary road. Here I want the low values to be red, instead of the high values, so I'm reversing the palette with the "-". I'm also including a lot more colors, otherwise it was difficult to see where the roads actually were. This makes the legend too large though so I'm hiding the legend for the raster.

# tm_shape(rsf.stack$prirds_dist,
#          name = "Distance to Primary Rd") +
#   tm_raster(palette = "-YlOrRd", n = 30,
#             alpha = 0.6,
#             legend.show = FALSE) +
#   tm_shape(WB.data.sf %>%
#              filter(id == "Noontare"),
#            name = "Noontare Locations") +
#   tm_dots(size = 0.01,
#           col = "gray")

# Check tmap_options()
# This is where the default basemaps are set.  This is why we see ESRI.WorldGrayCanvas, OpenStreetMap, and ESRI.WorldTopoMap


## ----Create Subset, message=FALSE, warning=FALSE, echo=TRUE------------------
# Create Subset of dataset, filtering to 1 animal
nt <- WB.data.sf %>% 
  mutate(x = st_coordinates(WB.data.sf)[ ,1],
         y = st_coordinates(WB.data.sf)[ ,2]) %>% 
  as_tibble() %>% 
  select(-geometry) %>%  
  filter(id == "Ntishya") %>% 
  mutate(id = droplevels(id))
  
# Create Movement Track (similar to ltraj function in adehabitatLT)
# We don't really need the time component for this initial analysis.  That would be a step-selection function (ssf) analysis.  But, we'll include here as an example.
nt.trk <- mk_track(nt, 
                   .x = x, 
                   .y = y,
                   .t = date, 
                   crs = crs(WB.data.sf), # Alternative just include crs = 9822.  9822 is the EPSG code for AEA.
                   order_by_ts = T,
                   id = id, 
                   sex = sex)

# Note that amt does not calculate steplengths or turning angles when you make the track.  To do so, you'd need to specify those variables:
# nt.trk.test <- nt.trk %>% 
#   mutate(
#     sl = step_lengths(.),
#     ta = direction_abs(.) # absolute and relative turning angles can be calculated (direction_rel())
#     )

# Remember also that the track should be kept regular, otherwise we have potential to bias our results (e.g., more points in the day than at night). Our track has already been resampled to a 3 hour interval in a previous lecture
# In AMT, we would do this by using the track_resample command:
nt.trk.rs <- track_resample(nt.trk,
                            rate = hours(3), 
                            tolerance = minutes(20))

# Let's thin the data as an example   
rcd.amt <- ceiling(nrow(nt.trk.rs)*0.20) # ceiling just rounds up to a whole number
nt.trk.rs <- nt.trk.rs[sample(nrow(nt.trk.rs),size = rcd.amt),c("id","x_","y_","id","sex")]


## ----Availability, message=FALSE, warning=FALSE, echo=TRUE-------------------
# Here we will generate 10 times the number of "Use" points to get our list of "available" points. We'll start with this level to begin to inspect our data and look for initial patterns. Before running our final models we'll assess what a robust minimum # of available points should be.

# Create simple boundary of points
hr <- hr_mcp(nt.trk.rs, levels = 1) # Amt has multiple options, including hr_akde().  Levels indicates the isopleth (Here: 100%).  We could also upload our own polygon.

# Generate random points within the defined polygon.  Important to include the type (random or regular) that we want to be generated.
nt.rsf.10 <- random_points(hr,
                           n = nrow(nt.trk.rs) * 10,
                           type ="regular",
                           presence = nt.trk.rs) 

# What does our dataset look like now?
head(nt.rsf.10)
# Cool, the function created a 'case_' field, tracking our 'Use' points (TRUE) and our 'Availability' points (FALSE)

# Summarize
table(nt.rsf.10$case_)
#class(nt.rsf.10)

# Plot Use/Availability
plot(nt.rsf.10) 


## ----Extract, message=FALSE, warning=FALSE, echo=TRUE------------------------
# Extract all the raster variables at Use/Available points.
nt.rsf.10 <- nt.rsf.10 %>% 
  extract_covariates(rsf.stack)

# Look at result
head(nt.rsf.10)
#names(nt.rsf.10)

# Now let's summarize these results to get a sense of how the values compare between 'use' and 'available' locations
nt.rsf.10 %>% 
  ggplot(aes(y = anth_risk,
             col = case_)) + 
  geom_boxplot() +
  labs(title = "Anthropogenic Risk")

# You may have many covariates, so doing this in a loop makes sense. Instead using a 'for' loop, we'll use the 'map' function.  See 'help(map)'.  First we need to make a vector of variable names that we want to plot.
vars <- nt.rsf.10 %>% 
  select(anth_risk:woody_dist) %>% 
  names()

# The, we simply use 'map' to create a plot for each var.  the '.x' is like our i in our 'for' loop.
box.plots <- map(vars, ~
                   nt.rsf.10 %>% 
                   select(case_, var = .x) %>% 
                   ggplot(aes(y = var,
                              col = case_)) + 
                   geom_boxplot() + 
                   labs(title = .x))

# We can put them all together using the cowplot
cowplot::plot_grid(plotlist = box.plots)


## ----Collinearity, message=FALSE, warning=FALSE, echo=TRUE-------------------
# Assess collinearity
cor(nt.rsf.10[,4:10])

# The table of correlation values indicates that high levels of correlation exist between fence distance and waterpoint distance (0.84).  We should remove one of these variables based on our research objectives.  Here, I'll keep fence distance because I am more interested in this question from a management standpoint.  Some other correlations are observed between woody distance and waterpoints and woody distance and primary roads.

# Before making any decisions, let's also check the Variance inflation Factor. 
vif(as.data.frame(nt.rsf.10[,4:10]))

# Let's remove waterpoints and see if doing so is helpful
vif(as.data.frame(nt.rsf.10[,c(4:8,10)]))

# This helps a alot, but woody distance is still quite high (VIF  = 5.7).  We could make a decision to remove entirely or if we think it could be an important factor (e.g., woody vegetation is likely avoided by wildebeest because of increased risk of predation), we could include it while making sure woody vegetation was not included in the same model with fences and primary roads.  We could then evaluate each model separately using AIC.  For now, we'll proceed.


## ----Sensitivity, message=FALSE, warning=FALSE, echo=TRUE--------------------
# Setup number of available locations to sample during each model fitting
n.frac <- c(1, 5, 20, 50, 100) 

# Total available locations to be generated, based on the number of Use points
n.pts <- nrow(nt.trk.rs) * n.frac

# Number of repetitions.  This is simulation so that we can evaluate the variability that exists in the coefficients
# For a publication, I would increase the n.rep to 100.  Here, 20 is fine.
n.rep <- 20

# Create a table which saves the settings of each scenario
# We then extract the covariates during each repetition from rsf.stack during each run (the points will vary)
# Then, fit a glm model for each available location (1,5,20,50,100) and for each replication (re-sampling)

# Run simulation and store results
# Commented out as the process is time consuming
# **********************************************
# wb.sim <- tibble(
#  n.pts = rep(n.pts, n.rep),
#  frac = rep(n.frac, n.rep),
#  result = map(
#     n.pts, ~
#       nt.trk.rs %>% random_points(n = .x) %>%
#       extract_covariates(rsf.stack) %>%
#       mutate(anth_risk = scale(anth_risk),
#              woody_dist = scale(woody_dist),
#              fence_dist = scale(fence_dist),
#              prirds_dist = scale(prirds_dist),
#              river_dist = scale(river_dist),
#              secrds_dist = scale(secrds_dist),
#              waterpts_dist = scale(waterpts_dist)) %>%
#       # Fit basic model
#       glm(case_ ~ anth_risk +
#             woody_dist +
#             fence_dist +
#             prirds_dist +
#             river_dist +
#             secrds_dist +
#             waterpts_dist,
#           data = ., family = binomial(link = "logit")) %>%
#       broom::tidy()))

# Save file so don't need to run everytime
# write_rds(wb.sim, file = "Data/nt.sim.rds")

# Read in the result
wb.sim <- read_rds("Data/nt.sim.rds")

# Look at the summary table of results
# This shows that we have 100 rows, summarizing the number of points (n.pts) and the availability fraction (frac)
# We have 100 rows because we have 5 different fractions (n.rep = 5 -> 1,5,20,50,100) and we ran the simulations (n.rep = 20) times.  5 x 20 = 100 rows of results
wb.sim 
# Look at the first simulation result
wb.sim$result[[1]]

# Visualize/Graph findings
# We must 'unnest' the contents in each nested dataframe so visualize the coefficient estimates from individual model fits
wb.sim %>% unnest(cols = result) %>% 
  mutate(term = recode(term, 
                       "(Intercept)" = "Intercept",
                       anth_risk = "Anthropogenic disturbance", 
                       woody_dist = "Woody Distance",
                       fence_dist = "Fence Distance", 
                       prirds_dist = "Pri Rds Distance",
                       river_dist = "River Distance",
                       secrds_dist = "Sec Rds Distance",
                       waterpts_dist = "Water Point Distance")) %>% 
  ggplot(aes(factor(frac), 
             y = estimate)) +
  geom_boxplot() + 
  facet_wrap(~ term, scale  ="free") +
  geom_jitter(alpha = 0.2) + 
  labs(x = "Available Points per Use Location", 
       y = "Estimate") +
  theme_light()

# This is a key tool in making a final assessment of the minimum # of available points needed for each use location. We can see little change in the coefficient values once we have 20 points per use location, but we'll use 50 points per 'Use' location to be safe.


## ----Model Fitting, message=FALSE, warning=FALSE, echo=TRUE------------------
# Generating 50 available locations per each used location within animal homerange
nt.rsf.50 <- random_points(hr,
                           n = nrow(nt.trk.rs) * 50,
                           type ="regular",  
                           presence = nt.trk.rs) 

# Extracting the covariates
nt.rsf.50 <- nt.rsf.50 %>% extract_covariates(rsf.stack)

# Fit a "full" model after removing waterpoints because of high collinearity. We also should avoid running fence distance or primary roads with woody distance, because they are highly correlated. As a result, let's fit three "full" models and then compare them.

# Fitting Model 1 - No woody distance
M1.NoWood <- glm(case_ ~ scale(anth_risk) + 
                   scale(fence_dist) + 
                   scale(prirds_dist) + 
                   scale(secrds_dist) + 
                   scale(river_dist), 
                 family = binomial(link="logit"), 
                 data = nt.rsf.50)

# Get Summary & calculate profile confidence intervals to see if coefficients overlap zero
summary(M1.NoWood)
# broom::tidy(M2.Wood) # Could also use the tidy summary, which is a little easier to read
confint(M1.NoWood)

# Fitting Model 2 - No Fencing or Primary Roads
M2.Wood <- glm(case_ ~ scale(anth_risk) + 
                 scale(secrds_dist) + 
                 scale(river_dist) +
                 scale(woody_dist), 
               family = binomial(link="logit"), 
               data = nt.rsf.50)

# Get Summary & calculate profile confidence intervals to see if coefficients overlap zero
summary(M2.Wood)
confint(M2.Wood)

# Compare models with AIC
AIC(M1.NoWood, M2.Wood)


## ----Model Interpretation, message=FALSE, warning=FALSE, echo=TRUE-----------
# Plot raw, non-transformed estimates (log-odds)
plot_model(M1.NoWood, transform = NULL) 

# This indicates that we have two predictors (river distance and secondary roads) that have positive effects and three variables (anthropogenic risk, fence distance, and primary roads) that have negative effects.  

# Specifically, we see that when we increase human risk (1 unit increase), we reduce (negative sign) the relative use probability.  This is what we would expect (i.e., more human population leads to less wildebeest).  For fences and primary roads, we also see a negative effect.  Here, a 1 unit increase in the distance to fences or primary roads results in a negative effect on the relative use probability of wildebeest.  That is, relative use probability of wildebeest is higher when closer to these features (a positive association), than further away.  This is an unexpected result.  We should remember, however, that we are only investigating 1 animal at this time.

# For secondary roads and rivers, we see positive effects (i.e., an increase in the distance to roads and rivers leads to an increase in the relative use probability of wildebeest).  That is, relative use probability is higher as we move away from these features (a negative association).  This is also unexpected, at least for secondary roads.

# Create an alternative plot of the raw coeeficients (log-odds)
coef_plot <- sjPlot::plot_model(M1.NoWood,
                   type = "est",
                   colors = c("cadetblue", "tomato"),
           transform = NULL,
           title = "Plot of Standardized Coefficents") 
# coef_plot

# We can keep the above plot, or improve it using ggplot.  Here we'll add a horizontal line, centered on 0, so that we can easily see the positive and negative effects.
coef_plot +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             col = "darkgray") +
  theme_classic() 

# Check out help(plot_model), as there are a number of customizable arguments.  You could, for example, use the "terms" argument to indicate which terms you want to show on the plot and also show the odds ratios instead of the log odds (raw coeffients is the defaul).  The odds ratios, however, are harder to interpret when the coefficients are negative, as we'll see below. 
# coef_plot.example <- sjPlot::plot_model(M1.NoWood,
#                    type = "est",
#                    colors = c("cadetblue", "tomato"),
#                    terms = c("scale(anth_risk)","scale(fence_dist)"), # Will only plot the predictors specified
#            transform = "exp", # Exponentiating to plot the odds ratios
#            title = "Plot of Standardized Coefficents") 
# 
# coef_plot.example

# Now let's convert the log odds to odds ratios to evaluate Relative Selection Strength.  We do this by exponentiating the coefficients
my_coefs <- coef(M1.NoWood)
exp(my_coefs)

# Remember: Don't interpret the intercept in RSF models.  This is meaningless.
# Any value below 1 has a negative impact on use probability; Any value above 1 is a positive impact. These values can be interpreted as the relative change in the probability of use with 1 unit of change in the predictor.  Easy to interpret for predictors with a positive sign.

# For negative coefficients, it's easier to include a negative sign when exponentiating
exp(-my_coefs[2:4])

# Thus, for example, we can say a location with a 1 unit increase in anth_risk is 1.20 times LESS likely to be used by this wildebeest.  We see the same general patterns for fence distance (1.72 times LESS likely) and primary road distance (2.55 times LESS likely).

# One last IMPORTANT note: We scaled our predictor variables.  As a result, a 1 unit change is actually a change in 1 STANDARD DEVIATION.  If we did not scale our predictor variables, a 1 unit change in the distance variables would be 1 meter, for example. Just something to be aware of. 


## ----Visualize, message=FALSE, warning=FALSE, echo=TRUE----------------------
# Use visreg to visualize
# visreg(M1.NoWood,
#        xvar = "anth_risk",
#        scale="response", 
#        ylab="Relative Use Intensity",
#        xlab="Anthropogenic risk",
#        partial=F,
#        rug=F,
#        line=list(col="black"), 
#        fill=list(col="light gray"))
# 
# # Let's loop over all the variables and plot them
# # Step 1: Create an object holding the variables
# var1 <- c("anth_risk",
#           "fence_dist",
#           "prirds_dist",
#           "secrds_dist",
#           "river_dist")
# 
# # Step 2: Create an object of the variable names (for plotting labels)
# var.names <- c("Anthropogenic Risk Index",
#                "Fence Distance (m)",
#                "Primary Road Distance (m)",
#                "Secondary Roads Distance (m)",
#                "River Distance (m)")
# 
# # Step 3: Set plotting window
# par(mfrow=c(3,2))
# 
# # Step 4: Loop over each variable and plot
# for (i in 1:length(var1)){
#   visreg(M1.NoWood,
#          xvar = var1[i],
#          scale="response", 
#          ylab="Relative Use Intensity",
#          xlab=var.names[i],
#          partial=F,
#          rug=F,
#          line=list(col="black"), 
#          fill=list(col="light gray"))
# }
# 
# # Step 5: Reset plotting window
# par(mfrow=c(1,1))

# Use log_rss()
# Step 1: Create a dataframe of a sequency of values to predict (seq(min to max)).  Note that I am holding all other variables constant (mean values)
df1 <- data.frame(anth_risk = seq(min(nt.rsf.50$anth_risk),
                                    max(nt.rsf.50$anth_risk),
                                    length.out = 100),
                    fence_dist = mean(nt.rsf.50$fence_dist),
                    river_dist = mean(nt.rsf.50$river_dist),
                    prirds_dist = mean(nt.rsf.50$prirds_dist),
                    secrds_dist = mean(nt.rsf.50$secrds_dist))

# Step 2: Create a dataframe of what are comparing to.  Here, I'm comparing to the mean anth_risk value
df2 <- data.frame(anth_risk = mean(nt.rsf.50$anth_risk),
                    fence_dist = mean(nt.rsf.50$fence_dist),
                    river_dist = mean(nt.rsf.50$river_dist),
                    prirds_dist = mean(nt.rsf.50$prirds_dist),
                    secrds_dist = mean(nt.rsf.50$secrds_dist))

# Step 3: Use the log_rss() function to predict across our sequence using coefficients from our model
logRSS_riskRange <- log_rss(object = M1.NoWood,
                            x1 = df1,
                            x2 = df2)

# Look at result
# logRSS_riskRange

# This output is set up nicely for plotting, since we have our input values of risk, and we have our output values for log-RSS. I'd like to plot the odds-ratio values though (as opposed to the log-odds) so in plotting, we'll take the exp() of those values.

# Step 4: Plot
ggplot(logRSS_riskRange$df, 
       aes(x = anth_risk_x1, 
           y = exp(log_rss))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = exp(0), 
             linetype = "dashed", 
             color = "gray30") +
  xlab("Anthropogenic Risk") +
  ylab("RSS vs Mean Risk") +
  theme_bw()
  
  # The RSS of 1.0 here crosses the line at the mean value of risk, since that is what this relative value (odds ratio) is being compared to. 


## ----Prediction, message=FALSE, warning=FALSE, echo=TRUE---------------------
# Extract the coefficient of each predictor from the model summary
coeff <- M1.NoWood$coefficients

# Then, use the logistic equation to generate predictions (no B_0)
# w*(x)=exp(β1x1 + β2x2 +.... + βnxn)/exp(1 + β1x1 + β2x2 +.... + βnxn)

# where w*(x) is the relative probability of selection, dependent upon covariates X1 through Xn, and their estimated regression coefficients β1 to βn, respectively.

# IMPORTANT: Since we scaled our parameters, we MUST also scale our raster layers.  The scaled raster layers will look the same, but the range of values being plotted will change.  Here, I'm scaling by subtracting the mean value of the parameter and dividing by the standard deviation.  This is what we did above by using the scale() function above.  See help(scale).

# See how this is written: ScaledRasterLayer = (unscaled raster layer - mean(dataset$predictor)) / sd(dataset$predictor)
anth.scale <- (rsf.stack$anth_risk - mean(nt.rsf.50$anth_risk)) / sd(nt.rsf.50$anth_risk)
fence.scale <- (rsf.stack$fence_dist - mean(nt.rsf.50$fence_dist)) / sd(nt.rsf.50$fence_dist)
prirds.scale <- (rsf.stack$prirds_dist - mean(nt.rsf.50$prirds_dist)) / sd(nt.rsf.50$prirds_dist)
secrds.scale <- (rsf.stack$secrds_dist - mean(nt.rsf.50$secrds_dist)) / sd(nt.rsf.50$secrds_dist)
river.scale <- (rsf.stack$river_dist - mean(nt.rsf.50$river_dist)) / sd(nt.rsf.50$river_dist)

# Prediction - Need to get the coefficient numbers [[i]] correct
pred <- exp(anth.scale*coeff[[2]] + 
              fence.scale*coeff[[3]] + 
              prirds.scale*coeff[[4]] + 
              secrds.scale*coeff[[5]] + 
              river.scale*coeff[[6]]) /
  (1+exp(anth.scale*coeff[[2]] + 
           fence.scale*coeff[[3]] + 
           prirds.scale*coeff[[4]] + 
           secrds.scale*coeff[[5]] + 
           river.scale*coeff[[6]]))

# Provide Spatial Prediction - Based off of the coefficients from this single animal
plot(pred)

# Let's mask the data with the Athi-Kaputiei Boundary file and then use tmap to plot
# We loaded a spatial polygon data layer at the start of this script named "Athi.Bound"
pred.nt <- mask(pred, Athi.Bound)

# Graph on tmap
tm_shape(pred.nt,
         name = "Habitat Suitability") +
  tm_raster(palette = "PuBuGn", n = 10,
            alpha = 0.6,
            title = "Habitat Suitability") +
  tm_shape(Athi.Bound,
           name = "Athi-Kaputiei Plains") +
  tm_polygons(alpha = 0, # Make polygon transparent
              border.col = "green") + 
  tm_layout("Ntishya Example")

# Save this file in your data directory
# writeRaster(pred.nt, filename ="./Output/Prediction.nt.tif", overwrite = TRUE)


## ----GLMER, message=FALSE, warning=FALSE, echo=TRUE--------------------------
# Grab the original dataset (before filtering for Ntishya)
nt <- WB.data.sf %>% 
  mutate(x = st_coordinates(WB.data.sf)[ ,1],
         y = st_coordinates(WB.data.sf)[ ,2]) %>% 
  as_tibble() %>% 
  select(-geometry)

# Grab the Ids from the file
WB.Id <- unique(nt$id)

# Creating a list that contain the dataframe for each individual
# Splitting by id
nt <- split(nt,nt$id)

# Creating track for each individual animal in the list
# Note, we've again included time, but we're not using it in the model structure
WB.trk.all <- lapply(nt, function(x) mk_track(x,
         .x=x,
         .y=y,
         .t=date,
         crs = crs(rsf.stack),
         order_by_ts = T,
         id = id,
         sex = sex))

# Let's reduce the file of each id subset, selecting only 20% of each
WB.trk.all <- lapply(WB.trk.all, 
                     function(x) x[sample(1:nrow(x), 
                                          size = ceiling(nrow(x)*0.20)),c("x_","y_","t_","id","sex")])

# Generate 50 random points for each animal 
# Using a list apply here to apply the function
# This will create the case_ column
WB.Athi <- lapply(WB.trk.all, 
                  function(x) random_points(x, 
                                            n=nrow(x) * 50,
                                            typ="regular"))


## ----GLMER Extract, message=FALSE, warning=FALSE, echo=TRUE------------------
# Extract all covariates
WB.Athi <- lapply(WB.Athi, function(x) extract_covariates(x, rsf.stack))

# Note, we lost the ID column, so need to put that back in (that's why we created above).  We'll use this as the random effect.
WB.Athi <- mapply(cbind, WB.Athi, "Animal_ID" = WB.Id, SIMPLIFY=F) 

# Binding all dataframes in the list into a single dataframe for analysis
WB.Athi <- do.call(rbind,WB.Athi)


## ----GLMER Fit, message=FALSE, warning=FALSE, echo=TRUE----------------------
# This model will be slow to execute
# mixed.model <- glmer(case_ ~ scale(anth_risk) + 
#                         scale(fence_dist) +
#                         scale(prirds_dist) + 
#                         scale(secrds_dist) +
#                         scale(river_dist) +
#                         (1 + scale(anth_risk) + 
#                            scale(fence_dist) + 
#                            scale(prirds_dist) +
#                            scale(secrds_dist) +
#                            scale(river_dist)
#                          |Animal_ID), 
#                       family = binomial(link = "logit"), 
#                       data = WB.Athi)

# Save/Load the model
# saveRDS(mixed.model, file = "Output/MixedModel.rds")
# Data load/import
mixed.model <- read_rds("Output/MixedModel.rds")

# Print Model Summary
summary(mixed.model)

# Look at the random effects
ranef(mixed.model)


## ----GLMER Response, message=FALSE, warning=FALSE, echo=TRUE-----------------
# Plot coefficients and response curves for inference at the population level
plot_model(mixed.model,transform = NULL) 

# We can see that the coefficient plot works the same, and it just shows the fixed effects. But we can see the huge impact it has on the variance estimation if we allow for random effects.

# Exponentiate the coefficients and plot the odds ratios
#tab_model(mixed.model)

# Graph response curves.  Use effect plot from jtools.  We could do this for all our parameters.
effect_plot(mixed.model,
            data = WB.Athi, 
            pred = anth_risk,
            interval = TRUE,
            x.label = "Anthropogenic risk",
            y.label = "Relative Selection Strength") + 
  theme_classic()


## ----GlMER Prediction, message=FALSE, warning=FALSE, echo=TRUE---------------
# Extract the coefficient of each predictor from the model summary
coef <- summary(mixed.model)
coeff <- coef$coefficients

# Scale rasters.  Need to do this again because the sample used is different than previously calculated:
anth.scale <- (rsf.stack$anth_risk - mean(WB.Athi$anth_risk)) / sd(WB.Athi$anth_risk)
fence.scale <- (rsf.stack$fence_dist - mean(WB.Athi$fence_dist)) / sd(WB.Athi$fence_dist)
prirds.scale <- (rsf.stack$prirds_dist - mean(WB.Athi$prirds_dist)) / sd(WB.Athi$prirds_dist)
secrds.scale <- (rsf.stack$secrds_dist - mean(WB.Athi$secrds_dist)) / sd(WB.Athi$secrds_dist)
river.scale <- (rsf.stack$river_dist - mean(WB.Athi$river_dist)) / sd(WB.Athi$river_dist)

# Prediction
pred.mixed <- exp(anth.scale*coeff[[2]] + 
                    fence.scale*coeff[[3]] + 
                    prirds.scale*coeff[[4]] + 
                    secrds.scale*coeff[[5]] + 
                    river.scale*coeff[[6]]) /
  (1+exp(anth.scale*coeff[[2]] + 
           fence.scale*coeff[[3]] + 
           prirds.scale*coeff[[4]] + 
           secrds.scale*coeff[[5]] + 
           river.scale*coeff[[6]]))

# Provide Spatial Prediction - Based off of the coefficients from this single animal
plot(pred.mixed)

# Mask result and plot using tmap with the Athi Bounday
pred.mixed <- mask(pred.mixed, Athi.Bound)

# Graph on tmap
tm_shape(pred.mixed,
         name = "Habitat Suitability") +
  tm_raster(palette = "PuBuGn", n = 10,
            alpha = 0.6,
            title = "Habitat Suitability") +
  tm_shape(Athi.Bound,
           name = "Athi-Kaputiei Plains") +
  tm_polygons(alpha = 0, # Make polygon transparent
              border.col = "green") + 
  tm_layout("Random Effects Example")

# Save this file in your data directory
writeRaster(pred.nt, filename ="Output/Prediction.mixed.tif", overwrite = TRUE)

