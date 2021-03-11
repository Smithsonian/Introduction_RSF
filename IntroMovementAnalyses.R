# Project: Wildebeest
# Description: Introduction to Animal Movement Analyses using wildebeest dataset as a template  
#               There's lots we could do here.  I'd like to start with meeting class objectives, but would eventually be interested in simulating animal movement based on iSSF model (new to me)
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

# I am going to put all of this in a Markdown document as it develops and will credit you for your assistance

# Remove anything in memory
rm(list=ls())

# Load Libraries
library(ctmm)
library(lubridate)
library(dplyr)
library(raster)

# Load directly if dataset is ready from Movebank
# Ex. yourAnimals <- as.telemetry("yourAnimalsMoveBank.csv")

# Read in file downloaded from Movebank to correct column issues (timestamp)
wild <- read.csv("./Data/White-bearded wildebeest in Kenya.csv", header=TRUE)
head(wild)

# Must fix the time, first by correcting the time in the Excel file
# I then needed to set the format of the time to POSIXct using lubridate
# Find the function that matches the format of your data
wild$timestamp <- mdy_hm(wild$study.local.timestamp)
#wild  <- wild %>% mutate(timestamp = mdy_hm(study.local.timestamp)) # Could use mutate to get the same result
head(wild)
str(wild)

# Create ID field and Re-arrange Columns
wild <- wild[,c(3:9,11,13:14,16:19)]
wild <- wild %>% 
  mutate(id = tag.local.identifier) %>% 
  relocate(id)

# Let's restrict to a few animals in the Maasai Mara.  Animals that I know that fit the home range assumption.
Mara.An <- c(2829, 2830, 2831, 2832, 2833, 2834, 2835, 2836, 2838, 2839, 2841, 2843, 2844, 2845, 2846)
#Ambo.An <- c(2840, 2842, 30068, 30070, 30071, 30072, 30074, 30077, 30079, 30082, 30084, 30086)
#Athi.An <- c(2837, 30069, 30073, 30075, 30076, 30078, 30081, 30083, 30085)

# Subset to Mara animals
wild2 <- wild[(wild$id %in% Mara.An),]

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
