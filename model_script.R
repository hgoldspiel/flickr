#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R CODE FOR FLICKR-NFR POINT PATTERN ANALYSES
# Harrison B Goldspiel | hbgoldspiel@gmail.com
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load("data/cluster_analysis_output.RData")

# Predictive models for NFR Flickr data
## 1. Link images to geographic raster data
## 2. Generate random background points
## 3. Create predictive models for focal areas (NFR, states, main parks)
## 4. Identify relevant geographic drivers of visitation
## 5. Use models to predict CES suitability surfaces for focal areas

# LOAD PACKAGES
library(beepr)
library(sf)
library(raster)
library(ranger)
library(randomForest)
library(Boruta)



# LOAD GIS DATA
rastlist <- list.files(path = "data/GIS/baselayers", patter = '.tif$',
                       all.files = TRUE, full.names = TRUE)
rastdata <- stack(rastlist, quick = TRUE)
names(rastdata) <- c("elev", "nlcd", "roads", "rough", "shores", 
                     "slope", "urblarge", "urbmed", "urbsmall")
ordered_names <- c("elev", "slope", "rough", "roads", "shores", 
                   "nlcd", "urbsmall", "urbmed", "urblarge")
rastdata <- rastdata[[ordered_names]]

# EXTRACT GIS DATA FOR EACH IMAGE



