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


