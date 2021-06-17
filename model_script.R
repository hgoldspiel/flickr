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

###############################################################################
# Data prep
###############################################################################

# LOAD PACKAGES
library(beepr)
library(sp)
library(sf)
library(raster)
library(ranger)
library(randomForest)
library(Boruta)
library(lubridate)

# LOAD GIS DATA
NFR <- st_read("data/GIS/baselayers/NFR.shp") # NFR extent
rural_zone <- st_read("data/GIS/rural_regions.shp") # rural area of NFR

rastfiles <- list.files(path = "data/GIS/baselayers", patter = '.tif$',
                        all.files = TRUE, full.names = TRUE)
rastlist <- lapply(rastfiles, raster)
names(rastlist) <- c("elev", "nlcd", "protected", "roads", "rough", 
                     "shores", "slope", "urblarge", "urbmed", "urbsmall")
ordered_names <- c("elev", "slope", "rough", "roads", "shores", 
                   "urbsmall", "urbmed", "urblarge", "nlcd", "protected")
rastlist <- rastlist[ordered_names]

# reclassify NLCD into broader eight groupings and split into distinct layers
reclass_m <- matrix(c(0, 1, NA,
                      11, 12, 1,
                      21, 24, 2,
                      31, 31, 3,
                      41, 43, 4,
                      51, 52, 5,
                      71, 74, 6,
                      81, 82, 7,
                      90, 95, 8),
                    ncol = 3, byrow = TRUE)
nlcd_classified <- reclassify(rastlist$nlcd, reclass_m, right = NA, include.lowest = TRUE)
raststack <- stack(rastlist$elev, rastlist$slope, rastlist$rough,
                   rastlist$roads, rastlist$shores, rastlist$urbsmall,
                   rastlist$urbmed, rastlist$urblarge, rastlist$protected, 
                   nlcd_classified)
names(raststack) <- c("elev", "slope", "rough", "roads", "shores", 
                      "urbsmall", "urbmed", "urblarge", "protected", "nlcd")

# EXTRACT GIS DATA FOR EACH IMAGE (SUMMER ONLY)
rural_photos_coords <-
  rural_photos_clustered %>%
  dplyr::filter(!is.na(longitude) & !is.na(latitude) &
                  month(date) >= 5 & month(date) <= 9) %>%
  dplyr::select(id, url, id, uniqueID, owner, date, 
                state = STUSPS, longitude, latitude, theme)

# create shapefile for images from XY coordinates
rural_images <- 
  sf::st_as_sf(rural_photos_coords, 
               coords = c("longitude", "latitude")) %>%
  st_set_crs(4326)

# attach extracted raster values to image dataset
rural_images_geo <- matrix(nrow = nrow(rural_images), 
                           ncol = length(names(raststack)),
                           dimnames = list(rownames(rural_images), 
                                           names(raststack)))
for(i in names(raststack)) {  
  
  rural_images_geo[,i] <- raster::extract(raststack[[i]], rural_images)
}  
rural_images_geo <- cbind(rural_images, rural_images_geo)
rural_images_geo$pa <- 1

# CREATE RANDOM BACKGROUND POINTS IN RURAL ZONE (and link spatial covariates)
states <- st_read("data/GIS/baselayers/states.shp") %>% st_set_crs(4326)
bg_pts <-st_as_sf(st_sample(rural_zone, 10000))
bg_pts_geo <- matrix(nrow = 10000, 
                     ncol = length(names(raststack)),
                     dimnames = list(rownames(bg_pts), names(raststack)))
for(i in names(raststack)) {  
  bg_pts_geo[,i] <- raster::extract(raststack[[i]], bg_pts)
}

bg_pts_state <- st_join(bg_pts, states)$STUSPS

library(dplyr)
bg_tibble <- as_tibble(bg_pts_geo) %>%
  mutate(state = bg_pts_state,
         pa = 0,
         theme = "all")

# COMBINE PRESENCE AND BACKGROUND POINTS IN ONE DATA FRAME
rural_data <- as_tibble(rural_images_geo) %>% 
  select(names(raststack), state, pa, theme) %>%
  bind_rows(bg_tibble)


###############################################################################
# Random forest model development
###############################################################################

# examine collinearity in variables with VIFs
rural_data$pa <- as.factor(rural_data$pa)
rural_data$protected <- as.factor(rural_data$protected)
rural_data$nlcd <- as.factor(rural_data$nlcd)
corvif(rural_data[,1:10])

# 1. NFR: all rural

## split into training (70%) and testing (30%) datasets
rural_mod_data <- na.omit(rural_data[rural_data$theme == "scenery" | 
                                       rural_data$theme == "all",])
dt <- sort(sample(nrow(rural_mod_data), round(nrow(rural_mod_data)*0.7)))
rural_train <- rural_mod_data[dt, c(1:10,12)]
rural_test <- rural_mod_data[-dt, c(1:10,12)]

## tune hyperparameters

library(ranger) # for tuning
library(randomForest) # for fitting tuned RF model
library(Boruta) # for evaluating parameter importance

hyper_grid <- expand.grid(
  mtry        = seq(3, length(names(raststack)), by = 1),
  sample_size = c(0.55, 0.60, 0.65, 0.70, 0.75, 0.80),
  OOB_RMSE    = 0
)

# tuning grid search wtih ranger package
for(i in 1:nrow(hyper_grid)) {
  
  # train model
  model <- ranger(
    formula         = pa ~ .,
    data            = rural_train, 
    num.trees       = 1000,
    mtry            = hyper_grid$mtry[i],
    sample.fraction = hyper_grid$sample_size[i],
    seed            = 123
  )
  
  # add OOB error to grid
  hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
}

hyper_grid <- hyper_grid %>% 
  dplyr::arrange(OOB_RMSE)

head(hyper_grid)

set.seed(123)
rural_rf <- randomForest(
  formula         = pa ~ ., 
  data            = rural_train, 
  ntree           = 1000,
  mtry            = hyper_grid$mtry[1],
  sampsize        = hyper_grid$sample_size[1]*nrow(rural_train),
  importance      = TRUE
)

# run boruta algorithm to distinguish between important, tentative, and non-important predictors
boruta.bank <- 
  Boruta(pa ~ ., 
         data = rural_train, 
         ntree = 1000, 
         mtry = hyper_grid$mtry[1], 
         sample.fraction = hyper_grid$sample_size[1],
         maxRuns = 500, doTrace = 2)
boruta.bank

# create table of importance statistics
boruta_df <- attStats(boruta.bank)
boruta_df$Variable <- rownames(boruta_df)
boruta_df$Scale <- "Northern Forest"
boruta_df$Rel_imp <- boruta_df$meanImp/max(boruta_df$maxImp)

# create list to store boruta results for each outcome
boruta_results <- list()
boruta_results[["Northern Forest"]] <- boruta_df

# plot boruta variable importance
library(reshape2)
library(ggplot2)
boruta.melt <- melt(boruta.bank$ImpHistory)
boruta.melt$rel_imp <- boruta.melt$value/max(boruta.melt$value)
colnames(boruta.melt) <- c("Run", "Variable", "Importance", "rel_imp")
boruta.melt %>%
  left_join(boruta_df, by = "Variable") %>%
  mutate(Decision = factor(decision, levels = c("Confirmed", "Rejected"))) %>%
  na.omit() %>% 
  filter(is.infinite(rel_imp) == FALSE) %>%
  ggplot(aes(x = reorder(Variable, rel_imp, mean), y = rel_imp, fill = Decision)) +
  geom_boxplot() + coord_flip() + 
  scale_fill_manual(values = c("#009E73", "gold", "#D55E00")) +
  labs(y = "Relative importance", x = "Variable") +
  theme_bw()



###############################################################################
# Model visualizations and predictions
###############################################################################

# plot partial response curves

# plot suitability for "scenery" engagement across region

# CREATE COARSE SPATIAL DF FOR PREDICTION
rast.agg.list <- list()
for(i in names(raststack)) {  
  if (i %in% c("protected", "nlcd")) {
    rast.agg.list[[i]] <- aggregate(projectRaster(raststack[[i]], 
                                                  crs = "EPSG:5070", 
                                                  method = "ngb"), 
                                    fact = 5)
  } else {
    rast.agg.list[[i]] <- aggregate(projectRaster(raststack[[i]], 
                                                  crs = "EPSG:5070",
                                                  method = "bilinear"), 
                                    fact = 5)
  }
}
raststack.agg <- stack(rast.agg.list)
geo_df <- data.frame(rasterToPoints(raststack.agg[[names(raststack)]]))

sessionInfo()