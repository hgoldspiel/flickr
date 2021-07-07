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

# note: the protected lands layer from PAD-US includes closed access areas (e.g., a lot in ME)
# maybe create new categorical layer of not protected, protected-closed, protected-restricted, and protected-open
# maybe also include a trails layer
rastfiles <- list.files(path = "data/GIS/baselayers", patter = '.tif$',
                        all.files = TRUE, full.names = TRUE)
rastlist <- lapply(rastfiles, raster)
names(rastlist) <- c("elev", "nlcd", "public", "roads", "rough", 
                     "shores", "slope", "urblarge", "urbmed", "urbsmall")
ordered_names <- c("elev", "slope", "rough", "roads", "shores", 
                   "urbsmall", "urbmed", "urblarge", "nlcd", "public")
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
                   rastlist$urbmed, rastlist$urblarge, rastlist$public, 
                   nlcd_classified)
names(raststack) <- c("elev", "slope", "rough", "roads", "shores", 
                      "urbsmall", "urbmed", "urblarge", "public", "nlcd")

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
  dplyr::select(names(raststack), state, pa, theme) %>%
  bind_rows(bg_tibble)


###############################################################################
# Random forest model development
###############################################################################

# examine collinearity in variables with VIFs
rural_data$pa <- as.factor(rural_data$pa)
rural_data$public <- as.factor(rural_data$public)
rural_data$nlcd <- as.factor(rural_data$nlcd)
corvif(rural_data[,1:10])

# 1. NFR: all rural

## A. Split into training (70%) and testing (30%) datasets
rural_mod_data <- na.omit(rural_data[rural_data$theme == "scenery" | 
                                       rural_data$theme == "all",])
dt <- sort(sample(nrow(rural_mod_data), round(nrow(rural_mod_data)*0.7)))
rural_train <- as.data.frame(rural_mod_data[dt, c(1:10,12)])
rural_test <- as.data.frame(rural_mod_data[-dt, c(1:10,12)])

## B. Tune hyperparameters

library(ranger) # for tuning
library(randomForest) # for fitting tuned RF model
library(Boruta) # for evaluating parameter importance

hyper_grid <- expand.grid(
  mtry        = seq(3, length(names(raststack)), by = 1),
  sample_size = c(0.55, 0.60, 0.65, 0.70, 0.75, 0.80),
  OOB_RMSE    = 0
)

## tuning grid search wtih ranger package
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

## C. Create rf model with tuned hyperparameters
set.seed(123)
rural_rf <- randomForest(
  formula         = pa ~ ., 
  data            = rural_train, 
  ntree           = 1000,
  mtry            = hyper_grid$mtry[1],
  sampsize        = hyper_grid$sample_size[1]*nrow(rural_train),
  importance      = TRUE
)

## out-of-bag (OOB) performance with the caret package
library(caret)
df_rf <- rural_train %>%
  mutate(predicted = predict(rural_rf))

confusionMatrix(df_rf$predicted, df_rf$pa, positive = "1")

# D. Run boruta algorithm to identify important predictors
boruta.bank <- 
  Boruta(pa ~ ., 
         data = rural_train, 
         ntree = 1000, 
         mtry = hyper_grid$mtry[1], 
         sample.fraction = hyper_grid$sample_size[1],
         maxRuns = 500, doTrace = 2)
boruta.bank

# E. Plot variable importance






# test for stack overflow
library(pdp)
set.seed(123)
N <- 40000
rf.data <- data.frame(
  Y = as.factor(rbinom(N, 1, 0.5)),
  X1 = rnorm(N, 800, 50),
  X2 = rnorm(N, 8, 2),
  X3 = rnorm(N, 1600, 500),
  X4 = rnorm(N, 600, 60),
  X5 = rnorm(N, 22, 10),
  X6 = rnorm(N, 20000, 200),
  X7 = rnorm(N, 60000, 600),
  X8 = rnorm(N, 150000, 1000),
  X9 = as.factor(rbinom(N, 1, 0.2)),
  X10 = as.factor(sample(c(1:8), N, replace = T))
)

rf.mod <- randomForest(
  formula         = Y ~ ., 
  data            = rf.data, 
  ntree           = 1000,
  mtry            = 9,
  sampsize        = 0.75*N,
  importance      = TRUE
)

library(h2o)
library(iml)
# 1. create a data frame with just the features
features <- as.data.frame(rf.data) %>% select(-Y)

# 2. Create a vector with the actual responses
response <- rf.data %>% pull(as.numeric(Y))

# 3. Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results[[3L]])
}

predictor.rf <- Predictor$new(
  model = rf.mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)
# compute feature importance with specified loss metric
imp.rf <- FeatureImp$new(predictor.rf, loss = "mse")

# plot output
p2 <- plot(imp.rf)

# pdp
rf.X1 <- Partial$new(predictor.rf, "X1", ice = TRUE, grid.size = 50)
rf.X1$center(min(features$X1))
p2 <- plot(rf.X1) + ggtitle("RF")

rf.mod %>%  # the %>% operator is read as "and then"
  partial(pred.var = "X1", rug = TRUE, grid.resolution = 30) %>%
  autoplot(smooth = TRUE, ylab = expression(f(X1))) +
  theme_light()

rf.mod %>%  # the %>% operator is read as "and then"
  partial(pred.var = "X10", prob = TRUE) %>%
  autoplot() +
  theme_light()

library(plotmo)
plotmo(rf.mod, pmethod="partdep", degree2 = FALSE) # plot partial dependencies

# Get variable importance measures
imp_df <- data.frame(importance(rf.mod, scale = FALSE, type = 1))
# Tidy up and sort the data frame
imp_df <- imp_df %>% 
  mutate(names = rownames(imp_df)) %>% 
  arrange(desc(MeanDecreaseAccuracy))

partial(rf.mod, pred.var = "X1", plot = TRUE, rug = TRUE, plot.engine = "ggplot2")







## create table of importance statistics
boruta_df <- attStats(boruta.bank)
boruta_df$Variable <- rownames(boruta_df)
boruta_df$Scale <- "Northern Forest"
boruta_df$Rel_imp <- boruta_df$meanImp/max(boruta_df$maxImp)

## create list to store boruta results for each outcome
boruta_results <- list()
boruta_results[["Northern Forest"]] <- boruta_df

## plot boruta variable importance
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
  stat_summary() + coord_flip() + 
  scale_fill_manual(values = c("#009E73", "gold", "#D55E00")) +
  labs(y = "Relative importance", x = "Variable") +
  theme_bw()

## E. Plot partial response curves of predictors


rural_rf %>%  # the %>% operator is read as "and then"
  partial(pred.var = "elev", rug = TRUE, grid.resolution = 30) %>%
  autoplot(smooth = TRUE) +
  theme_light() + labs(x = "Elevation (m)", y = "Probability of visitation")


library(plotmo)
plotmo(rural_rf, pmethod="partdep", degree2 = FALSE, nrug = TRUE) # plot partial dependencies



library(pdp)
library(edarf)
source("perfectPartialPlot.R")
pd_df <- partial_dependence(fit = rural_rf,
                            vars = colnames(rural_train)[1:10],
                            data = rural_train,
                            n = c(100, 200))



perfectPartialPlot(df = pd_df, x = colnames(rural_train)[1:10], y = "1")


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