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
library(ggpubr)
library(dplyr)
library(caret) # for evaluating model performance
library(pdp) # for making partial dependence plots
library(reshape2) 
library(ggplot2)
library(lisa)

# LOAD GIS DATA
NFR <- st_read("data/GIS/baselayers/NFR.shp") # NFR extent
rural_zone <- st_read("data/GIS/rural_regions.shp") # rural area of NFR
public_zone <- st_read("data/GIS/baselayers/PAD_US_NFR_open_WGS.shp") # public
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
bg_tibble <- as_tibble(bg_pts_geo) %>%
  mutate(state = bg_pts_state,
         pa = 0,
         theme = "all")

# CREATE RANDOM BACKGROUND POINTS IN PUBLIC ZONE (and link spatial covariates)
public_zone_sf <- st_as_sf(public_zone)
bg_pts_pub <-st_as_sf(st_sample(public_zone, 10000))

bg_pts_pub_geo <- matrix(nrow = 10000, 
                         ncol = length(names(raststack)),
                         dimnames = list(rownames(bg_pts_pub), names(raststack)))
for(i in names(raststack)) {  
  bg_pts_pub_geo[,i] <- raster::extract(raststack[[i]], bg_pts_pub)
}

bg_pts_pub_state <- st_join(bg_pts_pub, states)$STUSPS
bg_pub_tibble <- as_tibble(bg_pts_pub_geo) %>%
  select(-public) %>%
  mutate(state = bg_pts_pub_state,
         theme = "all")

# COMBINE PRESENCE AND BACKGROUND POINTS IN ONE DATA FRAME
rural_data <- as_tibble(rural_images_geo) %>% 
  dplyr::select(names(raststack), state, pa, theme) %>%
  bind_rows(bg_tibble) %>%
  mutate(pa = as.factor(pa),
         public = as.factor(public),
         nlcd = as.factor(nlcd))

public_data <- as_tibble(rural_images_geo) %>% 
  filter(public == 1) %>%
  dplyr::select(names(raststack), state, pa, theme) %>% 
  select(-public) %>%
  bind_rows(bg_pub_tibble) %>%
  mutate(pa = as.factor(if_else(is.na(pa), 0, pa)),
         nlcd = as.factor(nlcd))

###############################################################################
# Random forest model development
###############################################################################

# examine collinearity in variables with VIFs
corvif(rural_data[,1:10])

# function for random forest models
flickr.mod.fun <- function(data, num.vars) {
  models <- list()
  training <- list()
  testing <- list()
  diagnostics <- tibble(from = c(rep(c("NFR", "NY", "VT", "NH", "ME"),each = 5)),
                        to = c(rep(c("NFR", "NY", "VT", "NH", "ME"),5)),
                        accuracy = rep(NA, 25),
                        accuracy_lwr = rep(NA, 25),
                        accuracy_upr = rep(NA, 25),
                        kappa = rep(NA, 25),
                        sensitivity = rep(NA, 25),
                        specificity = rep(NA, 25))
  ## make empty output table for variable importance metrics
  boruta.results <- list()
  ## make empty list for partial dependence predictions
  pdp.data <- list()
  for(region in c("NFR", "NY", "VT", "NH", "ME")) {
    if(region != "NFR") {
      input_data <- data[data$state == region,] 
    } else {
      input_data <- data
    }
    mod_data <- na.omit(input_data[input_data$theme == "scenery" | 
                                     input_data$theme == "natural life" |
                                     input_data$theme == "all",])
    dt <- sort(sample(nrow(mod_data), round(nrow(mod_data)*0.7)))
    train <- as.data.frame(mod_data[dt, c(1:num.vars, num.vars + 2)])
    test <- as.data.frame(mod_data[-dt, c(1:num.vars, num.vars + 2)])
    
    # save training and testing data to list for export
    training[[region]] <- train
    testing[[region]] <- test
    
    ## B. Tune hyperparameters
    hyper_grid <- expand.grid(
      mtry        = seq(3, num.vars, by = 1),
      sample_size = c(0.55, 0.60, 0.65, 0.70, 0.75, 0.80),
      OOB_RMSE    = 0
    )
    
    ## tuning grid search with ranger package
    for(i in 1:nrow(hyper_grid)) {
      # train model
      model <- ranger(
        formula         = pa ~ .,
        data            = train, 
        num.trees       = 500,
        mtry            = hyper_grid$mtry[i],
        sample.fraction = hyper_grid$sample_size[i],
        seed            = 123
      )
      
      # add OOB error to grid
      hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
    }
    
    hyper_grid <- hyper_grid %>% 
      dplyr::arrange(OOB_RMSE)
    
    ## C. Create rf model with tuned hyperparameters
    set.seed(123)
    
    # build with random forest for getting the confusion matrix
    rf.mod <- randomForest(
      formula         = pa ~ ., 
      data            = train, 
      ntree           = 500,
      mtry            = hyper_grid$mtry[1],
      sampsize        = hyper_grid$sample_size[1]*nrow(train),
      seed            = 123,
      probability     = TRUE
    )
    
    # save model to list for export
    models[[region]] <- rf.mod
    
    # build with ranger for plotting PDPs
    ranger.mod <- ranger(
      formula         = pa ~ ., 
      data            = train, 
      num.trees       = 500,
      mtry            = hyper_grid$mtry[1],
      sample.fraction = hyper_grid$sample_size[1],
      seed            = 123,
      probability     = TRUE
    )
    
    ## cross-validation (down-scaling, up-scaling, and between states)
    for(to in c("NFR", "NY", "VT", "NH", "ME")) {
      if(to == region) {
        df_rf <- test %>% mutate(predicted = predict(rf.mod, newdata = test))
        diag.obj <- confusionMatrix(df_rf$predicted, df_rf$pa, positive = "1")
        diagnostics$accuracy[diagnostics$from == region & 
                               diagnostics$to == region] <- 
          diag.obj$overall["Accuracy"]
        diagnostics$accuracy_lwr[diagnostics$from == region & 
                                   diagnostics$to == region] <- 
          diag.obj$overall["AccuracyLower"]
        diagnostics$accuracy_upr[diagnostics$from == region & 
                                   diagnostics$to == region] <- 
          diag.obj$overall["AccuracyUpper"]
        diagnostics$kappa[diagnostics$from == region & 
                            diagnostics$to == region] <-  
          diag.obj$overall["Kappa"]
        diagnostics$sensitivity[diagnostics$from == region & 
                                  diagnostics$to == region] <- 
          diag.obj$byClass["Sensitivity"]
        diagnostics$specificity[diagnostics$from == region & 
                                  diagnostics$to == region] <- 
          diag.obj$byClass["Specificity"]
      } else {
        if(to != region) {
          testdat <- na.omit(
            data[data$state == to & 
                   c(data$theme == "scenery" | 
                       data$theme == "natural life" |
                       data$theme == "all"), c(1:num.vars, num.vars + 2)])
          pred <- predict(rf.mod, newdata = testdat)
          diag.obj <- confusionMatrix(pred$predicted, testdat$pa, positive = "1")
          diagnostics$accuracy[diagnostics$from == region & 
                                 diagnostics$to == to] <- 
            diag.obj$overall["Accuracy"]
          diagnostics$accuracy_lwr[diagnostics$from == region & 
                                     diagnostics$to == to] <- 
            diag.obj$overall["AccuracyLower"]
          diagnostics$accuracy_upr[diagnostics$from == region & 
                                     diagnostics$to == to] <- 
            diag.obj$overall["AccuracyUpper"]
          diagnostics$kappa[diagnostics$from == region & 
                              diagnostics$to == to] <-  
            diag.obj$overall["Kappa"]
          diagnostics$sensitivity[diagnostics$from == region & 
                                    diagnostics$to == to] <- 
            diag.obj$byClass["Sensitivity"]
          diagnostics$specificity[diagnostics$from == region & 
                                    diagnostics$to == to] <- 
            diag.obj$byClass["Specificity"]
        }
      }
    }
    
    # D. Run boruta algorithm to identify important predictors
    boruta.bank <- 
      Boruta(pa ~ ., 
             data = train, 
             ntree = 500, 
             mtry = hyper_grid$mtry[1], 
             sample.fraction = hyper_grid$sample_size[1],
             maxRuns = 500, doTrace = 2)
    
    ## create table of importance statistics
    boruta_df <- attStats(boruta.bank)
    boruta_df$Variable <- rownames(boruta_df)
    boruta_df$Scale <- region
    boruta_df$rel_imp <- boruta_df$meanImp/max(boruta_df$maxImp)
    
    ## create list to store boruta results for each outcome
    boruta.results[[region]] <- boruta_df
    
    ## E. Plot partial response curves of predictors
    # get partial dependence predictions
    pred.vars <- colnames(train)[1:num.vars]
    for(variable in pred.vars) {
      pdp.data[[region]][[variable]] <- 
        ranger.mod %>% 
        partial(pred.var = variable, rug = TRUE, grid.resolution = 30, 
                prob = TRUE, train = train[sample(nrow(train), 500),]) %>%
        as_tibble() %>%
        mutate(region = region)
    }
    print(paste0(region, " done!"))
  }
  # merge together predictions across regions
  pdp.data.merge <- list()
  for(variable in pred.vars) {
    pdp.data.merge[[variable]] <- bind_rows(pdp.data[["NFR"]][[variable]], 
                                            pdp.data[["NY"]][[variable]],
                                            pdp.data[["VT"]][[variable]],
                                            pdp.data[["NH"]][[variable]],
                                            pdp.data[["ME"]][[variable]])
    pdp.data.merge[[variable]]$Region <- 
      factor(pdp.data.merge[[variable]]$region, 
             levels = c("NFR", "NY", "VT", "NH", "ME"))
  }
  
  # combine Boruta results for plotting
  boruta.df <- do.call(rbind, boruta.results)
  boruta.plot <- boruta.df %>%
    mutate(Region = factor(Scale, levels = c("NFR", "NY", "VT", "NH", "ME"))) %>%
    ggplot(aes(x = reorder(Variable, rel_imp, mean), y = rel_imp, shape = Region)) +
    geom_boxplot(inherit.aes = FALSE, 
                 aes(x = reorder(Variable, rel_imp, mean), y = rel_imp), 
                 fill = "grey", color = "grey", alpha = 0.50) +
    geom_point() + coord_flip() + 
    labs(y = "Relative importance", x = "Variable") +
    theme_bw()
  
  # combine PDP results for plotting
  # create list for custom x-axis labels
  if(num.vars == 10) {
    xvars <- c("Elevation (m)", "Slope (%)", "Roughness", "Dist. to road (m)",
               "Dist. to shore (m)", "Small urban dist. (m)", 
               "Medium urban dist. (m)", "Large urban dist. (m)", 
               "Public land", "Land cover")
  } else {
    xvars <- c("Elevation (m)", "Slope (%)", "Roughness", "Dist. to road (m)",
               "Dist. to shore (m)", "Small urban dist. (m)", 
               "Medium urban dist. (m)", "Large urban dist. (m)", "Land cover")
  }
  names(xvars) <- colnames(train)[1:num.vars]
  
  # make PDPs for each x variable
  pdp.var.plots <- list()
  for(variable in pred.vars) {
    if(is.factor(as.data.frame(pdp.data.merge[[variable]])[,variable])) {
      pdp.var.plots[[variable]] <- 
        pdp.data.merge[[variable]] %>%
        ggplot(aes_string(x = variable, 
                          y = "yhat")) + 
        geom_point(aes(color = Region)) + ylim(0,1) +
        theme_light() + 
        scale_color_manual(values = c("black", lisa$GeneDavis[1:4])) +
        theme(legend.position = "none") +
        labs(x = xvars[variable], y = NULL)
    } else {
      pdp.var.plots[[variable]] <- 
        pdp.data.merge[[variable]] %>%
        ggplot(aes_string(x = variable, 
                          y = "yhat")) + 
        geom_line(aes(color = Region)) + 
        geom_smooth(aes(color = region), se = FALSE) + ylim(0,1) +
        theme_light() + 
        scale_color_manual(values = c("black", lisa$GeneDavis[1:4])) +
        theme(legend.position = "none") +
        labs(x = xvars[variable], y = NULL)
    }
  }
  
  if(num.vars == 10) {
    leg <- get_legend(pdp.var.plots$elev + theme(legend.position = "right"))
    pdp.out <- ggarrange(pdp.var.plots$elev, pdp.var.plots$slope, 
                         pdp.var.plots$rough, pdp.var.plots$roads, 
                         pdp.var.plots$shores, pdp.var.plots$urbsmall, 
                         pdp.var.plots$urbmed, pdp.var.plots$urblarge, 
                         pdp.var.plots$public, pdp.var.plots$nlcd, leg, 
                         nrow = 3, ncol = 4)
    pdp.out <- annotate_figure(pdp.out,
                               left = text_grob("Probability of nature recreation", 
                                                rot = 90))
  } else {
    pdp.out <- ggarrange(pdp.var.plots$elev, pdp.var.plots$slope, 
                         pdp.var.plots$rough, pdp.var.plots$roads, 
                         pdp.var.plots$shores, pdp.var.plots$urbsmall, 
                         pdp.var.plots$urbmed, pdp.var.plots$urblarge, 
                         pdp.var.plots$nlcd, 
                         nrow = 3, ncol = 3, 
                         common.legend = TRUE, legend = "bottom")
    pdp.out <- annotate_figure(pdp.out,
                               left = text_grob("Probability of nature recreation", 
                                                rot = 90))
  }
  return(list(data = list(training = training,
                          testing = testing),
              results = list(models = models, 
                             results = diagnostics, 
                             boruta = boruta.df),
              plots = list(boruta = boruta.plot, 
                           pdp = pdp.out)))
  beep(3)
}

rural.model.results <- flickr.mod.fun(data = rural_data, num.vars = 10)
public.model.results <- flickr.mod.fun(data = public_data, num.vars = 9)

  ###############################################################################
  # Model visualizations and predictions
  ###############################################################################
  
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