#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R CODE FOR FLICKR-NFR POINT PATTERN ANALYSES
# Harrison B Goldspiel | harrison.goldspiel@maine.edu
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load image data with assigned themes from cluster analysis
load("data/cluster_analysis_output.RData")

# Predictive models for NFR Flickr data
## 1. Link images to geographic raster data
## 2. Generate random background points
## 3. Create predictive models for focal areas (NFR, states)
## 4. Identify relevant geographic drivers of visitation
## 5. Use models to predict CES suitability surfaces for focal areas

# Data preparation ----------------------------------------------------------

# load packages
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
library(caret) 
library(pdp) 
library(reshape2) 
library(ggplot2)
library(lisa)

# load GIS data
NFR <- st_read("data/GIS/baselayers/NFR.shp") # NFR extent
rural_zone <- st_read("data/GIS/rural_regions.shp") # rural extent
urban_zone <- st_read("data/GIS/urban_regions.shp") # urban area
public_zone <- st_read("data/GIS/PADUS_NWI_merge_validate.shp") # public + water
rastfiles <- list.files(path = "data/GIS/baselayers", patter = '.tif$',
                        all.files = TRUE, full.names = TRUE)
rastlist <- lapply(rastfiles, raster)
names(rastlist) <- c("elev", "nlcd", "roads", "rough", "shores", 
                     "slope", "urblarge", "urbmed", "urbsmall")
ordered_names <- c("elev", "slope", "rough", "roads", "shores", 
                   "urbsmall", "urbmed", "urblarge", "nlcd")
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
nlcd_classified <- reclassify(rastlist$nlcd, reclass_m, right = NA, 
                              include.lowest = TRUE)
raststack <- stack(rastlist$elev, rastlist$slope, rastlist$rough,
                   rastlist$roads, rastlist$shores, rastlist$urbsmall,
                   rastlist$urbmed, rastlist$urblarge, nlcd_classified)
names(raststack) <- c("elev", "slope", "rough", "roads", "shores", 
                      "urbsmall", "urbmed", "urblarge", "nlcd")

# extract GIS data for each image (summer only)
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

# make empty sq km grid
NFR.albers <- NFR %>% st_transform(5070)
NFR.grid = st_make_grid(NFR.albers, c(1000, 1000), 
                        what = "polygons", square = TRUE) 
NFR.grid.centroids = 


%>% 
  st_transform(4326) %>%
  st_sf()

NFR.grid.centroids <- st_cast(plot_locations_df, "POLYGON")
plot_locations_df


NFR.grid.clip <- st_intersection(NFR.grid, NFR)



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
  mutate(state = as.factor(bg_pts_pub_state),
         nlcd = as.factor(nlcd),
         urbsmall = urbsmall/1000,
         urbmed = urbmed/1000,
         urblarge = urblarge/1000,
         pa = 0,
         pa = as.factor(pa),
         theme = "all")

# COMBINE PRESENCE AND BACKGROUND POINTS IN ONE DATA FRAME
rural_data <- as_tibble(rural_images_geo) %>% 
  dplyr::select(names(raststack), state, pa, theme, id) %>%
  bind_rows(bg_tibble) %>%
  mutate(pa = as.factor(pa),
         nlcd = as.factor(nlcd),
         urbsmall = urbsmall/1000,
         urbmed = urbmed/1000,
         urblarge = urblarge/1000)

public_data <- rural_data %>%
  filter(id %in% flickr_public$id) %>%
  dplyr::select(-id) %>%
  bind_rows(bg_pub_tibble) %>%
  mutate(rough = (rough - cellStats(rastlist$rough, "mean"))/cellStats(rastlist$rough, "sd"))

rural_data <- dplyr::select(rural_data, -id) %>%
  mutate(rough = (rough - cellStats(rastlist$rough, "mean"))/cellStats(rastlist$rough, "sd"))


# Random forest models ----------------------------------------------------

# examine collinearity in variables with VIFs
corvif(rural_data[,1:9])

# CES themes to subset image data with
CES.themes <- c("scenery", "biota", "aquatics")

# function for random forest models
flickr.mod.fun <- function(data) {
  diagnostics <- tibble(from = c(rep(c("NFR", "NY", "VT", "NH", "ME"),each = 5)),
                        to = c(rep(c("NFR", "NY", "VT", "NH", "ME"),5)),
                        accuracy = rep(NA, 25),
                        accuracy_lwr = rep(NA, 25),
                        accuracy_upr = rep(NA, 25),
                        kappa = rep(NA, 25),
                        sensitivity = rep(NA, 25),
                        specificity = rep(NA, 25))
  ## make empty list for model objects
  models <- list()
  ## make empty output table for variable importance metrics
  boruta.results <- list()
  ## make empty list for partial dependence predictions
  pdp.data <- list()
  ## make empty list for ice plots
  pdp.ice.plots <- list()
  data.IDs <- data %>% mutate(uniqueID = 1:nrow(data))
  for(location in c("NFR", "NY", "VT", "NH", "ME")) {
    if(location != "NFR") {
      input_data <- data.IDs[data.IDs$state == location,] 
    } else {
      input_data <- data.IDs
    }
    mod_data <- na.omit(input_data[input_data$theme %in% CES.themes |
                                     input_data$theme == "all",])
    dt <- sort(sample(nrow(mod_data), round(nrow(mod_data)*0.7)))
    training.IDs <- c(mod_data$uniqueID[dt])
    train <- as.data.frame(mod_data[dt, c(1:9, 11)])
    test <- as.data.frame(mod_data[-dt, c(1:9, 11)])
    
    ## B. Tune hyperparameters
    hyper_grid <- expand.grid(
      mtry        = seq(3, 9, by = 1),
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
        seed            = 123,
        verbose         = FALSE
      )
      # add OOB error to grid
      hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
    }
    
    hyper_grid <- hyper_grid %>% 
      dplyr::arrange(OOB_RMSE)
    
    ## C. Create rf model with tuned hyperparameters
    set.seed(123)
    
    # build model with ranger
    rf.mod <- ranger(
      formula         = pa ~ ., 
      data            = train, 
      num.trees       = 500,
      mtry            = hyper_grid$mtry[1],
      sample.fraction = hyper_grid$sample_size[1],
      seed            = 123,
      verbose         = FALSE
    )
    
    ## cross-validation (within-region, down-scaling, up-scaling, & b/w states)
    for(to in c("NFR", "NY", "VT", "NH", "ME")) {
      # within-region
      if(to == location) {
        pred <- predict(rf.mod, data = test)
        diag.obj <- confusionMatrix(pred$predictions, test$pa, positive = "1")
        diagnostics$accuracy[diagnostics$from == location & 
                               diagnostics$to == location] <- 
          diag.obj$overall["Accuracy"]
        diagnostics$accuracy_lwr[diagnostics$from == location & 
                                   diagnostics$to == location] <- 
          diag.obj$overall["AccuracyLower"]
        diagnostics$accuracy_upr[diagnostics$from == location & 
                                   diagnostics$to == location] <- 
          diag.obj$overall["AccuracyUpper"]
        diagnostics$kappa[diagnostics$from == location & 
                            diagnostics$to == location] <-  
          diag.obj$overall["Kappa"]
        diagnostics$sensitivity[diagnostics$from == location & 
                                  diagnostics$to == location] <- 
          diag.obj$byClass["Sensitivity"]
        diagnostics$specificity[diagnostics$from == location & 
                                  diagnostics$to == location] <- 
          diag.obj$byClass["Specificity"]
      } else if(location == "NFR" & to != "NFR") {
        # scaling model down
        testdat <- na.omit(
          data.IDs[data.IDs$state == to & 
                     (data.IDs$theme %in% CES.themes |
                        data.IDs$theme == "all") &
                     data.IDs$uniqueID != training.IDs, 
                   c(1:9, 11)])
        pred <- predict(rf.mod, data = testdat)
        diag.obj <- confusionMatrix(pred$predictions, testdat$pa, positive = "1")
        diagnostics$accuracy[diagnostics$from == location & 
                               diagnostics$to == to] <- 
          diag.obj$overall["Accuracy"]
        diagnostics$accuracy_lwr[diagnostics$from == location & 
                                   diagnostics$to == to] <- 
          diag.obj$overall["AccuracyLower"]
        diagnostics$accuracy_upr[diagnostics$from == location & 
                                   diagnostics$to == to] <- 
          diag.obj$overall["AccuracyUpper"]
        diagnostics$kappa[diagnostics$from == location & 
                            diagnostics$to == to] <-  
          diag.obj$overall["Kappa"]
        diagnostics$sensitivity[diagnostics$from == location & 
                                  diagnostics$to == to] <- 
          diag.obj$byClass["Sensitivity"]
        diagnostics$specificity[diagnostics$from == location & 
                                  diagnostics$to == to] <- 
          diag.obj$byClass["Specificity"]
      } else if(location != "NFR" & to == "NFR") {
        # scaling model up
        testdat <- na.omit(
          data.IDs[(data.IDs$theme %in% CES.themes |
                      data.IDs$theme == "all") &
                     data.IDs$uniqueID != training.IDs, 
                   c(1:9, 11)])
        pred <- predict(rf.mod, data = testdat)
        diag.obj <- confusionMatrix(pred$predictions, testdat$pa, positive = "1")
        diagnostics$accuracy[diagnostics$from == location & 
                               diagnostics$to == to] <- 
          diag.obj$overall["Accuracy"]
        diagnostics$accuracy_lwr[diagnostics$from == location & 
                                   diagnostics$to == to] <- 
          diag.obj$overall["AccuracyLower"]
        diagnostics$accuracy_upr[diagnostics$from == location & 
                                   diagnostics$to == to] <- 
          diag.obj$overall["AccuracyUpper"]
        diagnostics$kappa[diagnostics$from == location & 
                            diagnostics$to == to] <-  
          diag.obj$overall["Kappa"]
        diagnostics$sensitivity[diagnostics$from == location & 
                                  diagnostics$to == to] <- 
          diag.obj$byClass["Sensitivity"]
        diagnostics$specificity[diagnostics$from == location & 
                                  diagnostics$to == to] <- 
          diag.obj$byClass["Specificity"]
      } else {
        # between states
        testdat <- na.omit(
          data[data$state == to & 
                 (data$theme %in% CES.themes |
                    data$theme == "all"), c(1:9, 11)])
        pred <- predict(rf.mod, data = testdat)
        diag.obj <- confusionMatrix(pred$predictions, testdat$pa, positive = "1")
        diagnostics$accuracy[diagnostics$from == location & 
                               diagnostics$to == to] <- 
          diag.obj$overall["Accuracy"]
        diagnostics$accuracy_lwr[diagnostics$from == location & 
                                   diagnostics$to == to] <- 
          diag.obj$overall["AccuracyLower"]
        diagnostics$accuracy_upr[diagnostics$from == location & 
                                   diagnostics$to == to] <- 
          diag.obj$overall["AccuracyUpper"]
        diagnostics$kappa[diagnostics$from == location & 
                            diagnostics$to == to] <-  
          diag.obj$overall["Kappa"]
        diagnostics$sensitivity[diagnostics$from == location & 
                                  diagnostics$to == to] <- 
          diag.obj$byClass["Sensitivity"]
        diagnostics$specificity[diagnostics$from == location & 
                                  diagnostics$to == to] <- 
          diag.obj$byClass["Specificity"]
      }
    }
    # D. Run boruta algorithm to identify important predictors
    boruta.bank <- 
      Boruta(pa ~ ., 
             data = train, 
             ntree = 500, 
             mtry = hyper_grid$mtry[1], 
             sample.fraction = hyper_grid$sample_size[1],
             maxRuns = 500, doTrace = 0)
    
    ## create table of importance statistics
    boruta_df <- attStats(boruta.bank)
    boruta_df$Variable <- rownames(boruta_df)
    boruta_df$Scale <- location
    boruta_df$rel_imp <- boruta_df$meanImp/max(boruta_df$maxImp)
    
    ## create list to store boruta results for each outcome
    boruta.results[[location]] <- boruta_df
    
    ## E. Plot partial response curves of predictors
    # make rf model with "probability = TRUE" using ranger to get PDPs
    rf.mod.pred <- ranger(
      formula         = pa ~ ., 
      data            = train, 
      num.trees       = 500,
      mtry            = hyper_grid$mtry[1],
      sample.fraction = hyper_grid$sample_size[1],
      seed            = 123,
      probability     = TRUE,
      verbose         = FALSE
    )
    
    # save model for later prediction maps
    models[[location]] <- rf.mod.pred
    
    # get PDP curves (which.class = 2 gets the probability of visitation)
    pred.vars <- colnames(train)[1:9]
    for(variable in pred.vars) {
      # ICE + PDP plots
      pdp.ice.plots[[location]][[variable]] <-
        rf.mod.pred %>%
        pdp::partial(pred.var = variable, rug = TRUE, grid.resolution = 100, 
                     prob = TRUE, which.class = 2, ice = TRUE, alpha = 0.1,
                     plot = TRUE, train = train[sample(nrow(train), 500),],
                     plot.engine = "ggplot2")
      # PDP prediction data frame
      pdp.data[[location]][[variable]] <- 
        rf.mod.pred %>% 
        pdp::partial(pred.var = variable, rug = TRUE, grid.resolution = 100, 
                prob = TRUE, which.class = 2, 
                train = train[sample(nrow(train), 500),]) %>%
        as_tibble() %>%
        mutate(Region = location)
    }
    print(paste0(location, " done!"))
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
      factor(pdp.data.merge[[variable]]$location, 
             levels = c("NFR", "NY", "VT", "NH", "ME"))
  }
  
  # combine Boruta results for plotting
  boruta.df <- do.call(rbind, boruta.results)
  boruta.plot <- boruta.df %>%
    mutate(Region = factor(Scale, levels = c("NFR", "NY", "VT", "NH", "ME"))) %>%
    ggplot(aes(x = reorder(Variable, rel_imp, mean), y = rel_imp, 
               shape = Region, color = Region, fill = Region)) +
    geom_crossbar(inherit.aes = FALSE, 
                  aes(x = reorder(Variable, rel_imp, mean), y = rel_imp),
                  stat = "summary", fun.min = min, fun.max = max, fun = mean,
                  fill = "grey", color = "grey", alpha = 0.50) +
    geom_point(size = 2, alpha = 0.9) + coord_flip() + 
    scale_color_manual(values = c("black", lisa$GeneDavis[1:4])) +
    scale_fill_manual(values = c("black", lisa$GeneDavis[1:4])) +
    scale_shape_manual(values = c(21:25)) +
    labs(y = "Relative importance", x = "Variable") +
    theme_bw()
  
  # combine PDP results for plotting
  # create list for custom x-axis labels
  xvars <- c("Elevation (m)", "Slope (%)", "Roughness", "Dist. to road (m)",
             "Dist. to shore (m)", "Small urban dist. (km)", 
             "Medium urban dist. (km)", "Large urban dist. (km)", "Land cover")
  names(xvars) <- colnames(train)[1:9]
  
  # make PDPs for each x variable
  pdp.var.plots <- list()
  for(variable in pred.vars) {
    # include ifelse to make custom axis limits for rural and public data
    if(nrow(data) > 100000) {
      if(is.factor(as.data.frame(pdp.data.merge[[variable]])[,variable])) {
        pdp.var.plots[[variable]] <- 
          pdp.data.merge[[variable]] %>%
          ggplot(aes_string(x = variable, 
                            y = "yhat")) + 
          geom_point(aes(shape = Region, fill = Region), size = 2) +
          scale_color_manual(values = c(rep("black", 5))) +
          scale_fill_manual(values = col2[c(1:5)]) +
          scale_shape_manual(values = c(21:25)) +
          theme_light() + ylim(0.4,1) +
          labs(x = xvars[variable], y = NULL) +
          theme(legend.position = "none", 
                plot.margin=unit(c(0.1,0.4,0.1,0.1),"cm"))
      } else {
        pdp.var.plots[[variable]] <- 
          pdp.data.merge[[variable]] %>%
          ggplot(aes_string(x = variable, 
                            y = "yhat")) + 
          geom_line(aes(color = Region)) + 
          geom_smooth(aes(color = Region), se = FALSE) 
          theme_light() + ylim(0.4,1) +
        scale_color_manual(values = col2[c(1:5)]) +
        labs(x = xvars[variable], y = NULL) +
        theme(legend.position = "none", 
              plot.margin=unit(c(0.1,0.4,0.1,0.1),"cm"))
      }
    } else {
      if(is.factor(as.data.frame(pdp.data.merge[[variable]])[,variable])) {
        pdp.var.plots[[variable]] <- 
          pdp.data.merge[[variable]] %>%
          ggplot(aes_string(x = variable, 
                            y = "yhat")) + 
          geom_point(aes(shape = Region, fill = Region), size = 2) +
          scale_color_manual(values = c(rep("black", 5))) +
          scale_fill_manual(values = col2[c(1:5)]) +
          scale_shape_manual(values = c(21:25)) +
          theme_light() + ylim(0.25,1) +
          labs(x = xvars[variable], y = NULL) +
          theme(legend.position = "none", 
                plot.margin=unit(c(0.1,0.4,0.1,0.1),"cm"))
      } else {
        pdp.var.plots[[variable]] <- 
          pdp.data.merge[[variable]] %>%
          ggplot(aes_string(x = variable, 
                            y = "yhat")) + 
          geom_line(aes(color = Region)) + 
          geom_smooth(aes(color = Region), se = FALSE) + 
          theme_light() + ylim(0.25,1) +
          scale_color_manual(values = col2[1:5]) +
          labs(x = xvars[variable], y = NULL) +
          theme(legend.position = "none", 
                plot.margin=unit(c(0.1,0.4,0.1,0.1),"cm"))
      }
    }
  }
  pdp.out <- 
    ggarrange(pdp.var.plots$elev, 
              pdp.var.plots$slope + theme(axis.ticks.y = element_blank(), 
                                          axis.text.y = element_blank()),
              pdp.var.plots$rough + theme(axis.ticks.y = element_blank(), 
                                          axis.text.y = element_blank()),
              pdp.var.plots$urbsmall, 
              pdp.var.plots$urbmed + theme(axis.ticks.y = element_blank(), 
                                           axis.text.y = element_blank()),
              pdp.var.plots$urblarge + theme(axis.ticks.y = element_blank(), 
                                             axis.text.y = element_blank()),
              pdp.var.plots$roads, 
              pdp.var.plots$shores + theme(axis.ticks.y = element_blank(), 
                                           axis.text.y = element_blank()),
              pdp.var.plots$nlcd + theme(axis.ticks.y = element_blank(), 
                                         axis.text.y = element_blank()),
              nrow = 3, ncol = 3, align = "h",
              labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
              label.x = rep(c(0.12, 0.05, 0.05),3), label.y = 0.95, 
              font.label = list(size = 12), common.legend = TRUE, legend = "bottom")
  pdp.out <- annotate_figure(pdp.out,
                             left = text_grob(
                               "Probability of nature-based engagement", 
                               rot = 90))
  return(list(models = models,
              results = list(results = diagnostics, 
                             boruta = boruta.df,
                             pdp = pdp.data),
              plots = list(boruta = boruta.plot, 
                           pdp.all = pdp.var.plots,
                           pdp.ice = pdp.ice.plots,
                           pdp = pdp.out)))
}

rural.model.results <- flickr.mod.fun(data = rural_data)
public.model.results <- flickr.mod.fun(data = public_data)


# Model visualizations and predictions ------------------------------------

# print model validation tables
## rural models
rural.accuracy <- rural.model.results$results$results %>%
  select(from, to, accuracy) %>%
  mutate(from = factor(from, levels = c("NFR", "NY", "VT", "NH", "ME")),
         to = factor(to, levels = c("NFR", "NY", "VT", "NH", "ME"))) %>%
  dcast(from ~ to) %>%
  mutate(statistic = "accuracy")

rural.sensitivity <- rural.model.results$results$results %>%
  select(from, to, sensitivity) %>%
  mutate(from = factor(from, levels = c("NFR", "NY", "VT", "NH", "ME")),
         to = factor(to, levels = c("NFR", "NY", "VT", "NH", "ME"))) %>%
  dcast(from ~ to) %>%
  mutate(statistic = "sensitivity")

rural.specificity <- rural.model.results$results$results %>%
  select(from, to, specificity) %>%
  mutate(from = factor(from, levels = c("NFR", "NY", "VT", "NH", "ME")),
         to = factor(to, levels = c("NFR", "NY", "VT", "NH", "ME"))) %>%
  dcast(from ~ to) %>%
  mutate(statistic = "specificity")

rural.diagnostics <- rbind(rural.accuracy, rural.sensitivity, rural.specificity) 
write.csv(rural.diagnostics, "data/rural_model_diagnostics.csv", row.names = FALSE)

## public models
public.accuracy <- public.model.results$results$results %>%
  select(from, to, accuracy) %>%
  mutate(from = factor(from, levels = c("NFR", "NY", "VT", "NH", "ME")),
         to = factor(to, levels = c("NFR", "NY", "VT", "NH", "ME"))) %>%
  dcast(from ~ to) %>%
  mutate(statistic = "accuracy")

public.sensitivity <- public.model.results$results$results %>%
  select(from, to, sensitivity) %>%
  mutate(from = factor(from, levels = c("NFR", "NY", "VT", "NH", "ME")),
         to = factor(to, levels = c("NFR", "NY", "VT", "NH", "ME"))) %>%
  dcast(from ~ to) %>%
  mutate(statistic = "sensitivity")

public.specificity <- public.model.results$results$results %>%
  select(from, to, specificity) %>%
  mutate(from = factor(from, levels = c("NFR", "NY", "VT", "NH", "ME")),
         to = factor(to, levels = c("NFR", "NY", "VT", "NH", "ME"))) %>%
  dcast(from ~ to) %>%
  mutate(statistic = "specificity")

public.diagnostics <- rbind(public.accuracy, public.sensitivity, public.specificity)
write.csv(public.diagnostics, "data/public_model_diagnostics.csv", row.names = FALSE)

all.diagnostics <- data.frame(
  test = rep(c("accuracy", "sensitivity", "specificity"), each = 5),
  from = rep(c("NFR", "NY", "VT", "NH", "ME"),3),
  NFR = paste0(round(rural.diagnostics$NFR,3)," (", 
               round(public.diagnostics$NFR,3),")"),
  NY = paste0(round(rural.diagnostics$NY,3)," (", 
              round(public.diagnostics$NY,3),")"),
  VT = paste0(round(rural.diagnostics$VT,3)," (", 
              round(public.diagnostics$VT,3),")"),
  NH = paste0(round(rural.diagnostics$NH,3)," (", 
              round(public.diagnostics$NH,3),")"),
  ME = paste0(round(rural.diagnostics$ME,3)," (", 
              round(public.diagnostics$ME,3),")")
)
write.csv(all.diagnostics, "data/all_model_diagnostics.csv", row.names = FALSE)

# plot variable importance for different spatial extents
rural.boruta <- rural.model.results$results$boruta %>%
  mutate(Access = "A \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ All rural")

public.boruta <- public.model.results$results$boruta %>%
  mutate(Access = "B \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ Public access")

importance.results <- bind_rows(rural.boruta, public.boruta) %>%
  mutate(Access = as.factor(Access))

var.order <- rural.boruta %>%
  group_by(Variable) %>%
  summarize(mean.imp = mean(rel_imp)) %>%
  ungroup() %>%
  arrange(mean.imp) %>%
  pull(Variable)

importance.plot <- 
  importance.results %>%
  mutate(Scale = factor(Scale, levels = c("NFR", "NY", "VT", "NH", "ME"))) %>%
  ggplot(aes(x = factor(Variable, levels = var.order),
             y = rel_imp, 
             shape = Scale, color = Scale, fill = Scale)) +
  geom_crossbar(inherit.aes = FALSE, 
                aes(x = factor(Variable, levels = var.order), y = rel_imp),
                stat = "summary", fun.min = min, fun.max = max, fun = mean,
                fill = "grey", color = "grey", alpha = 0.50, width = 0.8) +
  geom_point(size = 2.5) + 
  coord_flip() + 
  scale_color_manual(values = c(rep("black", 5))) +
  scale_fill_manual(values = col2[c(1:5)]) +
  scale_shape_manual(values = c(21:25)) +
  labs(y = "Relative importance", x = "Variable") +
  facet_wrap(~ Access) +
  theme_bw() + mythemes + theme(strip.background = element_blank(),
                                strip.text = element_text(hjust = 0),
                                panel.spacing = unit(2, "lines"))

importance.plot
ggsave("figures/rf_importance.png", width = 10, height = 5, dpi = 600)

# PDPs combining across rural and public
pdp.data.merge1 <- list()
for(variable in c("elev", "slope", "rough", "shores", 
                    "roads", "urbsmall", "urbmed", "urblarge")) {
    pdp.data.merge1[[variable]] <- rural.model.results$results$pdp[["NFR"]][[variable]]
    pdp.data.merge1[[variable]]$Access <- "All rural"
}

pdp.data.merge2 <- list()
for(variable in c("elev", "slope", "rough", "shores", 
                    "roads", "urbsmall", "urbmed", "urblarge")) {
    pdp.data.merge2[[variable]] <- public.model.results$results$pdp[["NFR"]][[variable]]
    pdp.data.merge2[[variable]]$Access <- "Public"
}

xvars <- c("Elevation (m)", "Slope (%)", "Roughness", "Dist. to road (m)",
           "Dist. to shore (m)", "Small urban dist. (km)", 
           "Medium urban dist. (km)", "Large urban dist. (km)")
names(xvars) <- c("elev", "slope", "rough", "shores", 
                  "roads", "urbsmall", "urbmed", "urblarge")
pdp.NFR.plots <- list()
for(variable in c("elev", "slope", "rough", "shores", 
                  "roads", "urbsmall", "urbmed", "urblarge")) {
  pdp.all <- bind_rows(as.data.frame(pdp.data.merge1[variable]), 
                       as.data.frame(pdp.data.merge2[variable])) 
  colnames(pdp.all) <- c(variable, "yhat", "region", "Access")
  pdp.NFR.plots[[variable]] <- 
    ggplot(pdp.all, aes_string(x = variable, y = "yhat", color = "Access")) +
    geom_line() + 
    geom_smooth(se = FALSE) + 
    theme_light() + 
    scale_color_manual(values = c("grey", "forestgreen")) +
    labs(x = xvars[variable], y = NULL) + ylim(0.25, 1) +
    theme(plot.margin=unit(c(0.1,0.4,0.1,0.1),"cm"))
}

pdp.data.merge1 <- rural.model.results$results$pdp[["NFR"]][["nlcd"]]
pdp.data.merge1$Access <- "All rural"

pdp.data.merge2 <- public.model.results$results$pdp[["NFR"]][["nlcd"]]
pdp.data.merge2$Access <- "Public"
pdp.all <- bind_rows(pdp.data.merge1, pdp.data.merge2)

pdp.NFR.plots[["nlcd"]] <- 
  ggplot(pdp.all, aes(x = nlcd, y = yhat, fill = Access)) +
           geom_point(shape = 21, color = "black", size = 2) + 
           theme_light() + 
           scale_fill_manual(values = c("grey", "forestgreen")) +
           labs(x = "Land cover", y = NULL) + ylim(0.25, 1) +
           theme(plot.margin=unit(c(0.1,0.4,0.1,0.1),"cm"))

pdp.NFR.out <- 
  ggarrange(pdp.NFR.plots$elev, 
            pdp.NFR.plots$slope + theme(axis.ticks.y = element_blank(), 
                                        axis.text.y = element_blank()),
            pdp.NFR.plots$rough + theme(axis.ticks.y = element_blank(), 
                                        axis.text.y = element_blank()),
            pdp.NFR.plots$urbsmall, 
            pdp.NFR.plots$urbmed + theme(axis.ticks.y = element_blank(), 
                                         axis.text.y = element_blank()),
            pdp.NFR.plots$urblarge + theme(axis.ticks.y = element_blank(), 
                                           axis.text.y = element_blank()),
            pdp.NFR.plots$roads, 
            pdp.NFR.plots$shores + theme(axis.ticks.y = element_blank(), 
                                         axis.text.y = element_blank()),
            pdp.NFR.plots$nlcd + theme(axis.ticks.y = element_blank(), 
                                       axis.text.y = element_blank()),
            nrow = 3, ncol = 3, align = "h",
            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
            label.x = rep(c(0.12, 0.05, 0.05),3), label.y = 0.95, 
            font.label = list(size = 12), common.legend = TRUE, legend = "bottom")
pdp.NFR.out <- annotate_figure(pdp.NFR.out,
                           left = text_grob(
                             "Probability of nature-based engagement", 
                             rot = 90))

pdp.NFR.out
ggsave("figures/pdp_NFR.png", width = 7, height = 7, dpi = 600)

# plot PDPs for different spatial extents
## all rural (public + private)
rural.model.results$plots$pdp
ggsave("figures/pdp_rural.png", width = 7, height = 7, dpi = 600)

## public access rural areas
public.model.results$plots$pdp 
ggsave("figures/pdp_public.png", width = 7, height = 7, dpi = 600)

# plot suitability surfaces for different spatial extents 
## load grid
NFR_grid <- read.csv("data/GIS/baselayers/NFR_grid_covariates_spatial_scales.csv")
urban_cells <- NFR_grid[NFR_grid$urban != "",]
public_cells <- NFR_grid[!is.na(NFR_grid$public_land),]

## adjust NLCD levels to match the raster data in the models
NFR_grid_full <- mutate(NFR_grid, 
                        nlcd = case_when(nlcd >= 11 & nlcd <= 12 ~ 1,
                                         nlcd >= 21 & nlcd <= 24 ~ 2,
                                         nlcd == 31 ~ 3,
                                         nlcd >= 41 & nlcd <= 43 ~ 4,
                                         nlcd >= 51 & nlcd <= 52 ~ 5,
                                         nlcd >= 71 & nlcd <= 74 ~ 6,
                                         nlcd >= 81 & nlcd <= 82 ~ 7,
                                         nlcd >= 90 & nlcd <= 95 ~ 8),
                        nlcd = as.factor(nlcd),
                        rough = (rough - cellStats(rastlist$rough, "mean"))/
                          cellStats(rastlist$rough, "sd")) %>%
  dplyr::filter(urban == "" & !is.na(elev) & !is.na(nlcd) & !is.na(roads) & 
                  !is.na(rough) & !is.na(shore) & !is.na(slope) & !is.na(urblarge) &
                  !is.na(urbmed) & !is.na(urbsmall))

NFR_pred_grid <- NFR_grid_full %>%
  dplyr::select(elev, nlcd, roads, rough, shores = shore, 
                slope, urblarge, urbmed, urbsmall)

NFR_rural_pred <- predict(rural.model.results$models$NFR, 
                          data = NFR_pred_grid)

NFR_rural_pred_df <- as.data.frame(NFR_rural_pred$predictions)
NFR_rural_pred_df <- cbind(NFR_rural_pred_df, NFR_grid_full) %>%
  mutate(Lat = (EXT_MIN_Y + EXT_MAX_Y)/2,
         Lon = (EXT_MIN_X + EXT_MAX_X)/2)

# project state shapefile and add to plot
states.proj <- st_read("data/GIS/states_NFR_clipped_UTM.shp") 

## plot predictions from rural model
library(scico)
library(ggnewscale)
rur.pred.map <- 
  ggplot(NFR_rural_pred_df, aes(x = Lon, y = Lat, fill = `1`)) +
  geom_tile() +
  scale_fill_scico(palette = "nuuk") +
  labs(x = "", y = "Latitude", fill = "CES suitability") +
  theme_minimal() +
  ggnewscale::new_scale_fill() +
  geom_tile(data = urban_cells, inherit.aes = FALSE, 
            aes(x = (EXT_MIN_X + EXT_MAX_X)/2,
                y = (EXT_MIN_Y + EXT_MAX_Y)/2,
                fill = "Urban")) +
  scale_fill_manual(values = "black", name = NULL) +
  geom_sf(data = states.proj, inherit.aes = FALSE, 
          fill = "transparent", color = "grey20")
         
rur.pred.map

NFR_public_pred <- predict(public.model.results$models$NFR, 
                           data = NFR_pred_grid)

NFR_public_pred_df <- as.data.frame(NFR_public_pred$predictions)
NFR_public_pred_df <- cbind(NFR_public_pred_df, NFR_grid_full) %>%
  mutate(Lat = (EXT_MIN_Y + EXT_MAX_Y)/2,
         Lon = (EXT_MIN_X + EXT_MAX_X)/2)

## plot predictions from public model (just show public lands)
pub.pred.map <- 
  ggplot(NFR_public_pred_df, 
         aes(x = Lon, y = Lat, fill = `1`)) +
  scale_fill_scico(palette = "nuuk") +
  labs(x = "", y = "Latitude", fill = "CES suitability") +
  theme_minimal() +
  ggnewscale::new_scale_fill() +
  geom_tile(data = urban_cells, inherit.aes = FALSE, 
            aes(x = (EXT_MIN_X + EXT_MAX_X)/2,
                y = (EXT_MIN_Y + EXT_MAX_Y)/2,
                fill = "Urban")) +
  scale_fill_manual(values = "black", name = NULL) +
  geom_sf(data = states.proj, inherit.aes = FALSE, 
          fill = "transparent", color = "grey20")

pub.pred.map

ggarrange(rur.pred.map, pub.pred.map, nrow = 2, labels = c("A", "B"), 
          common.legend = T, legend = "right")

ggsave("figures/CES_suitability_maps.png", width = 6.5, height = 8, dpi = 600)

# R version and package info ----------------------------------------------
sessionInfo()

