#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R CODE FOR FLICKR-NFR CLUSTER ANALYSIS
# Harrison B Goldspiel | harrison.goldspiel@maine.edu
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load custom R functions and settings
source("custom_functions_settings.R")

# Image content clustering summary:
## 1. Mine Flickr data for NFR fro 2012-2016 from open-source API (previously done by AS)
## 2. Run Flickr images through Clarifai image classification algorithm (previously done by AS)
## 3. Perform cluster analysis on rural images to identify core themes of human engagement in rural parts of the NFR
## 4. Assign themes to images based on tag-theme composition

# LOAD PACKAGES
library(beepr)
library(lubridate)
library(tidyverse)
library(tidytext)
library(wordcloud)
library(igraph)
library(reshape2)
library(vegan)
library(recluster)
library(pvclust)
library(cluster)
library(ggpubr)


# Data preparation --------------------------------------------------------

# load raw Flickr data
flickr_all <- read.csv("data/flickr_NFR_all.csv") 
# omit corrupted data (six rows with misplaced or missing column values)
flickr_all <- flickr_all[!is.na(flickr_all$OBJECTID_12),]
# omit non-rural images
flickr_rural <- flickr_all[flickr_all$isUrban == 0,]
# load Flickr data in rural public zones
flickr_public <- read.csv("data/flickr_NFR_public_landwater.csv")
flickr_public$STUSPS <- flickr_all$STUSPS[flickr_all$id %in% flickr_public$id]
# extract rural images from private areas
flickr_private <- flickr_rural[flickr_rural$id %notin% flickr_public$id,]

dim(flickr_all)
# [1] 280503     25
dim(flickr_rural)
# [1] 194682     25
dim(flickr_private)
# [1] 118179     25
dim(flickr_public)
# [1] 81621    22

# omit images without any tags
flickr_rural_tagged <- flickr_rural %>%
  filter(tags1 %notin% " " & tags2 %notin% " " & tags3 %notin% " " &
           tags4 %notin% " " & tags5 %notin% " " & tags6 %notin% " " &
           tags7 %notin% " " & tags8 %notin% " " & tags9 %notin% " " &
           tags10 %notin% " ")

# get random sample of 500 images with both tags and user-provided captions for 
# manual human cross-validation of automated tags versus intended photo target
set.seed(131)
flickr_rural_tagged$caption[flickr_rural_tagged$caption == " "] <- NA
tag_cap_samp <- 
  sample_n(flickr_rural_tagged[!is.na(flickr_rural_tagged$caption),], 500)
write.csv(tag_cap_samp, "data/tag_cap_sample.csv", na = "", row.names = FALSE)

# create new datasets for clustering from the tidy datasets
# omit duplicate images from the Flickr dataset
flickr_rural_tidy <- flickr_rural %>%
  mutate(datetime = mdy_hm(datetaken),
         date = date(datetime),
         year = year(date),
         month = month(date),
         hour = hour(datetime),
         id = as.factor(id)) %>%
  distinct(id, .keep_all = TRUE)

flickr_private_tidy <- flickr_private %>%
  mutate(datetime = mdy_hm(datetaken),
         date = date(datetime),
         year = year(date),
         month = month(date),
         hour = hour(datetime),
         id = as.factor(id)) %>%
  distinct(id, .keep_all = TRUE)

flickr_public_tidy <- flickr_public %>%
  mutate(datetime = ymd_hms(datetaken),
         date = date(datetime),
         year = year(date),
         month = month(date),
         hour = hour(datetime),
         id = as.factor(id)) %>%
  distinct(id, .keep_all = TRUE)

dim(flickr_rural_tidy)
# [1] 148952    29
dim(flickr_private_tidy)
# [1] 93680    29
dim(flickr_public_tidy)
# [1] 58817    26



# Spatiotemporal trends  --------------------------------------------------

# visualize photography over months and time of day
## hourly image trends
hourly.trends.fun <- function(data) {
  states.df.list <- list()
  for(state in c("NY", "VT", "NH", "ME")) {
    states.df.list[[state]] <- data.frame(table(data$hour[data$STUSPS == state]))
    colnames(states.df.list[[state]]) <- c("Hour", "Freq")
  }
  states.hourly <- rbind(states.df.list[["NY"]], states.df.list[["VT"]],
                         states.df.list[["NH"]], states.df.list[["ME"]])
  states.hourly$State <- c(rep("NY", 24), rep("VT", 24), 
                                  rep("NH", 24), rep("ME", 24))
  states.hourly$Hour <- as.integer(states.hourly$Hour)-1
  states.hourly$State = factor(states.hourly$State, 
                               levels = c("NY", "VT", "NH", "ME"))
  return(states.hourly)
}

states.public.hourly <- hourly.trends.fun(flickr_public_tidy)
states.private.hourly <- hourly.trends.fun(flickr_private_tidy)
states.hourly <- hourly.trends.fun(flickr_rural_tidy)

# bar plots separated by state
colors <- c("All rural" = "grey", 
            "Public" = "forestgreen", 
            "Private" = "black")

hourly.trends.plot <-
  ggplot(states.public.hourly, aes(x = Hour, y = Freq)) + 
  geom_col(data = states.hourly, inherit.aes = FALSE,
           aes(x = Hour, y = Freq, fill = "All rural")) +
  geom_col(data = states.hourly, inherit.aes = FALSE,
           aes(x = Hour, y = -Freq, fill = "All rural")) +
  geom_col(aes(fill = "Public"), width = 0.5) + 
  geom_col(data = states.private.hourly, width = 0.5,
           inherit.aes = FALSE, aes(x = Hour, y = -Freq, fill = "Private")) +
  facet_grid(~State) +
  theme_bw() + mythemes + labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = 0:24, labels = c("0", as.character(1:24)), 
                     expand = c(.002,0)) +
  coord_polar(start = -.13, clip = "off") +
  theme(legend.position = "bottom", legend.text = element_text(size = 12),
        panel.border = element_blank(), strip.background = element_blank(), 
        axis.text.x = element_text(size = 11), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), panel.spacing = unit(1, "lines"))

## monthly image trends
monthly.trends.fun <- function(data) {
  states.df.list <- list()
  for(state in c("NY", "VT", "NH", "ME")) {
    states.df.list[[state]] <- data.frame(table(data$month[data$STUSPS == state]))
    colnames(states.df.list[[state]]) <- c("Month", "Freq")
  }
  states.monthly <- rbind(states.df.list[["NY"]], states.df.list[["VT"]],
                         states.df.list[["NH"]], states.df.list[["ME"]])
  states.monthly$State <- c(rep("NY", 12), rep("VT", 12), 
                           rep("NH", 12), rep("ME", 12))
  states.monthly$Month <- as.integer(states.monthly$Month)
  states.monthly$State = factor(states.monthly$State, 
                               levels = c("NY", "VT", "NH", "ME"))
  return(states.monthly)
}
states.public.monthly <- monthly.trends.fun(flickr_public_tidy)
states.private.monthly <- monthly.trends.fun(flickr_private_tidy)
states.monthly <- monthly.trends.fun(flickr_rural_tidy)

monthly.trends.plot <- 
  ggplot(states.public.monthly, aes(x = factor(Month), y = Freq)) + 
  geom_col(data = states.monthly, inherit.aes = FALSE,
           aes(x = factor(Month), y = Freq, fill = "All rural")) +
  geom_col(data = states.monthly, inherit.aes = FALSE,
           aes(x = factor(Month), y = -Freq, fill = "All rural")) +
  geom_col(aes(fill = "Public"), width = 0.5) + 
  geom_col(data = states.private.monthly, width = 0.5,
           inherit.aes = FALSE, aes(x = Month, y = -Freq, fill = "Private")) +
  facet_grid(~State) +
  scale_fill_manual(values = colors) +
  theme_bw() + mythemes + labs(x = NULL, y = NULL, fill = NULL) +
  scale_x_discrete(breaks = c(1:12), labels = month.abb) + 
  coord_polar(start = 6, clip = "off") +
  theme(legend.position = "bottom", legend.text = element_text(size = 12),
        panel.border = element_blank(), panel.spacing = unit(1, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank(), 
        axis.text.x = element_text(size = 11), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())


## annual image trends
yearly.trends.fun <- function(data) {
  states.df.list <- list()
  for(state in c("NY", "VT", "NH", "ME")) {
    states.df.list[[state]] <- data.frame(table(data$year[data$STUSPS == state & 
                                                            data$year <= 2017 & 
                                                            data$year >= 2012]))
    colnames(states.df.list[[state]]) <- c("Year", "Freq")
  }
  states.yearly <- rbind(states.df.list[["NY"]], states.df.list[["VT"]],
                          states.df.list[["NH"]], states.df.list[["ME"]])
  states.yearly$State <- c(rep("NY", 6), rep("VT", 6), 
                            rep("NH", 6), rep("ME", 6))
  states.yearly$State = factor(states.yearly$State, 
                                levels = c("NY", "VT", "NH", "ME"))
  return(states.yearly)
}
states.public.yearly <- yearly.trends.fun(flickr_public_tidy)
states.private.yearly <- yearly.trends.fun(flickr_private_tidy)
states.yearly <- yearly.trends.fun(flickr_rural_tidy)

yearly.trends.plot <- 
  ggplot(states.public.yearly, aes(x = as.numeric(as.character(Year)), y = Freq)) + 
  geom_point(aes(color = "Public"), size = 2) + 
  geom_line(aes(color = "Public")) +
  geom_point(data = states.yearly, inherit.aes = FALSE, size = 2,
             aes(x = as.numeric(as.character(Year)), y = Freq, color = "All rural")) +
  geom_line(data = states.yearly, inherit.aes = FALSE,
             aes(x = as.numeric(as.character(Year)), y = Freq, color = "All rural")) +
  geom_point(data = states.private.yearly, inherit.aes = FALSE, size = 2,
             aes(x = as.numeric(as.character(Year)), y = Freq, color = "Private")) +
  geom_line(data = states.private.yearly, inherit.aes = FALSE,
            aes(x = as.numeric(as.character(Year)), y = Freq, color = "Private")) +
  facet_grid(~State) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = c(2012, 2014, 2016)) +
  scale_y_continuous(breaks = c(2500, 5000, 7500, 10000)) +
  facet_grid(~State) +
  theme_bw() + mythemes + labs(x = NULL, y = "Images (n)", fill = NULL) +
  theme(legend.position = "none", legend.text = element_text(size = 12), 
        strip.background = element_blank(), strip.text.x = element_blank(), 
        axis.text.x = element_text(size = 11))


ggarrange(hourly.trends.plot, monthly.trends.plot, yearly.trends.plot,
          common.legend = TRUE, legend = "bottom",
          nrow = 3, labels = c("A", "B", "C"), align = "v",
          label.y = c(1,1,1.1))

ggsave("figures/rural_image_trends.png", width = 8, height = 7, dpi = 600)


# Organize image tags -----------------------------------------------------

# rank tags by percentile, select tags past certain threshold for clustering
# get list of images with all 10 tags
flickr_rural_tidy_tagged <- flickr_rural_tidy %>%
  filter(tags1 %notin% " " & tags2 %notin% " " & tags3 %notin% " " &
           tags4 %notin% " " & tags5 %notin% " " & tags6 %notin% " " &
           tags7 %notin% " " & tags8 %notin% " " & tags9 %notin% " " &
           tags10 %notin% " ")
flickr_private_tidy_tagged <- flickr_private_tidy %>%
  filter(tags1 %notin% " " & tags2 %notin% " " & tags3 %notin% " " &
           tags4 %notin% " " & tags5 %notin% " " & tags6 %notin% " " &
           tags7 %notin% " " & tags8 %notin% " " & tags9 %notin% " " &
           tags10 %notin% " ")
flickr_public_tidy_tagged <- flickr_public_tidy %>%
  filter(tags1 %notin% " " & tags2 %notin% " " & tags3 %notin% " " &
           tags4 %notin% " " & tags5 %notin% " " & tags6 %notin% " " &
           tags7 %notin% " " & tags8 %notin% " " & tags9 %notin% " " &
           tags10 %notin% " ")

# number of rural photos with fewer than 10 tags
nrow(flickr_rural_tidy)-nrow(flickr_rural_tidy_tagged)
# [1] 27653

# number of rural photos without any tags
flickr_rural_tidy %>%
  filter(tags1 %in% " " & tags2 %in% " " & tags3 %in% " " &
           tags4 %in% " " & tags5 %in% " " & tags6 %in% " " &
           tags7 %in% " " & tags8 %in% " " & tags9 %in% " " &
           tags10 %in% " ") %>%
  nrow()

# rank tags by percentile, select tags past certain threshold for clustering
rural_tidy_tags <- 
  list(flickr_rural_tidy_tagged$tags1, flickr_rural_tidy_tagged$tags2, 
       flickr_rural_tidy_tagged$tags3, flickr_rural_tidy_tagged$tags4, 
       flickr_rural_tidy_tagged$tags5, flickr_rural_tidy_tagged$tags6,
       flickr_rural_tidy_tagged$tags7, flickr_rural_tidy_tagged$tags8, 
       flickr_rural_tidy_tagged$tags9, flickr_rural_tidy_tagged$tags10)

private_tidy_tags <- 
  list(flickr_private_tidy_tagged$tags1, flickr_private_tidy_tagged$tags2, 
       flickr_private_tidy_tagged$tags3, flickr_private_tidy_tagged$tags4, 
       flickr_private_tidy_tagged$tags5, flickr_private_tidy_tagged$tags6,
       flickr_private_tidy_tagged$tags7, flickr_private_tidy_tagged$tags8, 
       flickr_private_tidy_tagged$tags9, flickr_private_tidy_tagged$tags10)

public_tidy_tags <- 
  list(flickr_public_tidy_tagged$tags1, flickr_public_tidy_tagged$tags2, 
       flickr_public_tidy_tagged$tags3, flickr_public_tidy_tagged$tags4, 
       flickr_public_tidy_tagged$tags5, flickr_public_tidy_tagged$tags6,
       flickr_public_tidy_tagged$tags7, flickr_public_tidy_tagged$tags8, 
       flickr_public_tidy_tagged$tags9, flickr_public_tidy_tagged$tags10)

rural_tidy_tag_freq <- as.data.frame(table(unlist(rural_tidy_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / length(rural_tidy_tags[[1]]),
         pct = percentile_rank(freq)) 

private_tidy_tag_freq <- as.data.frame(table(unlist(private_tidy_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / length(private_tidy_tags[[1]]),
         pct = percentile_rank(freq)) 

public_tidy_tag_freq <- as.data.frame(table(unlist(public_tidy_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / length(public_tidy_tags[[1]]),
         pct = percentile_rank(freq))

# tag rankings for all photos
rural_all_tags <- 
  list(flickr_rural_tidy$tags1, flickr_rural_tidy$tags2, 
       flickr_rural_tidy$tags3, flickr_rural_tidy$tags4, 
       flickr_rural_tidy$tags5, flickr_rural_tidy$tags6,
       flickr_rural_tidy$tags7, flickr_rural_tidy$tags8, 
       flickr_rural_tidy$tags9, flickr_rural_tidy$tags10)

private_all_tags <- 
  list(flickr_private_tidy$tags1, flickr_private_tidy$tags2, 
       flickr_private_tidy$tags3, flickr_private_tidy$tags4, 
       flickr_private_tidy$tags5, flickr_private_tidy$tags6,
       flickr_private_tidy$tags7, flickr_private_tidy$tags8, 
       flickr_private_tidy$tags9, flickr_private_tidy$tags10)

public_all_tags <- 
  list(flickr_public_tidy$tags1, flickr_public_tidy$tags2, 
       flickr_public_tidy$tags3, flickr_public_tidy$tags4, 
       flickr_public_tidy$tags5, flickr_public_tidy$tags6,
       flickr_public_tidy$tags7, flickr_public_tidy$tags8, 
       flickr_public_tidy$tags9, flickr_public_tidy$tags10)

rural_tag_freq <- as.data.frame(table(unlist(rural_all_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / length(rural_all_tags[[1]]),
         pct = percentile_rank(freq)) 

private_tag_freq <- as.data.frame(table(unlist(private_all_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / length(private_all_tags[[1]]),
         pct = percentile_rank(freq)) 
  
public_tag_freq <- as.data.frame(table(unlist(public_all_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / length(public_all_tags[[1]]),
         pct = percentile_rank(freq))

# export rankings
write.csv(rural_tidy_tag_freq, "data/rural_tidy_tags.csv")
write.csv(private_tidy_tag_freq, "data/private_tidy_tags.csv")
write.csv(public_tidy_tag_freq, "data/public_tidy_tags.csv")
write.csv(rural_tag_freq, "data/rural_all_tags.csv")
write.csv(private_tag_freq, "data/private_all_tags.csv")
write.csv(public_tag_freq, "data/public_all_tags.csv")

# list of uncommon (N < 5 images) tags
rare_tags <- as.character(rural_tag_freq$tag[rural_tag_freq$freq < 5])


# MCMC cluster analysis ---------------------------------------------------

## (1) create raw rural tag matrix
## (2) random sample 80% of rural tags over 1000 iterations, without replacement, 
## (3) and for each iteration, create tag co-occurence matrix and use Walktrap hierchical clustering algorithm, connecting clusters with Ward's distance, saving the modularity and cluster size
## (4) use those MC simulated modularities to pick optimal number of clusters (K);
## (5) run cluster analysis on full dataset of tags, using K clusters

## (1) create tag matrix
rural_network_in <- 
  flickr_rural_tidy %>%
  dplyr::select(id, tags1, tags2, tags3, tags4, tags5, 
                tags6, tags7, tags8, tags9, tags10) %>%
  melt(id.vars = "id", value.name = "tag") %>%
  # filter to omit rare tags (N < 5) and empty tags
  filter(tag %in% rural_tidy_tag_freq$tag[rural_tidy_tag_freq$freq >= 5] & tag != "") %>%
  # filter to omit overly common tags 
  filter(tag %in% rural_tidy_tag_freq$tag[rural_tidy_tag_freq$prop < 0.10]) %>%
  left_join(rural_tidy_tag_freq, by = "tag") %>%
  dplyr::select(id, tag, freq)

# create co-occurrence matrix
rural.tags.dt <- as.data.frame(crossprod(table(rural_network_in[1:2])))

# reduce matrix appropriately for undirected cluster analysis
rural.tags.v <- colnames(rural.tags.dt)
rural.tags.m <- (as.matrix(data.frame(rural.tags.dt)))
dimnames(rural.tags.m) <- list(rural.tags.v, rural.tags.v)
total.occur <- colSums(rural.tags.m)
rural.tags.m[lower.tri(rural.tags.m, diag=T)] <- 0

# function for MCMC walktrap algorithm for identifying optimal steps + clusters
mcmc.wt <- function(tag.matrix) {
  # MCMC sim code
  n.iter <- 100 # n MC iterations
  step.seq <- 3:8 # step range recommended by Pons & Lapaty (2005)
  cluster.seq <- 2:30 # cluster range
  tune_grid <- expand.grid(
    step       = step.seq,
    cluster    = cluster.seq,
    mod        = NA,
    mod_q025   = NA,
    mod_q975   = NA,
    mod_se     = NA)
  
  modularity.l <- list()
  for (step in step.seq) {
    modularity.m <- matrix(NA, nrow = n.iter, ncol = length(cluster.seq))
    subset.ratio <- 0.8 # 80% of the tag dataset
    for (iter.idx in 1:n.iter) {
      n.tags <- nrow(tag.matrix)
      subset.idx <- (sample(1:n.tags, size = floor(n.tags*subset.ratio), 
                            replace = F))
      
      # getting random subset of rural tag matrix
      tags.m.subset.tmp <- (tag.matrix[subset.idx, subset.idx])
      
      # create undirected adjacency matrix from random tag matrix
      fl.graph.subset.tmp <- graph_from_adjacency_matrix(tags.m.subset.tmp,
                                                         weighted=TRUE, 
                                                         mode="undirected",
                                                         diag=TRUE)
      
      # clustering based on the subsampled data
      cw.subset.tmp <- cluster_walktrap(fl.graph.subset.tmp, 
                                        weights = E(fl.graph.subset.tmp)$weight,
                                        membership = T,
                                        steps = step)
      
      # extract modularity score from algorithm run
      modularity.tmp.v <- sapply(cluster.seq, FUN = function(x)  
        modularity(fl.graph.subset.tmp, cut_at(cw.subset.tmp, no = x)))
      modularity.m[iter.idx, ] <- modularity.tmp.v
    }
    
    # extract modularity summary stats
    for (cluster in cluster.seq) {
      tune_grid$mod[tune_grid$step == step & tune_grid$cluster == cluster] <- 
        median(modularity.m[, cluster-1])
      tune_grid$mod_q025[tune_grid$step == step & tune_grid$cluster == cluster] <- 
        quantile(modularity.m[, cluster-1], 0.025)
      tune_grid$mod_q975[tune_grid$step == step & tune_grid$cluster == cluster] <- 
        quantile(modularity.m[, cluster-1], 0.975)
      tune_grid$mod_se[tune_grid$step == step & tune_grid$cluster == cluster] <- 
        sd(modularity.m[, cluster-1]/sqrt(nrow(tune_grid)))
    }
    modularity.l[[step]] <- modularity.m
    cat(paste0(step, " steps!", " (", 
               round((((step-min(step.seq))+1)/length(step.seq))*100,2), "%)"))
  }
  return(list(tune_grid = tune_grid, 
              modularity.l = modularity.l,
              cluster_seq = cluster.seq, 
              step_seq = step.seq))
}

set.seed(2000)
rural_mcmc_wt <- mcmc.wt(tag.matrix = rural.tags.m)

## (3) summarize MC modularity statistics

# top 10 step and cluster sizes based on modularity
rural_mcmc_wt$tune_grid %>% 
  dplyr::arrange(desc(mod)) %>%
  filter(cluster != 1) %>%
  head(10)

# optimal tuning parameters
opt.mod <- rural_mcmc_wt$tune_grid[
  which.max(rural_mcmc_wt$tune_grid$mod),]

# plot full 2D modularity grid
ggplot(rural_mcmc_wt$tune_grid, 
       aes(x = cluster, y = step, z = mod, fill = mod)) +
  geom_tile() + 
  geom_tile(data = opt.mod, color = "red", size = 2) +
  scale_x_continuous(breaks = seq(2, max(rural_mcmc_wt$cluster_seq), 2), 
                     expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(min(rural_mcmc_wt$step_seq), 
                                  max(rural_mcmc_wt$step_seq), 2), 
                     expand = c(0, 0)) +
  theme_bw() + mythemes +
  labs(x = "clusters (k)", y = "steps (n)", fill = "modularity (Q)")

ggsave("figures/modularity_tuning_grid.png", width = 10, height = 5, dpi = 600)

# modularity scores per cluster sizes for different steps
ggplot(rural_mcmc_wt$tune_grid, aes(x = cluster, y = mod)) +
  geom_pointrange(aes(ymin = mod_q025, ymax = mod_q975), col = "grey25") +
  facet_wrap(~step) + theme_light() + mythemes

ggsave("figures/modularity_tuning_facets.png")

# modularity scores w/ outer quantiles for different cluster sizes, using optimal step size
# suggests that you don't necessarily get a huge benefit from more than 9 clusters
rural_mcmc_wt$tune_grid %>%
  filter(step == opt.mod$step) %>%
  ggplot(aes(x = cluster, y = mod)) +
  geom_hline(aes(yintercept = opt.mod$mod_q025), lty = "dashed") +
  geom_errorbar(aes(ymin = mod_q025, ymax = mod_q975), width = 0) +
  geom_point(shape = 21, size = 2, fill = "white") + 
  theme_bw() + mythemes

ggsave("figures/modularity_quantiles_optimal_steps.png")

# boxplot of modularity for different cluster sizes, using optimal step size

png("figures/modularity_boxplot_optimal_steps.png", 
    width=12, height=8, units='in', res=300)

boxplot(rural_mcmc_wt$modularity.l[[opt.mod$step]],
        type="l", xlab= "clusters (k)", 
        ylab= "modularity (Q)", 
        main = paste("Modularity changing with k (Walktrap, steps = ", 
                     opt.mod$step, ")"))

dev.off()

## (4) run cluster analysis on full dataset using optimal steps and clusters
# create undirected adjacency matrix from full tag matrix
fl.graph <- graph_from_adjacency_matrix(rural.tags.m,
                                        weighted=TRUE, 
                                        mode="undirected",
                                        diag=TRUE)

# clustering from full matrix
set.seed(666)
cw <- cluster_walktrap(fl.graph, 
                       weights = E(fl.graph)$weight,
                       membership = T,
                       steps = opt.mod$step)

# optimal number of clusters
nclust <- opt.mod$cluster
cw.cut <- cut_at(cw, no = nclust)
modularity(fl.graph, membership = cw.cut)
table(cw.cut)

## eigenvector ranking
### Tag importance based on the whole network
ec <- eigen_centrality(fl.graph)$vector
cl <- closeness(fl.graph)
bt <- betweenness(fl.graph)
dg <- centr_degree(fl.graph)

### output table of tags, community identity, and eigenvector centrality (importance)
cw.comm <- data.frame(
  tag = cw$names,
  cluster = c(cw.cut),
  eigen_centrality = ec,
  betweenness = bt, 
  closeness = cl, 
  degree=dg
)

tag.ranks <- matrix(nrow = max(table(cw.cut)), ncol = nclust)
for (clust in 1:nclust) {
  cw.eigen.rank <- cw.comm %>%
    filter(cluster == clust) %>%
    arrange(desc(eigen_centrality))
  if (
    nrow(cw.eigen.rank) == max(table(cw.cut))
    )
    tag.ranks[,clust] <- cw.eigen.rank$tag
  else
    tag.ranks[,clust] <- 
      c(cw.eigen.rank$tag, 
        rep(NA, length = max(table(cw.cut))-nrow(cw.eigen.rank)))
  }

write.csv(tag.ranks, "data/rural_tag_clusters.csv", na = "", row.names = FALSE)


# Cluster (i.e., theme) assignment based on tag content --------------------

## aggregate clusters into core groups
## ** indicates a megacluster

## 1:  arts (7,19)                         [tech, music, photography]
## 2:  sports (12,16,21,23)                [team sports, boxing, shooting]
## 3:  scenery (10)**                      [mixture of natural and cultural landscape features]
## 4:  food/dining (4,5,6,8,9,25)          [food, dessert, dining, alcohol] 
## 5:  aquatic recreation (3,11)           [watercrafts, fishing, and watersports]
## 6:  biota (1,13,14,20,26)               [flora, fauna, and fungi (and some livestock, cat, and zoo animals)]
## 7:  equestrian (17)                     [horses]
## 8:  people (2)**                        [people]
## 9:  dogs (24)                           [dogs]
## 10: transport (18,22)                   [cars, trucks, buses, trains, bicycles, aircraft]
## 11: structures (15)                     [houses and other human structures]

cw.comm <- cw.comm %>%
  mutate(theme = case_when(cluster %in% c(7,19) ~ "arts",
                           cluster %in% c(12,16,21,23) ~ "sports",
                           cluster == 10 ~ "scenery",
                           cluster %in% c(4,5,6,8,9,25) ~ "food/dining",
                           cluster %in% c(3,11) ~ "aquatics",
                           cluster %in% c(1,13,14,20,26) ~ "biota",
                           cluster == 17 ~ "equestrian",
                           cluster == 2 ~ "people",
                           cluster == 24 ~ "dogs",
                           cluster %in% c(18,22) ~ "transport",
                           cluster == 15 ~ "structures"))

## create table for referencing top tags in each cluster
cw.comm.table <- cw.comm %>%
  group_by(theme) %>%
  arrange(desc(eigen_centrality), .by_group = TRUE) %>%
  ungroup()

write.csv(cw.comm.table, "data/rural_tag_clusters_agg.csv", 
          na = "", row.names = FALSE)

cw.comm.stats <- cw.comm %>%
  group_by(theme) %>%
  arrange(desc(eigen_centrality), .by_group = TRUE) %>%
  summarize(top10tags = c(paste0(tag[1:10], sep = ", ", collapse = "")),
            n_tags = n()) %>%
  ungroup() %>%
  arrange(theme)

write.csv(cw.comm.stats, "data/rural_tag_clusters_agg_sumstats.csv", 
          na = "", row.names = FALSE)

## assign themes to photos based on tags, randomly split ties
assign.theme <- function(photos, clusters) {
  message("Assigning themes...")
  out <- matrix(nrow = nrow(photos), ncol = 2, 
                dimnames = list(c(NULL, NULL), c("id", "theme")))
  for (i in 1:nrow(photos)) {
    tag.seq <- c(photos[i, c("tags1", "tags2", "tags3", "tags4", "tags5",  
                             "tags6", "tags7", "tags8","tags9", "tags10")])
    photo.id <- as.character(photos[i, "id"])
    # replace empty tags with NAs and remove from vector
    tag.seq[tag.seq == " "] <- NA
    tag.seq <- na.omit(tag.seq)
    # match clusters with tags
    match.theme <- function(x) {
      if (x %in% clusters$tag)
        tag.theme <- clusters$theme[clusters$tag == x]
      else 
        tag.theme <- "other"
    }
    theme.seq <- unlist(lapply(tag.seq, match.theme))
    # identify the dominant cluster (split ties randomly)
    photo.theme <- Mode(na.omit(theme.seq[theme.seq != "other"]))
    out[i, ] <- c(photo.id, photo.theme)
  } 
  out.df <- as.data.frame(out)
  colnames(out.df) <- c("id", "theme")
  out.themes <- left_join(photos, out.df, by = "id")
  message("done!")
  return(out.themes)
}

rural_photos_clustered <- assign.theme(photos = flickr_rural_tidy, 
                                       clusters = cw.comm)


# Summary statistics for theme composition --------------------------------

# summarize images by theme and region over year
library(janitor)

# total rural images per year
rural.image.trends <- 
  rural_photos_clustered %>%
  mutate(year = year(date)) %>% 
  filter(year >= 2012 & year <= 2017) %>%
  group_by(year, STUSPS) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  dcast(year ~ STUSPS) %>%
  adorn_totals("row") %>%
  adorn_totals("col")

# total rural images in public areas per year
public.image.trends <- 
  rural_photos_clustered %>%
  mutate(year = year(date),
         public = ifelse(id %in% flickr_public$id, "public", "private")) %>% 
  filter(year >= 2012 & year <= 2017 & public == "public") %>%
  group_by(year, STUSPS) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  dcast(year ~ STUSPS) %>%
  adorn_totals("row") %>%
  adorn_totals("col")

# total rural images in private areas per year
private.image.trends <- 
  rural_photos_clustered %>%
  mutate(year = year(date),
         public = ifelse(id %in% flickr_public$id, "public", "private")) %>% 
  filter(year >= 2012 & year <= 2017 & public == "private") %>%
  group_by(year, STUSPS) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  dcast(year ~ STUSPS) %>%
  adorn_totals("row") %>%
  adorn_totals("col")

# combined table of rural and public images
all.image.trends <- data.frame(
  Year = c("2012", "2013", "2014", "2015", "2016", "2017", "Total"),
  NY = paste0(rural.image.trends$NY, " (",
              public.image.trends$NY, ")"),
  VT = paste0(rural.image.trends$VT, " (",
              public.image.trends$VT, ")"),
  NH = paste0(rural.image.trends$NH, " (",
              public.image.trends$NH, ")"),
  ME = paste0(rural.image.trends$ME, " (",
              public.image.trends$ME, ")"),
  Total = paste0(rural.image.trends$Total, " (",
              public.image.trends$Total, ")")
)

# save all tables as CSV files
write.csv(rural.image.trends, "data/rural_images_by_year.csv", 
          row.names = FALSE)
write.csv(private.image.trends, "data/private_images_by_year.csv", 
          row.names = FALSE)
write.csv(public.image.trends, "data/public_images_by_year.csv", 
          row.names = FALSE)
write.csv(all.image.trends, "data/rural_public_images_by_year.csv", 
          row.names = FALSE)

# quick bar plot of cluster composition in public and private rural areas
public.private.themes <-
  rural_photos_clustered %>%
  mutate(public = ifelse(id %in% flickr_public$id, "public", "private")) %>%
  group_by(public, theme) %>%
  summarize(count = n()) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(public) %>%
  mutate(prop = round(count / sum(count), 3)) %>%
  ungroup()

public.private.themes

total.rural.themes <- 
  rural_photos_clustered %>%
  group_by(theme) %>%
  summarize(count = n()) %>%
  na.omit() %>%
  ungroup() %>%
  arrange(desc(count)) %>%
  mutate(prop = round(count / sum(count), 3))

total.rural.themes

total.rural.themes.by.state <- rural_photos_clustered %>%
  group_by(theme, STUSPS) %>%
  summarize(count = n()) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(STUSPS) %>%
  arrange(desc(count)) %>%
  mutate(prop = round(count / sum(count), 3)) %>%
  ungroup()

total.rural.themes.by.state %>%
  mutate(State = factor(STUSPS, levels = c("NY", "VT", "NH", "ME"))) %>%
  ggplot(aes(x = State, prop)) +
  geom_col(aes(fill = State)) +
  labs(x = "State", y = "Proportion of images") +
  scale_fill_manual(values = lisa$GeneDavis[1:4]) +
  facet_wrap(~factor(theme, levels = total.themes$theme)) +
  theme_bw() + mythemes + theme(legend.position = c(0.9, 0.15),
                                legend.text = element_text(size = 12),
                                legend.title = element_text(size = 14),
                                axis.text.x = element_text(size = 12), 
                                axis.text.y = element_text(size = 12),
                                axis.title.x = element_text(size = 14),
                                axis.title.y = element_text(size = 14))

ggsave("figures/rural_image_themes_states.png", width = 6, height = 6, dpi = 600)

write.csv(public.private.themes, "data/public_private_theme_props.csv", 
          row.names = FALSE)
write.csv(total.rural.themes, "data/rural_theme_props.csv", 
          row.names = FALSE)
write.csv(total.rural.themes.by.state, "data/rural_state_theme_props.csv", 
          row.names = FALSE)

# composition of themes in private areas
private.themes.plot <- 
  public.private.themes %>%
  filter(public == "private") %>%
  mutate(theme = factor(theme, levels = rev(total.themes$theme))) %>%
  ggplot(aes(x = theme, y = prop)) +
  scale_y_reverse(limits=c(1,0)) +
  geom_col(data = total.themes, 
           aes(x = factor(theme, levels = rev(theme)),
               y = prop), col = "grey", fill = "grey") +
  geom_col(fill = "black", col = "black", width = 0.3) +
  coord_flip() +
  theme_bw() + mythemes +
  labs(title = "Private", x = NULL, y = "Proportion of images") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

# theme rankings (text)
theme.cats.plot <- 
  ggplot(total.rural.themes, aes(x = theme, y = prop)) +
  geom_text(inherit.aes = FALSE, data = total.themes, 
            aes(x = factor(theme, levels = rev(theme)), 
                y = rep(0.5, nrow(total.themes)), label = theme), size = 5) +
  ylim(c(0,1)) +
  coord_flip() +
  theme_bw() + 
  theme(line = element_blank(), rect = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(title = " ", x = NULL, y = " ")

# composition of themes in public areas
public.themes.plot <- 
  public.priavte.themes %>%
  filter(public == "public") %>%
  mutate(theme = factor(theme, levels = rev(total.themes$theme))) %>%
  ggplot(aes(x = theme, y = prop)) +
  ylim(c(0,1)) +
  geom_col(data = total.themes, 
           aes(x = factor(theme, levels = rev(theme)),
               y = prop), col = "grey", fill = "grey") +
  geom_col(fill = "black", col = "black", width = 0.3) +
  coord_flip() +
  theme_bw() + mythemes +
  labs(title = "Public", x = NULL, y = "Proportion of images") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

# pull private/public/theme categories together in one figure
ggarrange(private.themes.plot, theme.cats.plot, public.themes.plot, 
          nrow = 1, widths = c(1,0.35,1), labels = c("A", "", "B"), align = "hv")

ggsave("figures/rural_image_composition.png", width = 10, height = 5, dpi = 600)

# save local environment to feed into model script
save.image(file = "data/cluster_analysis_output.RData")
