# R CODE FOR FLICKR-NFR CLUSTER ANALYSIS
# Harrison B Goldspiel | hbgoldspiel@gmail.com

source("custom_functions_settings.R")

# Analytic summary:
## 1. Mine Flickr data for NFR fro 2012-2016 from open-source API (DONE)
## 2. Run Flickr images through Clarifai image classification algorith (DONE)
## 3. Perform cluster analysis on rural images to identify core themes of human engagement in rural parts of the NFR
## 4. Assign themes to images based on tags
## 5. Visualize hotspots of images (PUD) (DONE for overall rural imagery)
## 6. Use random forest & Boruta to identify core drivers of different themes (do this for each theme/state and theme/NFR total)
## 7. Look at partial effects or GAMs of the different themes

# load packages
library(lubridate)
library(maptools)
library(spatstat)
library(ggmap)
library(rgdal)
library(tidyverse)
library(tidytext)
library(wordcloud)
library(igraph)
library(RColorBrewer)
library(reshape2)
library(rgexf)
library(UserNetR)
library(scales)
library(gplots)
library(raster)
library(xtable)
library(lsa)
library(intergraph)
library(sf)

# PREPARE DATA FOR CLUSTER ANALYSIS
# load raw Flickr data
flickr_all <- read.csv("data/flickr_NFR_all.csv") 
flickr_pub <- read.csv("data/flickr_NFR_public_land.csv")
dim(flickr_all)
# [1] 280509     25
dim(flickr_pub)
# [1] 28461    21

library(clarifai)
secret_id(c("flickr", "flickr-all-scopes"))
get_token()
tag_image_urls(flickr_rur_tidy$url[1])

# omit corrupted data (six rows with misplaced or missing column values)
flickr_all <- flickr_all[!is.na(flickr_all$OBJECTID_12),]
# omit non-rural images
flickr_rural <- flickr_all[flickr_all$isUrban == 0,]
# omit images without any tags
flickr_rural %>%
  filter(!(tags1 %in% " " & tags2 %in% " " & tags3 %in% " " &
             tags4 %in% " " & tags5 %in% " " & tags6 %in% " " &
             tags7 %in% " " & tags8 %in% " " & tags9 %in% " " &
             tags10 %in% " ")) -> flickr_rural_tags
flickr_pub %>%
  filter(!(tags1 %in% " " & tags2 %in% " " & tags3 %in% " " &
             tags4 %in% " " & tags5 %in% " " & tags6 %in% " " &
             tags7 %in% " " & tags8 %in% " " & tags9 %in% " " &
             tags10 %in% " ")) -> flickr_pub_tags
# get random sample of 500 images with both tags and user-provided captions for 
# manual human cross-validation of automated tags versus intended photo target
set.seed(131)
flickr_rural_tags$caption[flickr_rural_tags$caption == " "] <- NA
tag_cap_samp <- sample_n(flickr_rural_tags[!is.na(flickr_rural_tags$caption),], 500)
write.csv(tag_cap_samp, "data/tag_cap_sample.csv", na = "", row.names = FALSE)

# omit duplicate images from the Flickr dataset
flickr_rural %>%
  mutate(datetime = mdy_hm(datetaken),
         date = date(datetime),
         month = month(date),
         hour = hour(datetime)) %>%
  distinct(id, .keep_all = TRUE) -> 
  flickr_rur_tidy

flickr_pub %>%
  mutate(datetime = mdy_hm(datetaken),
         date = date(datetime),
         month = month(date),
         hour = hour(datetime)) %>%
  distinct(id, .keep_all = TRUE) -> 
  flickr_pub_tidy

dim(flickr_rur_tidy)
# [1] 149192     29
dim(flickr_pub_tidy)
# [1] 22325    25


# rank tags by percentile, select tags past certain threshold for clustering
rur_tags <- list(flickr_rur_tags$tags1, flickr_rur_tags$tags2, flickr_rur_tags$tags3, 
             flickr_rur_tags$tags4, flickr_rur_tags$tags5, flickr_rur_tags$tags6,
             flickr_rur_tags$tags7, flickr_rur_tags$tags8, flickr_rur_tags$tags9,
             flickr_rur_tags$tags10)

pub_tags <- list(flickr_pub_tags$tags1, flickr_pub_tags$tags2, flickr_pub_tags$tags3, 
                 flickr_pub_tags$tags4, flickr_pub_tags$tags5, flickr_pub_tags$tags6,
                 flickr_pub_tags$tags7, flickr_pub_tags$tags8, flickr_pub_tags$tags9,
                 flickr_pub_tags$tags10)

rur_tag_freq <- as.data.frame(table(unlist(rur_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / n_images,
         pct = p_rank(freq)) 
  
pub_tag_freq <- as.data.frame(table(unlist(pub_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / n_images,
         pct = p_rank(freq))

write.csv(rur_tag_freq, "data/rural_tags.csv")
write.csv(pub_tag_freq, "data/public_tags.csv")

# word cloud of frequently (i.e., N >= 5 images)  occurring tags in rural areas
png("figures/flickr_rural_toptags.png", 
    width=12, height=8, units='in', res=300)
set.seed(140)   
wordcloud(words = rur_tag_freq$tag, freq = rur_tag_freq$freq, 
          min.freq = 5, max.words = Inf, scale=c(8,0.2), 
          colors=brewer.pal(6, "Dark2"), random.order=FALSE, 
          rot.per=0.15)  
par(mar = rep(0, 4))
dev.off()

# word cloud of frequently (i.e., N >= 5 images)  occurring tags in public lands
png("figures/flickr_public_toptags.png", 
    width=12, height=8, units='in', res=300)
set.seed(140)   
wordcloud(words = pub_tag_freq$tag, freq = pub_tag_freq$freq, 
          min.freq = 5, max.words = Inf, scale=c(8,0.2), 
          colors=brewer.pal(6, "Dark2"), random.order=FALSE, 
          rot.per=0.15)  
par(mar = rep(0, 4))
dev.off()

# RUN MONTE CARLO (MC) CLUSTER ANALYSIS-----------------------------------------
## (1) create raw rural tag matrix
## (2) random sample 80% of rural tags over 1000 iterations, without replacement, 
## (3) and for each iteration, create tag co-occurence matrix and use Walktrap hierchical clustering algorithm, connecting clusters with Ward's distance, saving the modularity and cluster size
## (4) use those MC simulated modularities to pick optimal number of clusters (K);
## (5) run cluster analysis on full dataset of tags, using K clusters

## (1) create original tag matrix
flickr_rur_tags %>%
  dplyr::select(id, tags1, tags2, tags3, tags4, tags5, 
                tags6, tags7, tags8, tags9, tags10) %>%
  melt(id.vars = "id", value.name = "tag") %>%
  # filter to omit rare tags (N < 50)
  # filter(tag %in% pub_tag_freq$tag[pub_tag_freq$freq >= 50]) %>%
  left_join(pub_tag_freq, by = "tag") %>%
  dplyr::select(id, tag, freq) -> public_tags

# create co-occurance matrix
tags.dt <- as.data.frame(crossprod(table(public_tags[1:2])))

# reduce matrix appropriately for undirected cluster analysis
tags.v <- colnames(tags.dt)
tags.m <- (as.matrix(data.frame(tags.dt)))
dimnames(tags.m) <- list(tags.v, tags.v)
total.occur <- colSums(tags.m)
tags.m[lower.tri(tags.m, diag=T)] <- 0

## (2-3) run MC cluster algorithm to obtain robust modularity plot, iterating over a range of random step sizes
n.iter <- 100 # n MC iterations
step.seq <- 4:20 # step range
cluster.seq <- 1:20 # cluster range

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
  set.seed(2000)
  for (iter.idx in 1:n.iter) {
    n.tags <- nrow(tags.m)
    subset.idx <- (sample(1:n.tags, size = floor(n.tags*subset.ratio), 
                          replace = F))
    
    # getting random subset of rural tag matrix
    tags.m.subset.tmp <- (tags.m[subset.idx, subset.idx])
    
    # create undirected adjacency matrix from random tag matrix
    fl.graph.subset.tmp <- graph.adjacency(tags.m.subset.tmp,
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
  
  for (cluster in cluster.seq) {
      tune_grid$mod[tune_grid$step == step & 
                      tune_grid$cluster == cluster] <- 
        median(modularity.m[,cluster])
      tune_grid$mod_q025[tune_grid$step == step & 
                           tune_grid$cluster == cluster] <- 
        quantile(modularity.m[,cluster], 0.025)
      tune_grid$mod_q975[tune_grid$step == step & 
                           tune_grid$cluster == cluster] <- 
        quantile(modularity.m[,cluster], 0.975)
      tune_grid$mod_se[tune_grid$step == step & 
                         tune_grid$cluster == cluster] <- 
        sd(modularity.m[,cluster]/sqrt(nrow(tune_grid)))
  }
  modularity.l[[step]] <- modularity.m
  cat(paste0(step, " steps!", " (", 
             round((((step-min(step.seq))+1)/length(step.seq))*100,2), "%)"))
}

## (3) summarize MC modularity statistics

## top 10 step and cluster sizes based on modularity
tune_grid %>% 
  dplyr::arrange(desc(mod)) %>%
  filter(cluster != 1) %>%
  head(10)

## optimal tuning parameters
opt.mod <- tune_grid[which.max(tune_grid$mod),]

## plot full 2D modularity grid
ggplot(tune_grid, aes(x = cluster, y = step, z = mod, fill = mod)) +
  geom_tile() + 
  geom_tile(data = opt.mod, color = "red") +
  scale_x_continuous(breaks = seq(2, max(cluster.seq), 2), 
                     expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(min(step.seq), max(step.seq), 2), 
                     expand = c(0, 0)) +
  theme_bw() + mythemes +
  labs(x = "clusters (k)", y = "steps (n)", fill = "modularity (Q)")

ggsave("figures/modularity_tuning_grid.png")

## modularity scores per cluster sizes for different steps
ggplot(tune_grid, aes(x = cluster, y = mod)) +
  geom_pointrange(aes(ymin = mod_q025, ymax = mod_q975), col = "grey25") +
  facet_wrap(~step) + theme_light() + mythemes

ggsave("figures/modularity_tuning_facets.png")

## modularity scores w/ outer quantiles for different cluster sizes, using optimal step size
## suggests that you don't necessarily get a huge benefit from more than 9 clusters
tune_grid %>%
  filter(step == opt.mod$step) %>%
  ggplot(aes(x = cluster, y = mod)) +
  geom_hline(aes(yintercept = opt.mod$mod_q025), lty = "dashed") +
  geom_errorbar(aes(ymin = mod_q025, ymax = mod_q975), width = 0) +
  geom_point(shape = 21, size = 2, fill = "white") + 
  theme_bw() + mythemes

ggsave("figures/modularity_quantiles_optimal_steps.png")

## boxplot of modularity for different cluster sizes, using optimal step size

png("figures/modularity_boxplot_optimal_steps.png", 
    width=12, height=8, units='in', res=300)

boxplot(modularity.l[[opt.mod$step]],
        type="l", xlab= "clusters (k)", 
        ylab= "modularity (Q)", 
        main = paste("Modularity changing with k (Walktrap, steps = ", 
                     opt.mod$step, ")"))

dev.off()

## (4) run full cluster analysis on full dataset of tags using optimal number of steps and clusters
# create undirected adjacency matrix from full tag matrix
fl.graph <- graph.adjacency(tags.m,
                            weighted=TRUE, 
                            mode="undirected",
                            diag=TRUE)

# clustering from full matrix
cw <- cluster_walktrap(fl.graph, 
                       weights = E(fl.graph)$weight,
                       membership = T,
                       steps = opt.mod$step)

cw.cut <- cut_at(cw, no = 10)
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

summary(dg$res)

pdf("figures/centrality_eigenvalue.pdf", width = 12, height = 8)
barplot(sort(ec, decreasing = T)[1:30], las=2)
dev.off()

pdf("figures/centrality_closeness.pdf", width = 12, height = 8)
barplot(sort(cl, decreasing = T)[1:30], las=2)
dev.off()

pdf("figures/centrality_betweenness.pdf", width = 12, height = 8)
barplot(sort(bt, decreasing = T)[1:30], las=2)
dev.off()

tag.ranks <- matrix(nrow = 10, ncol = 10)
for (clust in 1:10) {
  cw.comm %>%
    filter(cluster == clust) %>%
    arrange(desc(eigen_centrality)) -> cw.eigen.rank
  if (nrow(cw.eigen.rank) < 10) {
    tag.ranks[,clust] <- c(cw.eigen.rank$tag, rep(NA, (10-nrow(cw.eigen.rank))))
  } else {
    cw.eigen.rank %>%
      slice(1:10) %>%
      dplyr::select(tag) %>% pull() -> tag.ranks[,clust]
  }
}

tag.ranks

# ASSIGN CLUSTERS TO IMAGES BASED ON TAG CONTENT -------------------------------

## assign themes to photos based on tags, randomly split ties
for (i in 1:nrow(flickr_rur_tidy)) {
  pic <- flickr_rur_tidy[i,]
  tag.seq <- c(pic$tags1, pic$tags2, pic$tags3, pic$tags4, pic$tags5, 
               pic$tags6, pic$tags7, pic$tags8, pic$tags9, pic$tags10)
  # replace empty tags with NAs and remove from vector
  tag.seq[tag.seq == " "] <- NA
  tag.seq <- na.omit(tag.seq)
  # match clusters with tags
  match.theme <- function(x) {
    if (x %in% cw.comm$tag)
      tag.theme <- cw.comm$cluster[cw.comm$tag == x]
    else 
      tag.theme <- "other"
  }
  theme.seq <- unlist(lapply(tag.seq, match.theme))
  # identify the dominant cluster (split ties randomly)
  photo.theme <- Mode(na.omit(theme.seq[theme.seq != "other"]))
  flickr_rur_tidy$cluster[i] <- photo.theme
}

barplot(table(flickr_rur_tidy$cluster))
