#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R CODE FOR FLICKR-NFR CLUSTER ANALYSIS
# Harrison B Goldspiel | hbgoldspiel@gmail.com
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREPARE IMAGE DATA FOR CLUSTER ANALYSIS ------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load raw Flickr data
flickr_all <- read.csv("data/flickr_NFR_all.csv") 
flickr_public <- read.csv("data/flickr_NFR_public_land.csv")

dim(flickr_all)
# [1] 280509     25
dim(flickr_public)
# [1] 61296    21

# omit corrupted data (six rows with misplaced or missing column values)
flickr_all <- flickr_all[!is.na(flickr_all$OBJECTID_12),]
# omit non-rural images
flickr_rural <- flickr_all[flickr_all$isUrban == 0,]

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
         month = month(date),
         hour = hour(datetime),
         id = as.factor(id)) %>%
  distinct(id, .keep_all = TRUE)

flickr_public_tidy <- flickr_public %>%
  mutate(datetime = ymd_hms(datetaken),
         date = date(datetime),
         month = month(date),
         hour = hour(datetime),
         id = as.factor(id)) %>%
  distinct(id, .keep_all = TRUE)

dim(flickr_rural_tidy)
# [1] 149192     29
dim(flickr_public_tidy)
# [1] 22325    25

flickr_rural_tidy_tagged <- flickr_rural_tidy %>%
  filter(tags1 %notin% " " & tags2 %notin% " " & tags3 %notin% " " &
             tags4 %notin% " " & tags5 %notin% " " & tags6 %notin% " " &
             tags7 %notin% " " & tags8 %notin% " " & tags9 %notin% " " &
             tags10 %notin% " ")
flickr_public_tidy_tagged <- flickr_public_tidy %>%
  filter(tags1 %notin% " " & tags2 %notin% " " & tags3 %notin% " " &
             tags4 %notin% " " & tags5 %notin% " " & tags6 %notin% " " &
             tags7 %notin% " " & tags8 %notin% " " & tags9 %notin% " " &
             tags10 %notin% " ")

# rank tags by percentile, select tags past certain threshold for clustering
rural_tags <- 
  list(flickr_rural_tidy_tagged$tags1, flickr_rural_tidy_tagged$tags2, 
       flickr_rural_tidy_tagged$tags3, flickr_rural_tidy_tagged$tags4, 
       flickr_rural_tidy_tagged$tags5, flickr_rural_tidy_tagged$tags6,
       flickr_rural_tidy_tagged$tags7, flickr_rural_tidy_tagged$tags8, 
       flickr_rural_tidy_tagged$tags9, flickr_rural_tidy_tagged$tags10)

public_tags <- 
  list(flickr_public_tidy_tagged$tags1, flickr_public_tidy_tagged$tags2, 
       flickr_public_tidy_tagged$tags3, flickr_public_tidy_tagged$tags4, 
       flickr_public_tidy_tagged$tags5, flickr_public_tidy_tagged$tags6,
       flickr_public_tidy_tagged$tags7, flickr_public_tidy_tagged$tags8, 
       flickr_public_tidy_tagged$tags9, flickr_public_tidy_tagged$tags10)

rural_tag_freq <- as.data.frame(table(unlist(rural_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / length(rural_tags[[1]]),
         pct = percentile_rank(freq)) 
  
public_tag_freq <- as.data.frame(table(unlist(public_tags))) %>%
  dplyr::select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / length(public_tags[[1]]),
         pct = percentile_rank(freq))

write.csv(rural_tag_freq, "data/rural_tags.csv")
write.csv(public_tag_freq, "data/public_tags.csv")

# list of uncommon (N < 5 images) tags
rare_tags <- as.character(rural_tag_freq$tag[rural_tag_freq$freq < 5])

# word cloud of frequently (i.e., N >= 5 images)  occurring tags in rural areas
png("figures/flickr_rural_toptags.png", 
    width=12, height=8, units='in', res=300)
set.seed(140)   
wordcloud(words = rural_tag_freq$tag, freq = rural_tag_freq$freq, 
          min.freq = 5, max.words = Inf, scale=c(8,0.2), 
          colors=brewer.pal(6, "Dark2"), random.order=FALSE, 
          rot.per=0.15)  
par(mar = rep(0, 4))
dev.off()

# word cloud of frequently (i.e., N >= 5 images)  occurring tags in public lands
png("figures/flickr_public_toptags.png", 
    width=12, height=8, units='in', res=300)
set.seed(140)   
wordcloud(words = public_tag_freq$tag, freq = public_tag_freq$freq, 
          min.freq = 5, max.words = Inf, scale=c(8,0.2), 
          colors=brewer.pal(6, "Dark2"), random.order=FALSE, 
          rot.per=0.15)  
par(mar = rep(0, 4))
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN MONTE CARLO (MC) CLUSTER ANALYSIS for rural and public land images------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## (1) create raw rural tag matrix
## (2) random sample 80% of rural tags over 1000 iterations, without replacement, 
## (3) and for each iteration, create tag co-occurence matrix and use Walktrap hierchical clustering algorithm, connecting clusters with Ward's distance, saving the modularity and cluster size
## (4) use those MC simulated modularities to pick optimal number of clusters (K);
## (5) run cluster analysis on full dataset of tags, using K clusters

## (1) create original tag matrix
rural_network_in <- 
  flickr_rural_tidy_tagged %>%
  dplyr::select(id, tags1, tags2, tags3, tags4, tags5, 
                tags6, tags7, tags8, tags9, tags10) %>%
  melt(id.vars = "id", value.name = "tag") %>%
  # filter to omit rare tags (N < 5)
  filter(tag %in% rural_tag_freq$tag[rural_tag_freq$freq >= 5] & tag != "") %>%
  # filter to omit overly common tags 
  filter(tag %in% rural_tag_freq$tag[rural_tag_freq$prop < 0.10]) %>%
  left_join(rural_tag_freq, by = "tag") %>%
  dplyr::select(id, tag, freq)

public_network_in <- 
  flickr_public_tidy_tagged %>%
  dplyr::select(id, tags1, tags2, tags3, tags4, tags5, 
                tags6, tags7, tags8, tags9, tags10) %>%
  melt(id.vars = "id", value.name = "tag") %>%
  # filter to omit rare tags (N < 5)
  filter(tag %in% public_tag_freq$tag[public_tag_freq$freq >= 5] & tag != "") %>%
  # filter to omit overly common tags
  filter(tag %in% public_tag_freq$tag[public_tag_freq$prop < 0.10]) %>%
  left_join(public_tag_freq, by = "tag") %>%
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
  geom_tile(data = opt.mod, color = "red") +
  scale_x_continuous(breaks = seq(2, max(rural_mcmc_wt$cluster_seq), 2), 
                     expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(min(rural_mcmc_wt$step_seq), 
                                  max(rural_mcmc_wt$step_seq), 2), 
                     expand = c(0, 0)) +
  theme_bw() + mythemes +
  labs(x = "clusters (k)", y = "steps (n)", fill = "modularity (Q)")

ggsave("figures/modularity_tuning_grid.png")

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

tag.ranks
write.csv(tag.ranks, "data/rural_tag_clusters.csv", na = "", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ASSIGN CLUSTERS TO IMAGES BASED ON TAG CONTENT -----------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## aggregate clusters into core groups

## 1:  arts (1, 16, 27)                              [tech, music, photography]
## 2:  sports (18, 21, 23, 26, 28, 29)*              [team sports, boxing, shooting, billiards]
## 3:  scenery (2)**                                 [megacluster, mixture of natural and cultural landscape features]
## 4:  food/dining (3, 6, 8, 14, 15, 22)             [food, dessert, dining, alcohol] 
## 5:  aquatic recreation (4, 5)                     [watercrafts, fishing, and watersports]
## 6:  natural life (9, 11, 17, 19, 25)**            [flora, fauna, and fungi (and some livestock and zoo animals)]
## 7:  equestrian (10)                               [horses]
## 8:  people (12)**                                 [people]
## 9:  pets (13)                                     [dogs and cats]
## 10: transportation (20, 24)                       [cars, trucks, buses, trains, bicycles, aircraft]

cw.comm <- cw.comm %>%
  mutate(theme = case_when(cluster %in% c(1,16,27) ~ "arts",
                           cluster %in% c(18,21,23,26,28,29) ~ "sports",
                           cluster == 2 ~ "scenery",
                           cluster %in% c(3,6,8,14,15,22) ~ "food",
                           cluster %in% c(4,5) ~ "aquatics",
                           cluster %in% c(9,11,17,19,25) ~ "natural life",
                           cluster == 10 ~ "equestrian",
                           cluster == 12 ~ "people",
                           cluster == 13 ~ "dogs/cats",
                           cluster %in% c(20,24) ~ "transportation"))


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

rural_photos_clustered <- assign.theme(photos = flickr_rural_tidy_tagged, 
                                       clusters = cw.comm)

# quick bar plot of cluster composition (mirror)
rural.themes <-
  rural_photos_clustered %>%
  mutate(public = ifelse(id %in% flickr_public$id, "public", "private")) %>%
  group_by(public, theme) %>%
  summarize(count = n()) %>%
  na.omit() %>%
  ungroup() %>%
  group_by(public) %>%
  mutate(prop = round(count / sum(count), 3)) %>%
  ungroup()

rural.themes.plot <-
  ggplot(rural.themes, aes(x = reorder(theme, count), y = prop)) +
  facet_wrap(~public, nrow = 2) +
  geom_col() +
  geom_text(aes(label = paste("n = ", count), y = prop + 0.1)) +
  geom_text(aes(label = prop, y = prop + 0.2)) +
  labs(x = "Theme", y = "Proportion of images") +
  theme_bw() + 
  labs(title = "Rural image themes")

rural.themes.plot
ggsave("figures/rural_themes.png", width = 10, height = 7, dpi = 600)


save.image(file = "data/cluster_analysis_output.RData")
