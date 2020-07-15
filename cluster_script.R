# R CODE FOR FLICKR-NFR CLUSTER ANALYSIS

source("custom_functions_settings.R")

# Analytic summary:
## 1. Mine Flickr data for NFR fro 2012-2016 from open-source API (DONE)
## 2. Run Flickr images through Clarifai image classification algorith (DONE)
## 3. Perform cluster analysis on rural images to identify core themes of 
##### human engagement in rural parts of the NFR
## 4. Assign themes to images based on tags
## 5. Visualize hotspots of images (PUD) (DONE for overall imagery)
## 6. Use random forest & boruta to identify core drivers of different themes
##### (do this for each theme/state and theme/NFR total)
## 7. Look at partial effects or GAMs of the different themes

# load packages
library(maptools)
library(spatstat)
library(ggmap)
library(rgdal)
library(tidyverse)
library(tidytext)
library(wordcloud)
library(igraph)
library(RColorBrewer)

# load raw Flickr data
flickr_all <- read.csv("data/flickr_NFR_all.csv") 

# omit non-rural images
flickr_rural <- flickr_all[flickr_all$isUrban == 0,]

# omit images without any tags
flickr_rural %>%
  filter(!(tags1 %in% " " & tags2 %in% " " & tags3 %in% " " &
             tags4 %in% " " & tags5 %in% " " & tags6 %in% " " &
             tags7 %in% " " & tags8 %in% " " & tags9 %in% " " &
             tags10 %in% " ")) -> flickr_rural

n_images <- nrow(flickr_rural)

# rank tags by percentile, select tags past certain threshold for clustering
tags <- list(flickr_rural$tags1, flickr_rural$tags2, flickr_rural$tags3, 
             flickr_rural$tags4, flickr_rural$tags5, flickr_rural$tags6,
             flickr_rural$tags7, flickr_rural$tags8, flickr_rural$tags9,
             flickr_rural$tags10)

tag_freq <- as.data.frame(table(unlist(tags))) %>%
  select(tag = Var1, freq = Freq) %>%
  filter(tag %notin% " ") %>%
  arrange(desc(freq)) %>%
  mutate(prop = freq / n_images,
         pct = p_rank(freq))

write.csv(tag_freq, "data/rural_tags.csv")

# word cloud of frequently occuring (N >= 50) tags
png("figures/flickr_rural_top100tags.png", width=12, height=8, units='in', res=300)
set.seed(140)   
wordcloud(words = tag_freq$tag, freq = tag_freq$freq, 
          min.freq = 50, max.words = Inf, scale=c(8,0.2), 
          colors=brewer.pal(6, "Dark2"), random.order=FALSE, 
          rot.per=0.15)  
par(mar = rep(0, 4))
dev.off()

# run MC hierarchical cluster analysis:
## (1) random sample 80% of tags over 10000 iterations, without replacement, 
## (2) and for each iteration, create tag co-occurence matrix and
## use Walktrap hierchical clustering algorithm, connecting clusters with Ward's 
## distance, saving the modularity and cluster size
## (3) use those MC simulated modularities to pick optimal number of clusters (K);
## (4) run cluster analysis on full dataset of tags, using K clusters




# save cluster types to tags

# assign themes to photos based on tags, randomly split ties



