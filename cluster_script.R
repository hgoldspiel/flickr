#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R CODE FOR FLICKR-NFR CLUSTER ANALYSIS
# Harrison B Goldspiel | hbgoldspiel@gmail.com
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("custom_functions_settings.R")

# Analytic summary:
## 1. Mine Flickr data for NFR fro 2012-2016 from open-source API (DONE)
## 2. Run Flickr images through Clarifai image classification algorith (DONE)
## 3. Perform cluster analysis on rural images to identify core themes of human engagement in rural parts of the NFR
## 4. Assign themes to images based on tags
## 5. Visualize hotspots of images (PUD) (DONE for overall rural imagery)
## 6. Use random forest & Boruta to identify core drivers of different themes (do this for each theme/state and theme/NFR total)
## 7. Look at partial effects or GAMs of the different themes

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
library(Rlda)
library(labdsv)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREPARE IMAGE DATA FOR CLUSTER ANALYSIS ------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load raw Flickr data
flickr_all <- read.csv("data/flickr_NFR_all.csv") 
flickr_public <- read.csv("data/flickr_NFR_public_land.csv")

dim(flickr_all)
# [1] 280509     25
dim(flickr_public)
# [1] 28461    21

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
# Obtain sample of 500 images for exploring numerical ecology approaches -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(131)
rural_sample <- sample_n(flickr_rural_tidy_tagged, 500)
write.csv(rural_sample, "data/rural_sample.csv", na = "", row.names = FALSE)

# convert sample dataset to community matrix with rows as photos and columns as tags
## create uniform levels for tags in the sample
tag_levels <- unique(c(rural_sample$tags1, rural_sample$tags2, rural_sample$tags3,
                       rural_sample$tags4, rural_sample$tags5, rural_sample$tags6,
                       rural_sample$tags7, rural_sample$tags8, rural_sample$tags9,
                       rural_sample$tags10))

tag_levels_all <- unique(c(flickr_rural_tidy_tagged$tags1, flickr_rural_tidy_tagged$tags2, 
                            flickr_rural_tidy_tagged$tags3, flickr_rural_tidy_tagged$tags4, 
                            flickr_rural_tidy_tagged$tags5, flickr_rural_tidy_tagged$tags6,
                            flickr_rural_tidy_tagged$tags7, flickr_rural_tidy_tagged$tags8, 
                            flickr_rural_tidy_tagged$tags9, flickr_rural_tidy_tagged$tags10))

## assign same levels to all tag columns
#tag.lev <- function(x) factor(as.factor(x), levels = tag_levels)
#rural_sample[,13:22] <- lapply(rural_sample[,13:22], tag.lev)

# convert to wide-form "community" matrix with each tag as its own binary col
resorted <- rural_sample[order(rural_sample$id),]
tagcols <- c("tags1", "tags2", "tags3", "tags4", "tags5", 
             "tags6", "tags7", "tags8", "tags9", "tags10")
suppressWarnings(
  rural_sample_mat <- 
    rural_sample %>%
    melt(id.vars = c("id", tagcols)) %>%
    dcast(id ~ c(tags1, tags2, tags3, tags4, tags5, 
                 tags6, tags7, tags8, tags9, tags10)) %>%
    dplyr::select(-c("id", " ")) %>%
    mutate_each(funs(replace(., . == 18, 1))) %>%
    mutate(user = resorted$owner, 
           id = resorted$id,
           state = resorted$STUSPS,
           date = resorted$date) %>%
    dplyr::select(id, user, state, date, c(1:length(tag_levels))))

## do same for full dataset
# tag.lev.full <- function(x) factor(as.factor(x), levels = tag_levels_full)
# flickr_rural_tidy[,13:22] <- lapply(flickr_rural_tidy[,13:22], tag.lev.full)

# convert to wide-form "community" matrix with each tag as its own binary col
resorted.tidy <- flickr_rural_tidy_tagged[order(flickr_rural_tidy_tagged$id),]
remove.rare.tags <- function(x) {x = ifelse(x %in% rare_tags, " ", x)}
suppressWarnings(
  rural_mat <-
    flickr_rural_tidy_tagged %>% 
    filter(tags1 %notin% " ") %>%
    mutate_at(c(tagcols), remove.rare.tags) %>%
    melt(id.vars = c("id", tagcols)) %>%
    dcast(id ~ c(tags1, tags2, tags3, tags4, tags5, 
                 tags6, tags7, tags8, tags9, tags10)) %>%
    dplyr::select(-c("id", " ")) %>%
    mutate_each(funs(replace(., . == 18, 1))) %>%
    mutate(user = resorted.tidy$owner, 
           id = resorted.tidy$id,
           state = resorted.tidy$STUSPS,
           date = resorted.tidy$date) %>%
    dplyr::select(id, user, state, date, c(1:length(tag_levels_all))))

# write new output csv file
write.csv(rural_sample_mat, "data/rural_sample_mat.csv", na = "", row.names = FALSE)
write.csv(rural_mat, "data/rural_mat.csv", na = "", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Try to break down images by hierarchical clustering approaches -------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# run different cluster algorithms on distance matrix
tag.mat <- rural_sample_mat[,-c(1:4)]
D.jac <- vegdist(tag.mat, "jaccard", binary = TRUE)
attr(D.jac, "labels") <- rural_sample_mat$id
D.jac.single <- hclust(D.jac, method = "single") # single linkage
D.jac.complete <- hclust(D.jac, method = "complete") # complete linkage
D.jac.UPGMA <- hclust(D.jac, method = "average") # UPGMA agglomerate
D.jac.ward <- hclust(D.jac, method = "ward.D2") # Ward's minimum variance

# the average (agglomerate, UPGMA) clustering method looks best (but it always will be, mathematically)
# look for interpretable clusters by graphing "fusion levels" (max 20 clusters)
par(mfrow = c(1,1))
plot(
  D.jac.UPGMA$height[481:499],
  20:2,
  type = "S",
  main = "Fusion levels - Chord - UPGMA",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(D.jac.UPGMA$height[481:499],
     20:2,
     20:2,
     col = "red",
     cex = 0.8) # looks like 10 clusters is best

# Now look at a silhouette plot for 10 clusters
# (this shows how well individual tags belong to their clusters)
# (more positive values indicate better coherence)
k = 10
cutc = cutree(D.jac.UPGMA,k=k)
sil = silhouette(cutc, D.jac)
sil.ord = sortSilhouette(sil)
rownames(sil.ord) = row.names(tag.mat)[attr(sil.ord,'iOrd')]
plot(sil.ord,main='Silhouette plot, UPGMA (Jaccard)',
     col=sil.ord[,1]+1,
     nmax.lab=100)

sil[1:500,1:3] %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  summarize(avg_width = mean(sil_width),
            n_k = n()) %>%
  ggplot(aes(x = as.factor(cluster), y = avg_width, label = n_k)) + 
  geom_col() + 
  geom_text(vjust = -0.15) + 
  theme_bw()
  
# most groups don't have a lot of coherence (clusters dominated by single few tags)
# average silhouette width is only 0.06


# what tags actually belong to each photo in the clusters
kdat <- rural_sample %>%
  mutate(k = cutc) %>%
  dplyr::select(k, tags1, tags2, tags3, tags4, tags5, 
                tags6, tags7, tags8, tags9, tags10, url) 

# plot top ten tags by frequency appearing in each cluster
tiff("figures/top_ten_tidy2_k10.tiff", width = 7, height = 6, 
     units = "in", res = 600, compression = "lzw")

par(mfrow = c(3,3))
for(k in c(1:5,7,9)){
  knew = kdat[kdat$k == k,]
  tags.l <- na.omit(list(knew[c("tags1", "tags2", "tags3", "tags4", "tags5", 
                              "tags6", "tags7", "tags8", "tags9", "tags10")]))
  tags.df <- as.data.frame(table(unlist(tags.l)))
  colnames(tags.df) <- c("tag", "freq")
  tags.df <- tags.df[order(tags.df$freq, decreasing = TRUE),][tags.df$tag != " ",]
  tags.df$prop <- tags.df$freq/nrow(knew)
  bp <- barplot(tags.df$prop[1:10], names.arg = "", 
                 ylim = c(0,max(tags.df$prop)), xlab = "", 
                 ylab = "Frequency", 
                main = paste("cluster ", k, " ( n = ", nrow(knew), ")"))
  text(bp, par("usr")[3], labels = tags.df$tag[1:10], 
       srt = 35, adj = c(1.1,1.1), xpd = TRUE)
  axis(2)
}

dev.off()

# use indicator value index to identify species specificity for each cluster
iva <- indval(tag.mat, cutc, numitr = 10000)

# Correction of the p-values for multiple testing
pval.adj <- p.adjust(iva$pval)

# Table of the significant indicator species
gr <- iva$maxcls[pval.adj <= 0.05]
iv <- iva$indcls[pval.adj <= 0.05]
pv <- iva$pval[pval.adj <= 0.05]
fr <- apply(tag.mat > 0, 2, sum)[pval.adj <= 0.05]
fidg <- data.frame(
  group = gr,
  indval = iv,
  pvalue = pv,
  freq = fr
)
fidg <- fidg[order(fidg$group, -fidg$indval), ]
fidg
# Export the result to a CSV file (to be opened in a spreadsheet)
write.csv(fidg, "data/cluster/IndVal-dfs.csv")




# Try silhouette function to find optimal number of clusters (up to 20)
sil.wid <- numeric(19)
# Calculate silhouette widths for each number of clusters,
# disregarding the trivial k = 1:
for(k in 2:20){
  tmp = silhouette(cutree(D.jac.UPGMA,k=k),D.jac)
  sil.wid[k] = summary(tmp)$avg.width
  }
# Best width
k.best = which.max(sil.wid)
# Plotting:
par(xpd=NA)
plot(1:20,sil.wid,type='h',
     main='Silhouette: optimal number of clusters, UPGMA',
     xlab='k number of clusters',ylab='Average silhouette width',cex.lab=1.25)
lines(rep(k.best,2),c(0,max(sil.wid)),col=2,cex=1.5,lwd=3)

# Silhouette plot of the 'optimal' partition:
k = k.best
cutc = cutree(D.jac.UPGMA,k=k)
sil = silhouette(cutc,D.jac)
sil.ord = sortSilhouette(sil)
rownames(sil.ord) = row.names(tag.mat)[attr(sil.ord,'iOrd')]
plot(sil.ord,main='Silhouette plot, UPGMA (Jaccard)',cex.names=0.8,
     col=sil.ord[,1]+1,nmax.lab=100, cex = 0.5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bayesian LDA for binary data to find clusters -----------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# try Bayesian LDA for binary data from the Rlda package
gamma <- 0.01
alpha0 <- 0.01
alpha1 <- 0.01
tag.mat$loc.id<-seq(1,nrow(tag.mat))
lda.out <- rlda.fastbernoulli(tag.mat,
                              loc.id = "loc.id",
                              n_community = 20, 
                              alpha0 = 0.01, 
                              alpha1 = 0.01,
                              gamma = 0.01, 
                              n_gibbs = 1000, 
                              ll_prior = TRUE,
                              display_progress = TRUE)

ll <- lda.out$logLikelihood
plot(ll, type = "l", xlab = "Iterations", ylab = "Log(likel.) + log(prior)")
abline(v = 500, col = 'grey')

summary(lda.out, burnin = 0.5)
plot(lda.out, burnin= 0.5)
Theta <- summary(lda.out, burnin= 0.5, silent = TRUE)$Theta


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN MONTE CARLO (MC) CLUSTER ANALYSIS---------------------------------------
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
  filter(tag %in% rural_tag_freq$tag[rural_tag_freq$freq >= 5]) %>%
  # filter to omit overly common tags 
  #filter(tag %in% rural_tag_freq$tag[rural_tag_freq$prop < 0.10]) %>%
  left_join(rural_tag_freq, by = "tag") %>%
  dplyr::select(id, tag, freq)

# create co-occurrence matrix
tags.dt <- as.data.frame(crossprod(table(rural_network_in[1:2])))

# reduce matrix appropriately for undirected cluster analysis
tags.v <- colnames(tags.dt)
tags.m <- (as.matrix(data.frame(tags.dt)))
dimnames(tags.m) <- list(tags.v, tags.v)
total.occur <- colSums(tags.m)
tags.m[lower.tri(tags.m, diag=T)] <- 0

## (2-3) run MC cluster algorithm to obtain robust modularity plot, iterating over a range of random step sizes
n.iter <- 100 # n MC iterations
step.seq <- 3:8 # step range (range recommended by Pons & Lapaty (2005))
cluster.seq <- 1:30 # cluster range

tune_grid <- expand.grid(
  step       = step.seq,
  cluster    = cluster.seq,
  mod        = NA,
  mod_q025   = NA,
  mod_q975   = NA,
  mod_se     = NA)

modularity.l <- list()
set.seed(2000)
for (step in step.seq) {
  modularity.m <- matrix(NA, nrow = n.iter, ncol = length(cluster.seq))
  subset.ratio <- 0.8 # 80% of the tag dataset
  for (iter.idx in 1:n.iter) {
    n.tags <- nrow(tags.m)
    subset.idx <- (sample(1:n.tags, size = floor(n.tags*subset.ratio), 
                          replace = F))
    
    # getting random subset of rural tag matrix
    tags.m.subset.tmp <- (tags.m[subset.idx, subset.idx])
    
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
  
  for (cluster in cluster.seq) {
      tune_grid$mod[
        tune_grid$step == step & 
          tune_grid$cluster == cluster] <- median(modularity.m[,cluster])
      tune_grid$mod_q025[
        tune_grid$step == step &
          tune_grid$cluster == cluster] <- quantile(modularity.m[,cluster], 0.025)
      tune_grid$mod_q975[
        tune_grid$step == step & 
          tune_grid$cluster == cluster] <- quantile(modularity.m[,cluster], 0.975)
      tune_grid$mod_se[
        tune_grid$step == step & 
          tune_grid$cluster == cluster] <- sd(modularity.m[,cluster]/
                                                sqrt(nrow(tune_grid)))
  }
  modularity.l[[step]] <- modularity.m
  cat(paste0(step, " steps!", " (", 
             round((((step-min(step.seq))+1)/length(step.seq))*100,2), "%)"))
} ; beep(3)

## (3) summarize MC modularity statistics

# top 10 step and cluster sizes based on modularity
tune_grid %>% 
  dplyr::arrange(desc(mod)) %>%
  filter(cluster != 1) %>%
  head(10)

# optimal tuning parameters
opt.mod <- tune_grid[which.max(tune_grid$mod),]

# plot full 2D modularity grid
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

# modularity scores per cluster sizes for different steps
ggplot(tune_grid, aes(x = cluster, y = mod)) +
  geom_pointrange(aes(ymin = mod_q025, ymax = mod_q975), col = "grey25") +
  facet_wrap(~step) + theme_light() + mythemes

ggsave("figures/modularity_tuning_facets.png")

# modularity scores w/ outer quantiles for different cluster sizes, using optimal step size
# suggests that you don't necessarily get a huge benefit from more than 9 clusters
tune_grid %>%
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

boxplot(modularity.l[[opt.mod$step]],
        type="l", xlab= "clusters (k)", 
        ylab= "modularity (Q)", 
        main = paste("Modularity changing with k (Walktrap, steps = ", 
                     opt.mod$step, ")"))

dev.off()

## (4) run full cluster analysis on full dataset of tags using optimal number of steps and clusters
# create undirected adjacency matrix from full tag matrix
fl.graph <- graph_from_adjacency_matrix(tags.m,
                                        weighted=TRUE, 
                                        mode="undirected",
                                        diag=TRUE)

# clustering from full matrix
set.seed(666)
cw <- cluster_walktrap(fl.graph, 
                       weights = E(fl.graph)$weight,
                       membership = T,
                       steps = opt.mod$step)

# optimal number of clusters (based on 95% quantile range of top modularity)
nclust <- min(tune_grid$cluster[
  tune_grid$step == opt.mod$step & tune_grid$mod_q975 >= 
    tune_grid[which.max(tune_grid$mod),]$mod_q025])

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

# ASSIGN CLUSTERS TO IMAGES BASED ON TAG CONTENT -------------------------------

## aggregate clusters into core groups

## 1:  arts/sports (7, 10, 16, 22, 24, 26, 29)*      [tech, music, sports, photography]
## 2:  scenery (2)**                                 [megacluster, mixture of natural and cultural landscape features]
## 3:  food/dining (5, 6, 9, 13, 17, 20, 23, 25, 27) [food, dessert, dining, alcohol] 
## 4:  fishing (16)                                  [watercrafts, fishing, and watersports]
## 5:  natural life (3, 4, 18, 19, 21)**             [flora, fauna, and fungi (and some livestock and zoo animals)]
## 6:  equestrian (8)                                [horses]
## 7:  people (1)**                                  [people]
## 8:  pets (15)                                     [dogs and cats]
## 9:  transportation (11, 28)                       [cars, trucks, buses, trains, bicycles, aircraft]
## 10: firearms and shooting (12, 14, 30)            [guns, hunting, and shooting ranges]

cw.comm <- cw.comm %>%
  mutate(theme = case_when(cluster %in% c(7,10,16,22,24,26,29) ~ "arts/sports",
                           cluster == 2 ~ "scenery",
                           cluster %in% c(5,6,9,13,17,20,23,25,27) ~ "food/dining",
                           cluster == 16 ~ "fishing",
                           cluster %in% c(3,4,18,19,21) ~ "natural life",
                           cluster == 8 ~ "equestrian",
                           cluster == 1 ~ "people",
                           cluster == 15 ~ "dogs/cats",
                           cluster %in% c(11,28) ~ "transportation",
                           cluster %in% c(12,14,30) ~ "firearms/hunting"))


## assign themes to photos based on tags, randomly split ties
assign.theme <- function(photos, clusters) {
  message("Assigning themes...")
  out <- matrix(nrow = nrow(photos), ncol = ncol(photos) + 1,
                dimnames = list(c(1:nrow(photos)),
                                c(colnames(photos), "theme")))
  for (i in 1:nrow(photos)) {
    tag.seq <- c(photos[i, c("tags1", "tags2", "tags3", "tags4", "tags5",  
                             "tags6", "tags7", "tags8","tags9", "tags10")])
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
    out[i, "theme"] <- photo.theme
  } 
  return(out)
  beep(3)
  message("done!")
}

rural_photos_clustered <- assign.theme(photos = flickr_rural_tidy_tagged, 
                                       clusters = cw.comm)

# quick bar plot of cluster composition
rural.themes.plot <-
  ggplot(as.data.frame(sort(table(rural_photos_clustered[,"theme"]))),
       aes(x = Var1, y = Freq)) +
  geom_col() + 
  geom_text(aes(label = round(Freq/sum(Freq), 3), y = Freq + 1000, angle = 90)) +
  geom_text(aes(label = paste("\ \ \ \ \ n = ", Freq), y = Freq + 5000), col = "grey50") +
  labs(x = "theme", y = "n (images)") + 
  theme_classic() + coord_flip() +
  labs(title = "Rural image themes")

rural.themes.plot




# repeat exercise above but after omitting tags that appear in over 10% of images
## (1) create original tag matrix
rural_network_in2 <- 
  flickr_rural_tidy_tagged %>%
  dplyr::select(id, tags1, tags2, tags3, tags4, tags5, 
                tags6, tags7, tags8, tags9, tags10) %>%
  melt(id.vars = "id", value.name = "tag") %>%
  # filter to omit rare tags (N < 5)
  filter(tag %in% rural_tag_freq$tag[rural_tag_freq$freq >= 5]) %>%
  # filter to omit overly common tags 
  filter(tag %in% rural_tag_freq$tag[rural_tag_freq$prop < 0.10]) %>%
  left_join(rural_tag_freq, by = "tag") %>%
  dplyr::select(id, tag, freq)

# create co-occurrence matrix
tags.dt2 <- as.data.frame(crossprod(table(rural_network_in2[1:2])))

# reduce matrix appropriately for undirected cluster analysis
tags.v2 <- colnames(tags.dt2)
tags.m2 <- (as.matrix(data.frame(tags.dt2)))
dimnames(tags.m2) <- list(tags.v2, tags.v2)
total.occur <- colSums(tags.m2)
tags.m2[lower.tri(tags.m2, diag=T)] <- 0

## (2-3) run MC cluster algorithm to obtain robust modularity plot, iterating over a range of random step sizes
n.iter <- 100 # n MC iterations
step.seq <- 3:8 # step range (range recommended by Pons & Lapaty (2005))
cluster.seq <- 1:30 # cluster range

tune_grid2 <- expand.grid(
  step       = step.seq,
  cluster    = cluster.seq,
  mod        = NA,
  mod_q025   = NA,
  mod_q975   = NA,
  mod_se     = NA)

modularity.l2 <- list()
set.seed(2000)
for (step in step.seq) {
  modularity.m2 <- matrix(NA, nrow = n.iter, ncol = length(cluster.seq))
  subset.ratio <- 0.8 # 80% of the tag dataset
  for (iter.idx in 1:n.iter) {
    n.tags <- nrow(tags.m2)
    subset.idx <- (sample(1:n.tags, size = floor(n.tags*subset.ratio), 
                          replace = F))
    
    # getting random subset of rural tag matrix
    tags.m.subset.tmp <- (tags.m2[subset.idx, subset.idx])
    
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
    
    modularity.m2[iter.idx, ] <- modularity.tmp.v
  }
  
  for (cluster in cluster.seq) {
    tune_grid2$mod[
      tune_grid2$step == step & 
        tune_grid2$cluster == cluster] <- median(modularity.m2[,cluster])
    tune_grid2$mod_q025[
      tune_grid2$step == step &
        tune_grid2$cluster == cluster] <- quantile(modularity.m2[,cluster], 0.025)
    tune_grid2$mod_q975[
      tune_grid2$step == step & 
        tune_grid2$cluster == cluster] <- quantile(modularity.m2[,cluster], 0.975)
    tune_grid2$mod_se[
      tune_grid2$step == step & 
        tune_grid2$cluster == cluster] <- sd(modularity.m2[,cluster]/
                                              sqrt(nrow(tune_grid2)))
  }
  modularity.l2[[step]] <- modularity.m2
  cat(paste0(step, " steps!", " (", 
             round((((step-min(step.seq))+1)/length(step.seq))*100,2), "%)"))
} ; beep(3)

## (3) summarize MC modularity statistics

# top 10 step and cluster sizes based on modularity
tune_grid2 %>% 
  dplyr::arrange(desc(mod)) %>%
  filter(cluster != 1) %>%
  head(10)

# optimal tuning parameters
opt.mod2 <- tune_grid2[which.max(tune_grid2$mod),]

# plot full 2D modularity grid
ggplot(tune_grid2, aes(x = cluster, y = step, z = mod, fill = mod)) +
  geom_tile() + 
  geom_tile(data = opt.mod, color = "red") +
  scale_x_continuous(breaks = seq(2, max(cluster.seq), 2), 
                     expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(min(step.seq), max(step.seq), 2), 
                     expand = c(0, 0)) +
  theme_bw() + mythemes +
  labs(x = "clusters (k)", y = "steps (n)", fill = "modularity (Q)")

ggsave("figures/modularity_tuning_grid.png")

# modularity scores per cluster sizes for different steps
ggplot(tune_grid2, aes(x = cluster, y = mod)) +
  geom_pointrange(aes(ymin = mod_q025, ymax = mod_q975), col = "grey25") +
  facet_wrap(~step) + theme_light() + mythemes

ggsave("figures/modularity_tuning_facets.png")

# modularity scores w/ outer quantiles for different cluster sizes, using optimal step size
# suggests that you don't necessarily get a huge benefit from more than 9 clusters
tune_grid2 %>%
  filter(step == opt.mod2$step) %>%
  ggplot(aes(x = cluster, y = mod)) +
  geom_hline(aes(yintercept = opt.mod2$mod_q025), lty = "dashed") +
  geom_errorbar(aes(ymin = mod_q025, ymax = mod_q975), width = 0) +
  geom_point(shape = 21, size = 2, fill = "white") + 
  theme_bw() + mythemes

ggsave("figures/modularity_quantiles_optimal_steps.png")

# boxplot of modularity for different cluster sizes, using optimal step size

png("figures/modularity_boxplot_optimal_steps.png", 
    width=12, height=8, units='in', res=300)

boxplot(modularity.l2[[opt.mod$step]],
        type="l", xlab= "clusters (k)", 
        ylab= "modularity (Q)", 
        main = paste("Modularity changing with k (Walktrap, steps = ", 
                     opt.mod$step, ")"))

dev.off()

## (4) run full cluster analysis on full dataset of tags using optimal number of steps and clusters
# create undirected adjacency matrix from full tag matrix
fl.graph2 <- graph_from_adjacency_matrix(tags.m,
                                        weighted=TRUE, 
                                        mode="undirected",
                                        diag=TRUE)

# clustering from full matrix
set.seed(666)
cw2 <- cluster_walktrap(fl.graph2, 
                       weights = E(fl.graph2)$weight,
                       membership = T,
                       steps = opt.mod2$step)

# optimal number of clusters (based on 95% quantile range of top modularity)
# nclust <- min(tune_grid$cluster[
  # tune_grid$step == opt.mod$step & tune_grid$mod_q975 >= 
    # tune_grid[which.max(tune_grid$mod),]$mod_q025])

nclust2 <- opt.mod2$cluster

cw.cut2 <- cut_at(cw, no = nclust2)
modularity(fl.graph2, membership = cw.cut2)
table(cw.cut2)

## eigenvector ranking
### Tag importance based on the whole network
ec2 <- eigen_centrality(fl.graph2)$vector
cl2 <- closeness(fl.graph2)
bt2 <- betweenness(fl.graph2)
dg2 <- centr_degree(fl.graph2)

### output table of tags, community identity, and eigenvector centrality (importance)
cw.comm2 <- data.frame(
  tag = cw2$names,
  cluster = c(cw.cut2),
  eigen_centrality = ec2,
  betweenness = bt2, 
  closeness = cl2, 
  degree=dg2
)

summary(dg2$res)

pdf("figures/centrality_eigenvalue.pdf", width = 12, height = 8)
barplot(sort(ec2, decreasing = T)[1:30], las=2)
dev.off()

pdf("figures/centrality_closeness.pdf", width = 12, height = 8)
barplot(sort(cl2, decreasing = T)[1:30], las=2)
dev.off()

pdf("figures/centrality_betweenness.pdf", width = 12, height = 8)
barplot(sort(bt2, decreasing = T)[1:30], las=2)
dev.off()

tag.ranks.omitted <- matrix(nrow = max(table(cw.cut2)), ncol = nclust2)
for (clust in 1:nclust2) {
  cw.eigen.rank <- cw.comm2 %>%
    filter(cluster == clust) %>%
    arrange(desc(eigen_centrality))
  if (
    nrow(cw.eigen.rank) == max(table(cw.cut2))
  )
    tag.ranks.omitted[,clust] <- cw.eigen.rank$tag
  else
    tag.ranks.omitted[,clust] <- 
      c(cw.eigen.rank$tag, 
        rep(NA, length = max(table(cw.cut2))-nrow(cw.eigen.rank)))
}

tag.ranks.omitted
write.csv(tag.ranks.omitted, "data/rural_tag_clusters_toptagsomit.csv", 
          na = "", row.names = FALSE)

# ASSIGN CLUSTERS TO IMAGES BASED ON TAG CONTENT -------------------------------

## aggregate clusters into core groups

## 1:  arts/sports (1,16,18,26,28)*                  [tech, music, sports, photography]
## 2:  scenery (2)**                                 [megacluster, mixture of natural and cultural landscape features]
## 3:  food/dining (3,6,7,8,14,15)                   [food, dessert, dining, alcohol] 
## 4:  aquatic recreation (4,5)                      [watercrafts, fishing, and watersports]
## 5:  natural life (9,11,17,19,22,25)**             [flora, fauna, and fungi (and some livestock and zoo animals)]
## 6:  equestrian (10)                               [horses]
## 7:  people (12)**                                 [people]
## 8:  pets (13)                                     [dogs and cats]
## 9:  transportation (20,24)                        [cars, trucks, buses, trains, bicycles, aircraft]
## 10: firearms and shooting (21,23,29)              [guns, hunting, and shooting ranges]

cw.comm2 <- cw.comm2 %>%
  mutate(theme = case_when(cluster %in% c(1,16,18,26,28) ~ "arts/sports",
                           cluster == 2 ~ "scenery",
                           cluster %in% c(3,6,7,8,14,15) ~ "food/dining",
                           cluster %in% c(4,5) ~ "aquatic recreation",
                           cluster %in% c(9,11,17,19,22,25) ~ "natural life",
                           cluster == 10 ~ "equestrian",
                           cluster == 12 ~ "people",
                           cluster == 13 ~ "dogs/cats",
                           cluster %in% c(20,24) ~ "transportation",
                           cluster %in% c(21,23,29) ~ "firearms/hunting"))


## assign themes to photos based on tags, randomly split ties
assign.theme <- function(photos, clusters) {
  message("Assigning themes...")
  out <- matrix(nrow = nrow(photos), ncol = ncol(photos) + 1,
                dimnames = list(c(1:nrow(photos)),
                                c(colnames(photos), "theme")))
  for (i in 1:nrow(photos)) {
    tag.seq <- c(photos[i, c("tags1", "tags2", "tags3", "tags4", "tags5",  
                             "tags6", "tags7", "tags8","tags9", "tags10")])
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
    out[i, "theme"] <- photo.theme
  } 
  return(out)
  beep(3)
  message("done!")
}

rural_photos_clustered2 <- assign.theme(photos = flickr_rural_tidy_tagged, 
                                       clusters = cw.comm2)

# quick bar plot of cluster composition
# quick bar plot of cluster composition
rural.themes.plot2 <-
  ggplot(as.data.frame(sort(table(rural_photos_clustered2[,"theme"]))),
         aes(x = Var1, y = Freq)) +
  geom_col() + 
  geom_text(aes(label = round(Freq/sum(Freq), 3), y = Freq + 1000, angle = 90)) +
  geom_text(aes(label = paste("\ \ \ \ \ n = ", Freq), y = Freq + 5000), col = "grey50") +
  labs(x = "theme", y = "n (images)") + 
  theme_classic() + coord_flip() +
  labs(title = "Rural image themes (common tags [>10% of photos] removed)")

rural.themes.plot2

library(ggpubr)
ggarrange(rural.themes.plot, rural.themes.plot2, nrow = 2)


ggsave("figures/photo_themes_contr.png", 
       width = 8, height = 8, units = "in", dpi = 600)

beep(3)




