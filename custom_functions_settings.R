# R CUSTOM FUNCTIONS / SETTINGS FOR ANALYSES

library(dplyr)
library(beepr)

## is.nan function for data frames
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

## notin
`%notin%` <- Negate(`%in%`)

## mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## percentile rank
p_rank <- function(x) ifelse(is.na(x),NA,rank(x)/length(x[is.na(x) == FALSE]))

# function to get closest quantile value
percentile_level <- function(x){
  probs_vec <- quantile(x, probs = seq(0,1,0.001))
  out <- data.frame(orig = x, perc = 0)
  for(i in 1:length(x)){
    abs.diff = abs(x[i] - probs_vec)
    prob.ind = which.min(abs.diff)
    perc.val = as.numeric(str_replace_all(names(probs_vec[prob.ind]), "%", ""))/100
    out[i,2] = perc.val
  }
  return(out[,2])
}

## floor to any interval
floor_any = function(x, accuracy, f=floor){f(x/ accuracy) * accuracy}

# VIF from Zuur (2009)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

# Support function for corvif. Will not be called by the user.
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

# GRAPHICS
## color blind friendly palette (with grey):
col1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## color_blind_friendly palette (with black): 
col2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## larger text size for ggplot
library(ggplot2)
mythemes <- list(
  theme(strip.background = element_rect(color = NA),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16)))