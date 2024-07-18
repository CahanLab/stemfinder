# Source: https://www.r-bloggers.com/2016/11/calculating-auc-the-area-under-a-roc-curve/#google_vignette 

#' Approximate AUC as probability that positive case outranks negative
#' Used to quantify the discrimination accuracy for distinguishing the most vs. least differentiated cells in the dataset
#'
#' @param labels numeric vector of ground truth differentiation state of most vs. least differentiated cells in dataset 
#' @param scores numeric vector of stemFinder (or other tool) scores for the cells in 'labels' 
#' 
#' @return area under the curve (AUC) ranging from 0 to 1
#' 
auc_probability <-
function(labels, scores, N=1e7){
  pos <- sample(scores[labels], N, replace=TRUE)
  neg <- sample(scores[!labels], N, replace=TRUE)
  (sum(pos > neg) + sum(pos == neg)/2) / N # give partial credit for ties
}