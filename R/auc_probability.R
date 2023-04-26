auc_probability <-
function(labels, scores, N=1e7){
  pos <- sample(scores[labels], N, replace=TRUE)
  neg <- sample(scores[!labels], N, replace=TRUE)
  (sum(pos > neg) + sum(pos == neg)/2) / N # give partial credit for ties
}
