\name{compute_performance_single}
\alias{compute_performance_single}
\title{Quantify performance of stemFinder and optionally one competitor method}
\description{
Compute single-cell Spearman correlation, phenotypic Spearman correlation, and AUC for stemFinder results on a query dataset. Optionally compute same performance metrics on a competitor method's results.}
\usage{
list_all = compute_performance_single(adata, competitor = F)
}
\arguments{
  \item{adata}{
Seurat object containing results of stemFinder metadata columns, ground truth differentiation degree, and optionally the results of a competitor method: stemFinder, Ground_truth, and (optional) competitor}
  \item{competitor}{
Logical indicating whether the user wants to compute performance of a competitor method of potency calculation }
  \item{comp.inverted}{
Logical indicating if competitor scores are inverted to correspond to pseudotime}
}

\value{
a list (list_all) containing 3 performance metrics for stemFinder on query data, and optionally performance metrics for a competitor method on query data.
}

\author{
K. Noller
}