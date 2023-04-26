\name{compute_performance_single}
\alias{compute_performance_single}
\title{Quantify performance of stemFinder and optionally one competitor method}
\description{
Compute single-cell Spearman correlation, phenotypic Spearman correlation, and AUC for stemFinder results on a query dataset. Optionally compute same performance metrics on a competitor method's results (CytoTRACE or CCAT).}
\usage{
list_all = compute_performance_single(adata, competitor = T, comp_id = "CytoTRACE")
}
\arguments{
  \item{adata}{
Seurat object containing results of stemFinder metadata columns and optionally the results of a competitor method: CytoTRACE_invert or ccat_invert, and stemFinder_invert}
  \item{competitor}{
Logical indicating whether the user wants to compute performance of a competitor method of potency calculation }
  \item{comp_id}{
String describing competitor method: "CytoTRACE" or "CCAT"}
}

\value{
a list (list_all) containing 3 performance metrics for stemFinder on query data, and optionally performance metrics for a competitor method on query data.
}

\author{
K. Noller
}