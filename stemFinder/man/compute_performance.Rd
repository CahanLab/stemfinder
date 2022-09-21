\name{compute_performance}
\alias{compute_performance}
\title{Quantify performance of potency calculators}
\description{
Compute single-cell Spearman correlation, phenotypic Spearman correlation, and AUC for potency calculator results}
\usage{
compute_performance(adata, d = 1, id, list_all = c())
}
\arguments{
  \item{adata}{
Seurat object containing results of 3 potency calculators in metadata columns: CytoTRACE_invert, ccat_invert, and stemFinder_invert}
  \item{d}{
Index value for list entry}
  \item{id}{
String describing dataset}
  \item{list_all}{
empty list, or list containing performance information from prior runs (if d > 1)}
}

\value{
a list containing performance metrics for 3 potency calculator methods
}

\author{
K. Noller
}