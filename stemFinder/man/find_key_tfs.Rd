\name{find_key_tfs}
\alias{find_key_tfs}
\title{
Find differentially expressed TFs with high stochasticity in one cluster and low stochasticity in a separate cluster}
\description{
Identify key TFs with putative roles in lineage development}
\usage{
find_key_tfs(adata, TFs = mmTFs, stoch_thresh_low = 0.3, stoch_thresh_high = 0.7, min.pct = 0.4, log.fc = 0.75)
}
\arguments{
  \item{adata}{
Seurat object containing scaled expression matrix }
  \item{TFs}{
reference list of TFs }
  \item{stoch_thresh_low}{
Quantile threshold to select TFs with low stochasticity in one cluster}
  \item{stoch_thresh_high}{
Quantile threshold to select TFs with high stochasticity in one cluster}
  \item{min.pct}{
For differential expression in Seurat; threshold for minimum percent cells in given cluster expressing a given DE gene}
  \item{log.fc}{
For differential expression in Seurat; threshold for log fold change}
}
\value{
List of TF names + equal length cell cycle genes for later input into run_stemFinder
}
\author{
K. Noller }