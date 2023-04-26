\name{find_key_tfs}
\alias{find_key_tfs}
\title{
Identify key TFs with putative roles in lineage development}
\description{
Compute Gini index of each TF in reference TF list for a given species.
Return TFs with high Gini score in one cluster, low score in another cluster, which overlap with DE genes. 
}
\usage{
finals = find_key_tfs(adata, TFs, lowthresh_tf = 0.4, highthresh_tf = 0.6, min.pct = 0.3, logfc.threshold = 0.5)
}
\arguments{
  \item{adata}{
Seurat object containing scaled gene expression matrix }
  \item{TFs}{
reference list of TFs for proper species }
  \item{lowthresh_tf}{
Quantile-based threshold to select TFs with low stochasticity in one cluster}
  \item{highthresh_tf}{
Quantile-based threshold to select TFs with high stochasticity in one cluster}
  \item{min.pct}{
For differential expression in Seurat; threshold for minimum percent cells in given cluster expressing a given DE gene}
  \item{logfc.threshold}{
For differential expression in Seurat; threshold for log fold change expression}
}
\value{
List of TF names 
}
\author{
K. Noller }