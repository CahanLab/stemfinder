\name{run_stemFinder}
\alias{run_stemFinder}
\title{
Compute single-cell differentiation time with stemFinder}
\description{
Compute single-cell differentiation time from scRNA-seq data by computing heterogeneity of cell cycle marker gene expression}
\usage{
adata = run_stemFinder(adata, nn = knn, k = k, thresh = 0, markers = markers)
}
\arguments{
  \item{adata}{
Seurat object containing scaled gene expression matrix}
  \item{nn}{
K nearest neighbors matrix, computed using FindNearestNeighbors() in Seurat}
  \item{k}{
Numeric input to K nearest neighbors}
  \item{thresh}{
Numeric threshold for binarizing gene expression matrix}
  \item{markers}{
Character vector of gene names (default: S and G2M genes from Scanpy gene set scoring function)}
}

\value{
Returns adata containing 2 columns:
1. stemFinder: raw score
2. stemFinder_invert: inverted score, comparable to pseudotime (lowest score = lowest degree of differentiation or highest fate potential)
}

\author{
K. Noller }