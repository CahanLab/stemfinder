\name{run_stemFinder}
\alias{run_stemFinder}
\title{
Compute single-cell differentiation time with stemFinder}
\description{
Compute single-cell differentiation time from scRNA-seq data by computing heterogeneity of cell cycle marker gene expression}
\usage{
adata = run_stemFinder(adata, nn = knn, k = k, thresh = 0, markers = markers, method = 'gini')
}
\arguments{
  \item{adata}{
Seurat object containing log-normalized, scaled gene expression matrix}
  \item{nn}{
K nearest neighbors matrix, computed using FindNearestNeighbors() in Seurat}
  \item{k}{
Numeric input to K nearest neighbors}
  \item{thresh}{
Numeric threshold for binarizing gene expression matrix}
  \item{markers}{
Character vector of gene names found in query data (default: S and G2M genes from Scanpy gene set scoring function)}
  \item{method}{
String denoting method of quantifying gene expression heterogeneity (default: 'gini', other: 'stdev' or 'variance')}
  }
}

\value{
Returns adata containing 2 columns:
1. stemFinder_raw: raw score
2. stemFinder: inverted score, comparable to pseudotime (lowest score = lowest degree of differentiation or highest fate potential)
}

\author{
K. Noller }