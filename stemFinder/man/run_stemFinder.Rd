\name{run_stemFinder}
\alias{run_stemFinder}
\title{
Compute single-cell potency with stemFinder}
\description{
Compute single-cell potency from scRNA-seq data by computing stochasticity of cell cycle marker gene expression}
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
Returns adata containing 3 columns:
1. stemFinder: raw potency score
2. stemFinder_invert: inverted potency score, comparable to pseudotime (lowest score = most potent cell)
3. stemFinder_comp: potency score that is comparable across different datasets 
}

\author{
K. Noller }