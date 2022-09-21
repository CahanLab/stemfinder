\name{run_stemFinder}
\alias{run_stemFinder}
\title{
Compute single-cell potency with stemFinder}
\description{
Compute single-cell potency from scRNA-seq data by computing stochasticity of cell cycle marker gene expression}
\usage{
run_stemFinder(adata, expDat, nn = knn, k = k, thresh = 0, markers = markers)
}
\arguments{
  \item{adata}{
Seurat object}
  \item{expDat}{
matrix of scaled gene expression}
  \item{nn}{
K nearest neighbors matrix}
  \item{k}{
input to K nearest neighbors}
  \item{thresh}{
threshold for binarizing gene expression matrix}
  \item{markers}{
cell cycle genes}
}

\value{
Returns adata containing 3 columns:
1. stemFinder: raw potency score
2. stemFinder_invert: inverted potency score, comparable to pseudotime (lowest score = most potent cell)
3. stemFinder_comp: potency score that is comparable across different datasets 
}

\author{
K. Noller }