\name{gene_set_score}
\alias{gene_set_score}
\title{Mean expression of set of genes}
\description{
Computes mean expression of provided gene set}
\usage{
adata = gene_set_score(gene.set, adata)
}
\arguments{
  \item{adata}{
Seurat object}
  \item{gene.set}{
Character vector of gene names of interest}
}
\value{
Returns adata with metadata column gene.set.score
}
\author{
K. Noller }
