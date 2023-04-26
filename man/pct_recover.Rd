\name{pct_recover}
\alias{pct_recover}
\title{Percent highly potent cells recovered by stemFinder }
\description{
Computes percent highly potent cells (lowest Ground_truth value) with stemFinder_invert scores below a quantile-based threshold.}
\usage{
pct.recov = pct_recover(adata)
}
\arguments{
  \item{adata}{
Seurat object containing metadata column with ground truth potency ('Ground_truth', numeric) and metadata column with inverted stemFinder scores ('stemFinder_invert') from run_stemFinder function}
}
\value{
pct.recov 
}
\author{
K. Noller }
