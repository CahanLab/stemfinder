\name{process_cyto_adata}
\alias{process_cyto_adata}
\title{Process datasets from CytoTRACE website}
\description{Code to process datasets from CytoTRACE websites (as noted in Table 1)}
\usage{
adata = process_cyto_adata(adata, mitotag = "^mt-", ribotag = "^Rp[sl]", mtmax = mtmax, s_genes = s_genes, g2m_genes = g2m_genes)
}
\arguments{
  \item{adata}{
Downloaded .rds }
  \item{mitotag}{
string to find mitochondrial genes}
  \item{ribotag}{
string to find ribosomal genes}
  \item{mtmax}{
max percent.mt for filtering}
  \item{s_genes}{
S phase cell cycle genes in the query data expression matrix}
  \item{g2m_genes}{
G2M phase cell cycle genes in the query data expression matrix}
}
\details{
https://cytotrace.stanford.edu/
}
\references{
Gulati, G. S., Sikandar, S. S., Wesche, D. J., Manjunath, A., Bharadwaj, A., Berger, M. J., Ilagan, F., Kuo, A. H., Hsieh, R. W., Cai, S., Zabala, M., Scheeren, F. A., Lobo, N. A., Qian, D., Yu, F. B., Dirbas, F. M., Clarke, M. F., & Newman, A. M. (2020). Single-cell transcriptional diversity is a hallmark of developmental potential. Science (New York, N.Y.), 367(6476), 405–411. https://doi.org/10.1126/science.aax0249}
\author{
K. Noller }