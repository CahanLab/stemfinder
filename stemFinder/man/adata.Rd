\name{adata}
\alias{adata}
\docType{data}
\title{
stemFinder example dataset
}
\description{
Adult mouse Bone marrow, 10X, GSE109774
}
\usage{data("adata")}
\format{
Seurat package
Assay: RNA
Metadata columns: 
orig.ident, nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, S.Score, G2M.Score, Phase
Phenotype, Ground_truth (numeric assignment, lowest value = highest potency)
Potency scores: CytoTRACE and CytoTRACE_invert, ccat and ccat_invert, stemFinder and stemFinder_invert, stemFinder_comp
\source{
GSE109774
}
\references{
https://pubmed.ncbi.nlm.nih.gov/30283141/
}
