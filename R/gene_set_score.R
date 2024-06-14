# Kathleen Noller (C) 2024
# katkats1@jh.edu

#' Obtain mean expression of a provided gene list in a scRNA-seq query dataset
#' 
#' @param gene.set character vector of gene names 
#' @param adata Seurat object containing normalized, log-transformed gene expression data 
#' 
#' @return Seurat object with gene set score added as metadata column ("gene.set.score")
gene_set_score <- function(gene.set, adata){
  
  gene.set = gene.set[gene.set %in% rownames(adata)]
  mean.exp <- colMeans(x = adata@assays$RNA@data[gene.set, ], na.rm = TRUE)
  adata@meta.data$gene.set.score <- mean.exp
  
  return(adata)
}
