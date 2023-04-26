
gene_set_score <- function(gene.set, adata){
  
  gene.set = gene.set[gene.set %in% rownames(adata)]
  mean.exp <- colMeans(x = adata@assays$RNA@data[gene.set, ], na.rm = TRUE)
  adata@meta.data$gene.set.score <- mean.exp
  
  return(adata)
}
