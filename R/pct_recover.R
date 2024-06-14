# Kathleen Noller (C) 2024
# katkats1@jh.edu

#' Compute percentage of cells with low ground truth degree of differentiation which are correctly given a low inverted stemFinder scores
#' Requires knowledge of ground truth degree of differentiation (stored in metadata as 'Ground_truth' numeric vector)
#' 
#' Also computes relative abundance of cells with lowest ground truth degree of differentiation for a given dataset
#' 
#' @param adata Seurat object containing 'Ground_truth' and 'stemFinder_invert' numeric metadata columns  
#' 
#' @return numeric value of percentage of cells with low ground truth degree of differentiation identified by stemFinder
pct_recover <-
function(adata){
  
  num.stem.gt = nrow(adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),])
  pct.recov = (nrow(adata@meta.data[adata$stemFinder_invert < quantile(adata$stemFinder_invert, (1/length(unique(adata$Ground_truth)))) & adata$Ground_truth == min(adata$Ground_truth),])/num.stem.gt) * 100
  print(paste("Percentage of cells with low degree of differentiation recovered by stemFinder:", pct.recov, sep = " "))
  
  rel.abund = (nrow(adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),]) / nrow(adata@meta.data)) * 100
  print(paste("Relative abundance of cells with low degree of differentiation:", rel.abund, sep = " "))

  return(pct.recov)
}
