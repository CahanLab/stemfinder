# Kathleen Noller (C) 2024
# katkats1@jh.edu

#' Computes percent highly potent cells (lowest Ground_truth value) with stemFinder scores below a quantile-based threshold
#' Also computes relative abundance of cells with lowest ground truth differentiation time for a given dataset
#' 
#' Requires knowledge of ground truth differentiation time (stored in metadata as 'Ground_truth' numeric vector)
#' 
#' @param adata Seurat object containing 'Ground_truth' and 'stemFinder' numeric metadata columns  
#' 
#' @return numeric value of percentage of cells with low ground truth differentiation time identified by stemFinder
pct_recover <-
function(adata){
  
  num.stem.gt = nrow(adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),])
  pct.recov = (nrow(adata@meta.data[adata$stemFinder < quantile(adata$stemFinder, (1/length(unique(adata$Ground_truth)))) & adata$Ground_truth == min(adata$Ground_truth),])/num.stem.gt) * 100
  print(paste("Percentage of cells with low degree of differentiation recovered by stemFinder:", pct.recov, sep = " "))
  
  rel.abund = (nrow(adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),]) / nrow(adata@meta.data)) * 100
  print(paste("Relative abundance of cells with low degree of differentiation:", rel.abund, sep = " "))

  return(pct.recov)
}
