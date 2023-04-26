pct_recover <-
function(adata){
  
  num.stem.gt = nrow(adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),])
  pct.recov = (nrow(adata@meta.data[adata$stemFinder_invert < quantile(adata$stemFinder_invert, (1/length(unique(adata$Ground_truth)))) & adata$Ground_truth == min(adata$Ground_truth),])/num.stem.gt) * 100
  print(paste("Percentage highly potent cells recovered by stemFinder:", pct.recov, sep = " "))
  
  rel.abund = (nrow(adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),]) / nrow(adata@meta.data)) * 100
  print(paste("Relative abundance of highly potent cells:", rel.abund, sep = " "))

  return(pct.recov)
}
