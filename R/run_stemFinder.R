run_stemFinder <-
function(adata, nn = knn, k = k, thresh = 0, markers = markers){
  
  gini_agg = data.frame("cell" = colnames(nn), "gini_index_agg" = rep(NA, ncol(nn)))
  rownames(gini_agg) = colnames(nn)
  
  expDat = as.matrix(adata@assays$RNA@scale.data)[markers,]
  
  for (cell in colnames(nn)){
    neigh = names(nn[cell,][nn[cell,]== 1]) 
    exp = expDat[markers,neigh] > thresh 
    exp_match = exp == exp[,cell] 
    n_match = apply(exp_match, 1, sum) - 1  
    p_g = n_match/(k - 1)
    gini_g = p_g * (1 - p_g) 
    
    gini_agg[cell,]$gini_index_agg = sum(gini_g) 
  }
  
  adata$stemFinder = gini_agg$gini_index_agg 
  adata$stemFinder_invert = 1 - (adata$stemFinder)/max(adata$stemFinder)
  adata$stemFinder_comp = adata$stemFinder/length(markers) #can compare this score across datasets
  
  return(adata)
}