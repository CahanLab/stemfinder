run_stemFinder <-
function(adata, expDat, nn = knn, k = k, thresh = 0, markers = markers){
  
  gini_agg = data.frame("cell" = colnames(nn), "gini_index_agg" = rep(NA, ncol(nn)))
  
  for (cell in rownames(nn)){
    neigh = names(nn[cell,][nn[cell,]== 1]) 
    exp = expDat[markers,neigh] > thresh 
    exp_match = exp == exp[,cell] 
    n_match = apply(exp_match, 1, sum) - 1  
    p_g = n_match/k 
    gini_g = p_g * (1 - p_g) 
    
    gini_agg[gini_agg$cell == cell,]$gini_index_agg = sum(gini_g) 
  }
  
  adata$stemFinder = gini_agg$gini_index_agg 
  adata$stemFinder_invert = 1 - (adata$stemFinder)/max(adata$stemFinder)
  adata$stemFinder_comp = adata$stemFinder/length(markers) #can compare this score across datasets
  
  return(adata)
}
