# Kathleen Noller (C) 2024
# katkats1@jh.edu

#' Compute stemFinder scores for a query scRNA-seq dataset
#' 
#' 
#' 
#' @param adata Seurat object containing normalized, log-transformed, and scaled gene expression data
#' @param nn matrix of K-nearest neighbor assignments
#' @param k number of neighbors identified per single cell
#' @param thresh threshold value for binarizing scaled gene expression matrix 
#' @param markers character vector with names of marker genes (cell cycle genes)
#' @param method string denoting method of computing cell cycle gene expression heterogeneity ('gini','stdev' or 'variance')
#' 
#' @return Seurat object with two new metadata columns: 
#' stemFinder: single-cell stemFinder score, where lower values correspond to relatively less differentiated cells within the query dataset and vice versa (corresponding to pseudotime)
#' stemFinder_raw: raw score before the inversion step, where lower values correspond to relatively more differentiated cells in the query dataset and vice versa
#' 
run_stemFinder <-
function(adata, nn = knn, k = k, thresh = 0, markers = markers, method = 'gini'){
  
  if(method == 'gini'){
    
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
    adata$stemFinder_raw = gini_agg$gini_index_agg 
    
  }else if(method == 'stdev'){
    
    sd_agg = data.frame("cell" = colnames(nn), "sd_index_agg" = rep(NA, ncol(nn)))
    rownames(sd_agg) = colnames(nn)
    expDat = as.matrix(adata@assays$RNA@data)[markers,] # use log-transformed, normalized counts
    
    for (cell in colnames(nn)){
      neigh = names(nn[cell,][nn[cell,]== 1])
      sd = apply(expDat[,neigh], 1, sd)
      sd_agg[cell,]$sd_index_agg = sum(sd)
      }
    adata$stemFinder_raw = sd_agg$sd_index_agg
    
  }else if(method == 'variance'){
    
    var_agg = data.frame("cell" = colnames(nn), "var_index_agg" = rep(NA, ncol(nn)))
    rownames(var_agg) = colnames(nn)
    expDat = as.matrix(adata@assays$RNA@data)[markers,] # use log-transformed, normalized counts
    
    for (cell in colnames(nn)){
      neigh = names(nn[cell,][nn[cell,]== 1])
      var = apply(expDat[,neigh], 1, var)
      var_agg[cell,]$var_index_agg = sum(var)
    }
    adata$stemFinder_raw = var_agg$var_index_agg
    
  }
  
  adata$stemFinder = 1 - (adata$stemFinder_raw)/max(adata$stemFinder_raw) #inverted to correspond to pseudotime

  return(adata)
}