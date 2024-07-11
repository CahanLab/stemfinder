# Kathleen Noller (C) 2024
# katkats1@jh.edu

#' Quantify stemFinder performance on a single scRNA-seq dataset compared with ground truth differentiation time
#' Optional: quantify performance of a competitor tool versus ground truth
#' 
#' Note that competitor method scores must be arranged so that the lowest score = the cell with highest "potency", however defined by competitor
#'
#' Comparisons are made using: single-cell Spearman correlation, AUC, and phenotypic Spearman correlation
#' Single-cell Spearman correlation: correlation between single-cell values of ground truth vs. inverted score
#' Phenotypic Spearman correlation: correlation between mean values of ground truth vs. inverted score per phenotype-defined cluster
#' 
#' @param adata Seurat object containing metadata with columns "Ground_truth" (numeric), "stemFinder_invert" (numeric), and "competitor_invert" if competitor = T
#' @param competitor logical denoting whether comparison between competitor method and ground truth is desired
#' 
#' 
#' @return list containing values from quantitative comparison for the query dataset
#' 
compute_performance_single <- function(adata, competitor = F){
  
  #Single-cell Spearman correlation
  spear_all_sf = cor.test(x = adata$stemFinder_invert, y = adata$Ground_truth, method = "spearman", exact = F)$estimate
  print(paste("Single-cell Spearman Correlation, stemFinder:", round(spear_all_sf, 2), sep = " "))
  
  #AUC
  mostpotent = adata@meta.data[adata$Ground_truth == max(adata$Ground_truth),]
  leastpotent = adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),]
  category = c(rep(1, nrow(mostpotent)), rep(0, nrow(leastpotent))) 
  prediction.sf = c(mostpotent$stemFinder_invert, leastpotent$stemFinder_invert)
  auc_sf = auc_probability(category, prediction.sf)
  print(paste("AUC, stemFinder:", round(auc_sf, 2)), sep = " ")
  
  #Phenotypic Spearman correlation 
  clusters = as.character(unique(adata$Phenotype))
  meanpotency = data.frame("cluster" = clusters, "sf_potency" = rep(NA, length(clusters)), "ground_truth" = rep(NA, length(clusters)))
  for (i in clusters){
    meanpotency$sf_potency[which(clusters == i)] = mean(adata$stemFinder_invert[adata$Phenotype == i])
    meanpotency$ground_truth[which(clusters == i)] = mean(adata$Ground_truth[adata$Phenotype == i])
  }
  spear_pheno_sf = cor.test(x = meanpotency$sf_potency, y = meanpotency$ground_truth)$estimate
  print(paste("Phenotypic Spearman correlation, stemFinder:", round(spear_pheno_sf, 2)), sep = " ")
  
  sf_performance = c(spear_all_sf, spear_pheno_sf, auc_sf)
  names(sf_performance) = c('Spearman_SingleCell','Spearman_Pheno','AUC')

  #Competitor performance
  
  if (competitor == T){ 
    comp_scores = adata[['competitor_invert']]
 
    spear_all_comp = cor.test(x = comp_scores[[1]], y = adata$Ground_truth, method = "spearman", exact = F)$estimate
    prediction.comp = c(comp_scores[rownames(mostpotent),], comp_scores[rownames(leastpotent),])
    auc_comp = auc_probability(category, prediction.comp)
    
    meanpotency = data.frame("cluster" = clusters, "comp_potency" = rep(NA, length(clusters)), "ground_truth" = rep(NA, length(clusters)))
    for (i in clusters){
      meanpotency$comp_potency[which(clusters == i)] = mean(comp_scores[rownames(adata@meta.data[adata$Phenotype == i,]),])
      meanpotency$ground_truth[which(clusters == i)] = mean(adata$Ground_truth[adata$Phenotype == i])
    }
    spear_pheno_comp = cor.test(x = meanpotency$comp_potency, y = meanpotency$ground_truth)$estimate
    
    comp_performance = c(spear_all_comp, spear_pheno_comp, auc_comp)
    names(comp_performance) = c('Spearman_SingleCell','Spearman_Pheno','AUC')
  
  }else {comp_performance = NA} 
  
  list_all = list("stemFinder results" = sf_performance, "Competitor results" = comp_performance)
  return(list_all)
}
