compute_performance <-
function(adata, d = 1, id, list_all = c()){
  
  #Single-cell Spearman correlation
  spear_all_sf = cor.test(x = adata$stemFinder_invert, y = adata$Ground_truth, method = "spearman", exact = F)$estimate
  spear_all_ccat = cor.test(x = adata$ccat_invert, y = adata$Ground_truth, method = "spearman", exact = F)$estimate
  spear_all_cyto = cor.test(x = adata$CytoTRACE_invert, y = adata$Ground_truth, method = "spearman", exact = F)$estimate
  print(paste("Single-cell Spearman Correlation, stemFinder:", spear_all_sf, sep = " "))
  print(paste("Single-cell Spearman Correlation, CCAT:", spear_all_ccat, sep = " "))
  print(paste("Single-cell Spearman Correlation, CytoTRACE:", spear_all_cyto, sep = " "))
  
  #AUC
  mostpotent = adata@meta.data[adata$Ground_truth == max(adata$Ground_truth),]
  leastpotent = adata@meta.data[adata$Ground_truth == min(adata$Ground_truth),]
  category = c(rep(1, nrow(mostpotent)), rep(0, nrow(leastpotent))) 
  prediction.sf = c(mostpotent$stemFinder_invert, leastpotent$stemFinder_invert)
  prediction.cyto = c(mostpotent$CytoTRACE_invert, leastpotent$CytoTRACE_invert)
  prediction.ccat = c(mostpotent$ccat_invert, leastpotent$ccat_invert)#calculator results
  auc_sf = auc_probability(category, prediction.sf)
  auc_ccat = auc_probability(category, prediction.ccat)
  auc_cyto = auc_probability(category, prediction.cyto)
  
  #Phenotypic Spearman correlation 
  clusters = as.character(unique(adata$Phenotype))
  meanpotency = data.frame("cluster" = clusters, "sf_potency" = rep(NA, length(clusters)), "ccat_potency" = rep(NA, length(clusters)), "cyto_potency" = rep(NA, length(clusters)), "ground_truth" = rep(NA, length(clusters)))
  for (i in clusters){
    meanpotency$sf_potency[which(clusters == i)] = mean(adata$stemFinder_invert[adata$Phenotype == i])
    meanpotency$cyto_potency[which(clusters == i)] = mean(adata$CytoTRACE_invert[adata$Phenotype == i])
    meanpotency$ccat_potency[which(clusters == i)] = mean(adata$ccat_invert[adata$Phenotype == i])
    meanpotency$ground_truth[which(clusters == i)] = mean(adata$Ground_truth[adata$Phenotype == i])
  }
  spear_pheno_sf = cor.test(x = meanpotency$sf_potency, y = meanpotency$ground_truth)$estimate
  spear_pheno_ccat = cor.test(x = meanpotency$ccat_potency, y = meanpotency$ground_truth)$estimate
  spear_pheno_cyto = cor.test(x = meanpotency$cyto_potency, y = meanpotency$ground_truth)$estimate
  
  #Store in list
  list_all$dataset[d] = id
  list_all$spear_all_sf[d] = spear_all_sf
  list_all$spear_all_ccat[d] = spear_all_ccat
  list_all$spear_all_cyto[d] = spear_all_cyto

  list_all$auc_sf[d] = auc_sf
  list_all$auc_ccat[d] = auc_ccat
  list_all$auc_cyto[d] = auc_cyto
  
  list_all$spear_pheno_sf[d] = spear_pheno_sf
  list_all$spear_pheno_ccat[d] = spear_pheno_ccat
  list_all$spear_pheno_cyto[d] = spear_pheno_cyto
  
  return(list_all)
}
