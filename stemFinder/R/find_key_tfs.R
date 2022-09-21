find_key_tfs <-
function(adata, TFs = mmTFs, stoch_thresh_low = 0.3, stoch_thresh_high = 0.7, min.pct = 0.4, log.fc = 0.75){

  markers = TFs[TFs %in% rownames(adata)]
  expDat = as.matrix(adata@assays$RNA@scale.data)
  sampTab = adata@meta.data
  idents = unique(adata$Phenotype)

  gini_agg = data.frame("cluster" = rep(idents, length(markers)), "gini_sum" = rep(NA, length(idents)*length(markers)), "TF" = rep("", length(idents)*length(markers)))

  #Compute stochasticity for each TF in each phenotype-defined cluster
  
  for (i in idents){
    gini_agg[gini_agg$cluster ==i,]$TF = markers
    neigh = rownames(sampTab[sampTab$Phenotype == i,]) 
    exp = expDat[markers,neigh] > 0 
    
    gini_allsingles = c()
    for (c in neigh){
      gini_singles = data.frame("cell" = rep(c, length(markers)), "gini" = rep(NA, length(markers)), "TF" = markers)
      
      exp_match = exp == exp[,c] 
      n_match = apply(exp_match, 1, sum) - 1 
      p_g = n_match/length(neigh) 
      gini_g = p_g * (1 - p_g) 
      
      gini_singles$gini = gini_g 
      gini_allsingles = rbind(gini_allsingles, gini_singles)
    }
    
    for (m in markers){
      gini_agg[gini_agg$cluster == i & gini_agg$TF == m,]$gini_sum = sum(gini_allsingles[gini_allsingles$TF == m,]$gini)
    }
  }

  #Select TFs with high stochasticity in one cluster and low stochasticity in another
  finals_down = gini_agg[gini_agg$gini_sum < quantile(gini_agg$gini_sum, stoch_thresh_low),]$TF
  finals_up = gini_agg[gini_agg$gini_sum > quantile(gini_agg$gini_sum, stoch_thresh_high),]$TF
  finals_stoch = unique(finals_down[finals_down %in% finals_up])
  
  #Narrow down to differentially expressed TFs
  Idents(adata) = 'Phenotype'
  markers = FindAllMarkers(adata, min.pct = min.pct, logfc.threshold = log.fc, only.pos=T, test.use="MAST")
  markers = markers[markers$p_val_adj < 0.05,]$gene
  finals = finals_stoch[finals_stoch %in% markers]
  print(paste("Selected TFs:", finals, sep = " "))
  if (length(finals) == 0){
    print("No TFs selected; select less stringent thresholds.")
  }
  
  #Return list of selected TFs + equal length cell cycle genes
  cc_random = sample(cell_cycle_genes, size = length(finals), replace = F)
  markers = c(finals, cc_random)
  
  return(markers)
  
}
