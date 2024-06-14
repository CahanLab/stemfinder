# Kathleen Noller
# katkats1@jh.edu

#' Identify transcription factors with high Gini index in one phenotype-defined cluster and low Gini index in another phenotype-defined cluster
#' Selects transcription factors that are differentially expressed in 1+ phenotype-defined clusters
#' 
#' @param adata Seurat object containing scaled gene expression matrix and metadata column with cell type annotations ("Phenotype") 
#' @param lowthresh_tf lower value for quantile-based thresholding of TFs by Gini index
#' @param highthresh_tf upper value for quantile-based thresholding of TFs by Gini index
#' @param min.pct value for differential expression analysis in Seurat; minimum of fraction of cells expressing a given gene in a given cluster
#' @param logfc.threshold value for DE analysis in Seurat; log fold change threshold
#' 
#' @return character vector of gene names meeting selected criteria
#' 
#' 
find_key_tfs <- function(adata, TFs, lowthresh_tf = 0.4, highthresh_tf = 0.6, min.pct = 0.3, logfc.threshold = 0.5){

  markers = TFs[TFs %in% rownames(adata)]
  expDat = as.matrix(adata@assays$RNA@scale.data)
  sampTab = adata@meta.data
  idents = unique(adata$Phenotype)
  idents = names(which(table(adata$Phenotype) != 1)) 

  gini_agg = data.frame("cluster" = rep(idents, length(markers)), "gini_sum" = rep(NA, length(idents)*length(markers)), "TF" = rep("", length(idents)*length(markers)))
  
  for (i in idents){
    gini_agg[gini_agg$cluster ==i,]$TF = markers
    
    neigh = rownames(sampTab[sampTab$Phenotype == i,]) 
    exp = expDat[markers,neigh] > 0 
    
    gini_allsingles = c()
    
    for (c in neigh){
      gini_singles = data.frame("cell" = rep(c, length(markers)), "gini" = rep(NA, length(markers)), "TF" = markers)
      
      exp_match = exp == exp[,c] 
      n_match = apply(exp_match, 1, sum) - 1
      p_g = n_match/(length(neigh)-1) 
      gini_g = p_g * (1 - p_g) 
      
      gini_singles$gini = gini_g
      gini_allsingles = rbind(gini_allsingles, gini_singles)
    }
    
    for (m in markers){
      gini_agg[gini_agg$cluster == i & gini_agg$TF == m,]$gini_sum = sum(gini_allsingles[gini_allsingles$TF == m,]$gini)
    }
  }
  
  finals_down = gini_agg[gini_agg$gini_sum < quantile(gini_agg$gini_sum, lowthresh_tf),]$TF
  finals_up = gini_agg[gini_agg$gini_sum > quantile(gini_agg$gini_sum, highthresh_tf),]$TF
  finals_stoch = unique(finals_down[finals_down %in% finals_up])
  if(length(finals_stoch)==0){
    stop(print("No sporadically expressed TFs found with given thresholds. Please select more permissive thresholds."))
  }
  
  Idents(adata) = 'Phenotype'
  markers_de = FindAllMarkers(adata, min.pct = min.pct, logfc.threshold = logfc.threshold, only.pos=T, test.use="wilcox")
  markers_de = markers_de[markers_de$p_val_adj < 0.05,]$gene
  finals = finals_stoch[finals_stoch %in% markers_de]
  print("Sporadically expressed TFs overlapping with DE genes:")
  print(finals)
  if(length(finals) == 0 & length(finals_stoch > 0)){
    stop(print("No overlapping TFs found with given thresholds for differential gene expression. Please select more permissive thresholds."))
  }
  
  return(finals)
}
