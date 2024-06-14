# Kathleen Noller (C) 2024
# katkats1@jh.edu

#' Convert single-cell transcriptomic data from CytoTRACE website into format needed for running stemFinder
#' Performs basic QC, filtering, normalization and log transformation, HVG identification, and scaling
#' 
#' @param adata Data downloaded from CytoTRACE, containing two slots: 'exprMatrix' (counts data, features x cells) and 'output' (metadata, which should contain UMAP coordinates)  
#' @param mitotag string denoting prefix for mitochondrial genes
#' @param ribotag string denoting prefix for ribosomal genes
#' @param mtmax threshold value for quality control filtering by percent mitochondrial counts per cell
#' @param s_genes character vector with names of S-phase cell cycle marker genes
#' @param g2m_genes character vector with names of G2M-phase cell cycle marker genes
#' 
#' 
#' @return Seurat object  
process_cyto_adata <-
function(adata, mitotag = "^mt-", ribotag = "^Rp[sl]", mtmax = mtmax, s_genes = s_genes, g2m_genes = g2m_genes) {
  
  #Load downloaded data
  expDat = adata$exprMatrix
  meta = adata$output
  adata = CreateSeuratObject(counts = expDat, project = id, min.cells=5, min.features=200, meta.data = meta)
  adata = subset(adata, cells = rownames(adata@meta.data[!is.na(adata$CytoTRACE),]))
  
  #Filtering
  adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = mitotag)
  adata[["percent.ribo"]] <- PercentageFeatureSet(adata, pattern = ribotag)
  ribogenes <- grep("^Rp[sl]", rownames(adata), value=TRUE, ignore.case = T)
  mtgenes <- grep("^mt--*", rownames(adata), value=TRUE, ignore.case = T)
  adata <- adata[!(rownames(adata) %in% c(mtgenes, "Malat1","MALAT1")),] #removed ribogenes from this expression
  #fmax <- quantile(adata$nFeature_RNA, 0.95)
  #adata <- subset(adata, subset = nFeature_RNA > 200 & nFeature_RNA < fmax & percent.mt < mtmax)
  
  #Normalize, compute HVG, scale
  adata<- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor=10000)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2500)
  VariableFeatures(adata) = VariableFeatures(adata)[!(VariableFeatures(adata) %in% cell_cycle_genes)] #prior was cell cycle genes
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- CellCycleScoring(adata, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
  
  #Use UMAP coordinates given by CytoTRACE
  umap = as.matrix(meta[,c("UMAP1","UMAP2")])
  rownames(umap) = rownames(meta)
  umap = umap[rownames(adata@meta.data),]
  adata@reductions[['umap']] = CreateDimReducObject(embeddings = umap, key = "UMAP_")
  
  return(adata)
}
