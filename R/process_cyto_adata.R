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
