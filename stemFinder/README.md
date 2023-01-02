---
title: "stemFinder vignette"
author: "Kathleen Noller"
date: "01/01/2023"
output: github_document
knit: (function(inputFile, encoding) {
        Sys.setenv(RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/MacOS/quarto/bin');
        rmarkdown::render(inputFile,
                        encoding=encoding, 
                        output_file=file.path(dirname(inputFile), "README.md")) })
---
## Setup

```r
#install.packages("devtools")
library(devtools)
#devtools::install_github("pcahan1/stemfinder")
library(stemFinder)
```

## Load example data - Mouse BMMC 10X dataset

#### Note: example data has already been filtered, normalized, and scaled
#### Query data must have two metadata columns: Phenotype (string of cell type annotations) and Ground_truth (numeric, ascending ground truth potency values)

[Query data: Murine bone marrow, 10X platform, available on S3](https://cahanlab.s3.amazonaws.com/kathleen.noller/stemFinder/MurineBoneMarrow10X_GSE109774.rds)


```r
adata = readRDS("MurineBoneMarrow10X_GSE109774.rds")
DimPlot(adata, group.by = 'Phenotype')
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

```r
head(adata,2)
```

```
##                            orig.ident nCount_RNA nFeature_RNA            Phenotype
## X10X_P7_3_AAACCTGAGCATCATC       X10X      12994         3468 Monocyte_progenitors
## X10X_P7_3_AAACCTGCAGAGTGTG       X10X       5437         1764            Monocytes
##                            Ground_truth CytoTRACE percent.mt percent.ribo    S.Score
## X10X_P7_3_AAACCTGAGCATCATC            2      2645          0     22.62583  0.3360127
## X10X_P7_3_AAACCTGCAGAGTGTG            3      1520          0     24.03899 -0.1894597
##                             G2M.Score Phase      ccat CytoTRACE_invert ccat_invert
## X10X_P7_3_AAACCTGAGCATCATC  0.3821197   G2M 0.3818031        0.2281879   0.2388558
## X10X_P7_3_AAACCTGCAGAGTGTG -0.3641231    G1 0.2712764        0.5564634   0.4591965
```

## Run stemFinder


```r
#PCA
adata <- RunPCA(adata)
```

```
## PC_ 1 
## Positive:  Ptma, Npm1, Rps2, Snrpf, Hmgb1, Erh, Hsp90aa1, Rpl3, Rpl36a, Serbp1 
## 	   Rps17, Ncl, Rpl12, Rpl14, Rpl15, Hsp90ab1, Set, Ran, Rpl10a, Snrpg 
## 	   Ranbp1, Rps8, Eef1g, Gnb2l1, Rps18, Cdk4, Eif4a1, Rplp1, Rps12, Slc25a5 
## Negative:  Ifitm6, Lcn2, Ngp, Cd177, Ltf, Camp, 1100001G20Rik, Pglyrp1, Mmp9, AA467197 
## 	   Anxa1, Lrg1, Retnlg, Fpr2, Mmp8, Chi3l1, St3gal5, Adpgk, Chi3l3, Itgb2l 
## 	   Gadd45a, Mapk13, G0s2, Actn1, Cd9, Trem3, Fam101b, Cxcr2, Ceacam10, Thbs1 
## PC_ 2 
## Positive:  Car2, Rhd, Ctse, Slc4a1, Gypa, Hemgn, Hmbs, Cldn13, Ermap, Fam55b 
## 	   Aqp1, Slc25a37, Tspan33, Cpox, Tspo2, Ubac1, Klf1, Rhag, Fech, Epb4.2 
## 	   Urod, Blvrb, Trim10, Dhrs11, Slc43a1, Spna1, Pdzk1ip1, Mgst3, Hebp1, Alas2 
## Negative:  Cst3, Plac8, Ly6c2, Ctss, Lgals3, Pld4, Ly86, Clec12a, Emb, Ahnak 
## 	   Nme2, Actb, Cstb, Zyx, Ms4a6c, Ccr2, Ybx1, Rap1b, Rnase6, Cd68 
## 	   H2-DMa, Ramp1, Pgam1, Ccl9, S100a4, Ctsc, Cmtm7, Ifi27l2a, Igsf6, Tgfbi 
## PC_ 3 
## Positive:  Ctss, Rnase6, Pld4, Ly86, Bst2, Ppfia4, Ahnak, Irf8, H2-DMa, Cd74 
## 	   Tcf4, Ctsl, Rgs10, Cd68, Bmyc, Siglech, Irf7, Cd7, Upb1, H2-DMb1 
## 	   Itgb7, Tsc22d1, Clec10a, Lgals1, Cox6a2, Ccnd1, Pltp, Ifi203, Ifi30, Ccr9 
## Negative:  Hmgn2, Trem3, 1100001G20Rik, Pglyrp1, Serpinb1a, Cebpe, Dstn, Ngp, Camp, Lcn2 
## 	   Anxa1, Cd63, Lta4h, Ltf, Lrg1, Mgst2, Fam101b, Adpgk, Cd9, Gca 
## 	   Itgb2l, Cd177, Ap3s1, Fcnb, Chi3l1, Gsr, 1700020L24Rik, Anxa3, Chi3l3, Tmem216 
## PC_ 4 
## Positive:  Prtn3, Ctsg, Mpo, Elane, Gstm1, Nkg7, Cst7, Mt1, Igfbp4, Ms4a3 
## 	   Tuba4a, Gsto1, Gatm, S100a1, Ap3s1, Gsr, Prss57, Dmkn, Lmo2, Cebpa 
## 	   Emb, Svip, Srm, Igsf6, Alas1, Atp8b4, Tfec, Ccl9, Ly6c2, Tspan4 
## Negative:  Pafah1b3, Blnk, Vpreb3, Cnp, Fcrla, Il7r, Ptprcap, Pou2af1, Cd79a, Spib 
## 	   Cd79b, Akap12, Cd72, Myl4, Bach2, 2010001M09Rik, Dnajc7, Atp1b1, Fam129c, Chchd10 
## 	   Bcl7a, Snn, Ly6d, Cecr2, Sdc4, Stambpl1, Cplx2, Cd19, Plekha2, Arl5c 
## PC_ 5 
## Positive:  Siglech, Spint2, Cox6a2, Lair1, Ccr9, Tsc22d1, Lag3, Runx2, Pacsin1, Cd7 
## 	   Prkca, Cd300c, Klk1, Ccnd1, Cybasc3, Tcf4, Ly6d, Fyn, Ctsl, Lefty1 
## 	   Atp2a1, Bcl11a, Paqr5, Nedd4, Il21r, Cmah, Nucb2, Klrd1, Pltp, P2ry14 
## Negative:  S100a4, Tgfbi, Ms4a6c, F13a1, Fn1, Ccr2, Ms4a4c, Ccdc109b, Ms4a6b, Lgals3 
## 	   Gm11428, Ctsc, Ccl9, Ifi30, Csf1r, Naaa, Clec4a3, Clec7a, Klf4, Ccl6 
## 	   Cd302, Mnda, Crip1, Clec4a1, Ifi204, Cx3cr1, Emb, AI607873, Stxbp6, Ifi27l2a
```

```r
ElbowPlot(adata, ndims=50)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

```r
pcs = 32

#K nearest neighbors
k = round(sqrt(ncol(adata)))
adata = FindNeighbors(adata, dims = 1:pcs, k.param = k)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
knn = adata@graphs$RNA_nn
expDat = as.matrix(adata@assays$RNA@scale.data)
cell_cycle_genes = c(s_genes_mouse, g2m_genes_mouse)[c(s_genes_mouse, g2m_genes_mouse) %in% rownames(expDat)]

#Compute single-cell potency
adata = run_stemFinder(adata, k = k, nn = knn, thresh = 0, markers = cell_cycle_genes)
head(adata,2) 
```

```
##                            orig.ident nCount_RNA nFeature_RNA            Phenotype
## X10X_P7_3_AAACCTGAGCATCATC       X10X      12994         3468 Monocyte_progenitors
## X10X_P7_3_AAACCTGCAGAGTGTG       X10X       5437         1764            Monocytes
##                            Ground_truth CytoTRACE percent.mt percent.ribo    S.Score
## X10X_P7_3_AAACCTGAGCATCATC            2      2645          0     22.62583  0.3360127
## X10X_P7_3_AAACCTGCAGAGTGTG            3      1520          0     24.03899 -0.1894597
##                             G2M.Score Phase      ccat CytoTRACE_invert ccat_invert stemFinder
## X10X_P7_3_AAACCTGAGCATCATC  0.3821197   G2M 0.3818031        0.2281879   0.2388558  17.450357
## X10X_P7_3_AAACCTGCAGAGTGTG -0.3641231    G1 0.2712764        0.5564634   0.4591965   5.141795
##                            stemFinder_invert stemFinder_comp
## X10X_P7_3_AAACCTGAGCATCATC         0.1417314      0.18967779
## X10X_P7_3_AAACCTGCAGAGTGTG         0.7471088      0.05588908
```
###3 columns are added to metadata: 
####1.) Raw potency score ("stemFinder")
####2.) Inverted potency score to correspond with pseudotime / ground truth ("stemFinder_invert")
####3.) Comparable potency score across datasets ("stemFinder_comp")

[Check against previously-computed stemFinder results on this dataset](https://cahanlab.s3.amazonaws.com/kathleen.noller/stemFinder/bmmc_sF_results.csv)
      
## Quantify stemFinder performance


```r
# Compute stemFinder performance metrics
list_all = compute_performance_single(adata, competitor = F)
```

```
## [1] "Single-cell Spearman Correlation, stemFinder: 0.74"
## [1] "AUC, stemFinder: 0.97"
## [1] "Phenotypic Spearman correlation, stemFinder: 0.89"
```

```r
pct.recov = pct_recover(adata)
```

```
## [1] "Percentage highly potent cells recovered by stemFinder: 82.9573934837093"
## [1] "Relative abundance of highly potent cells: 11.6428362999708"
```

```r
#Optional: compute performance of stemFinder and an alternative potency calculation method already run on query data 
list_all_withcomp = compute_performance_single(adata, competitor = T, comp_id = 'CytoTRACE')
```

```
## [1] "Single-cell Spearman Correlation, stemFinder: 0.74"
## [1] "AUC, stemFinder: 0.97"
## [1] "Phenotypic Spearman correlation, stemFinder: 0.89"
```

```r
print(list_all_withcomp)
```

```
## $`stemFinder results`
## Spearman_SingleCell      Spearman_Pheno                 AUC 
##           0.7362904           0.8866655           0.9698956 
## 
## $`Competitor results`
## Spearman_SingleCell      Spearman_Pheno                 AUC 
##           0.5712992           0.6307109           0.8847005
```

## Plot single-cell potency scores


```r
FeaturePlot(adata, features = c('Ground_truth','stemFinder_invert','CytoTRACE_invert','ccat_invert'), cols = c('blue','red'))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)
