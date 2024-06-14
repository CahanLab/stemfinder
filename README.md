stemFinder vignette
================
Kathleen Noller
06/14/2024

# stemFinder

### Single-cell estimation of the extent of differentiation from scRNA-seq data

# 

# 

## Setup

``` r
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("devtools")
```

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/hb/b7nzqfss2_l63s3qz23cqftr0000gp/T//RtmpKztfmS/downloaded_packages

``` r
devtools::install_github("pcahan1/stemfinder")
```

    ## Skipping install of 'stemFinder' from a github remote, the SHA1 (84b438e0) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
library(stemFinder, verbose = F)
```

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: Seurat

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, t

    ## Loading required package: ggplot2

    ## Warning: replacing previous import 'dplyr::select' by 'MASS::select' when
    ## loading 'stemFinder'

    ## Warning: replacing previous import 'dplyr::union' by 'graph::union' when
    ## loading 'stemFinder'

    ## Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
    ## 'stemFinder'

    ## Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
    ## loading 'stemFinder'

## Load query data - Bone marrow from Tabula Muris

#### Query data should be a Seurat object containing a scaled single-cell gene expression matrix

#### Query data must have two metadata columns:

##### Phenotype (character vector of cell type annotations) and Ground_truth (numeric vector of ascending ground truth values denoting extent of differentiation)

###### Note: example data has already been filtered, normalized, and scaled

[Download query data: Tabula Muris bone marrow, 10X
platform](https://cnobjects.s3.amazonaws.com/stemFinder/MurineBoneMarrow10X_GSE109774.rds)

``` r
adata = readRDS("MurineBoneMarrow10X_GSE109774.rds")
head(adata,2)
```

    ##                            orig.ident nCount_RNA nFeature_RNA
    ## X10X_P7_3_AAACCTGAGCATCATC       X10X      12994         3468
    ## X10X_P7_3_AAACCTGCAGAGTGTG       X10X       5437         1764
    ##                                       Phenotype Ground_truth percent.mt
    ## X10X_P7_3_AAACCTGAGCATCATC Monocyte_progenitors            2          0
    ## X10X_P7_3_AAACCTGCAGAGTGTG            Monocytes            3          0
    ##                            percent.ribo    S.Score  G2M.Score Phase
    ## X10X_P7_3_AAACCTGAGCATCATC     22.62583  0.3360127  0.3821197   G2M
    ## X10X_P7_3_AAACCTGCAGAGTGTG     24.03899 -0.1894597 -0.3641231    G1

<img src="figure/stemFinder-unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

## Prepare inputs to stemFinder

``` r
#PCA
adata <- RunPCA(adata, verbose = F)
p1 <- ElbowPlot(adata, ndims = 50)
```

<img src="figure/stemFinder-unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

``` r
#Select PCs based on elbow plot
pcs = 32

#Perform K nearest neighbors
k = round(sqrt(ncol(adata))) #default value of k parameter
adata = FindNeighbors(adata, dims = 1:pcs, k.param = k, verbose = F)
knn = adata@graphs$RNA_nn #KNN matrix

#Select input cell cycle marker gene list
cell_cycle_genes = c(s_genes_mouse, g2m_genes_mouse)[c(s_genes_mouse, g2m_genes_mouse) %in% rownames(adata)] #default G2M + S cell cycle gene list
```

## Run stemFinder

#### Inputs:

##### adata: Seurat object containing scaled gene expression data (features x cells)

##### k: number of nearest neighbors

##### nn: KNN matrix (cells x cells)

##### thresh: threshold for binarizing gene expression data (default = 0)

##### markers: character vector of cell cycle gene names

``` r
adata = run_stemFinder(adata, k = k, nn = knn, thresh = 0, markers = cell_cycle_genes)

head(adata,5) 
```

    ##                            orig.ident nCount_RNA nFeature_RNA
    ## X10X_P7_3_AAACCTGAGCATCATC       X10X      12994         3468
    ## X10X_P7_3_AAACCTGCAGAGTGTG       X10X       5437         1764
    ## X10X_P7_3_AAACCTGGTCGAACAG       X10X       4466         1526
    ## X10X_P7_3_AAACCTGTCACTTCAT       X10X      23852         4043
    ## X10X_P7_3_AAACGGGAGAAGGTTT       X10X       4375          977
    ##                                       Phenotype Ground_truth percent.mt
    ## X10X_P7_3_AAACCTGAGCATCATC Monocyte_progenitors            2          0
    ## X10X_P7_3_AAACCTGCAGAGTGTG            Monocytes            3          0
    ## X10X_P7_3_AAACCTGGTCGAACAG Monocyte_progenitors            2          0
    ## X10X_P7_3_AAACCTGTCACTTCAT     Stem_Progenitors            1          0
    ## X10X_P7_3_AAACGGGAGAAGGTTT         Granulocytes            3          0
    ##                            percent.ribo     S.Score  G2M.Score Phase stemFinder
    ## X10X_P7_3_AAACCTGAGCATCATC    22.625827  0.33601275  0.3821197   G2M  17.450357
    ## X10X_P7_3_AAACCTGCAGAGTGTG    24.038992 -0.18945969 -0.3641231    G1   5.141795
    ## X10X_P7_3_AAACCTGGTCGAACAG    33.631885  0.30172632 -0.1413534     S  15.401308
    ## X10X_P7_3_AAACCTGTCACTTCAT    33.104142 -0.01163238 -0.3062905    G1  18.712842
    ## X10X_P7_3_AAACGGGAGAAGGTTT     2.537143 -0.15402552 -0.1239491    G1   2.988407
    ##                            stemFinder_invert stemFinder_comp
    ## X10X_P7_3_AAACCTGAGCATCATC         0.1417314      0.18967779
    ## X10X_P7_3_AAACCTGCAGAGTGTG         0.7471088      0.05588908
    ## X10X_P7_3_AAACCTGGTCGAACAG         0.2425106      0.16740552
    ## X10X_P7_3_AAACCTGTCACTTCAT         0.0796380      0.20340045
    ## X10X_P7_3_AAACGGGAGAAGGTTT         0.8530199      0.03248268

### The following 3 columns are added to metadata:

##### -Raw stemFinder score (“stemFinder”)

##### -Inverted stemFinder score, corresponding to pseudotime / ground truth (“stemFinder_invert”)

##### -Comparable stemFinder score across datasets (“stemFinder_comp”)

[Check against previously-computed stemFinder results on this
dataset](https://cnobjects.s3.amazonaws.com/stemFinder/bmmc_sF_results.csv)

``` r
sF_scores = read.csv("bmmc_sF_results.csv", row.names = 1)
head(sF_scores,5)
```

    ##                            stemFinder stemFinder_invert stemFinder_comp
    ## X10X_P7_3_AAACCTGAGCATCATC  17.193335        0.15611517      0.18688408
    ## X10X_P7_3_AAACCTGCAGAGTGTG   5.174950        0.74600265      0.05624945
    ## X10X_P7_3_AAACCTGGTCGAACAG  15.202528        0.25382815      0.16524487
    ## X10X_P7_3_AAACCTGTCACTTCAT  18.691755        0.08256958      0.20317125
    ## X10X_P7_3_AAACGGGAGAAGGTTT   3.007756        0.85237303      0.03269300

## Quantify stemFinder performance relative to ground truth

``` r
# Compute stemFinder performance metrics
list_all = compute_performance_single(adata, competitor = F)
```

    ## [1] "Single-cell Spearman Correlation, stemFinder: 0.74"
    ## [1] "AUC, stemFinder: 0.97"
    ## [1] "Phenotypic Spearman correlation, stemFinder: 0.89"

``` r
pct.recov = pct_recover(adata)
```

    ## [1] "Percentage highly potent cells recovered by stemFinder: 82.9573934837093"
    ## [1] "Relative abundance of highly potent cells: 11.6428362999708"

## Optional: compare stemFinder performance to another method

[CytoTRACE and CCAT scores for BMMC query
data](https://cnobjects.s3.amazonaws.com/stemFinder/bmmc_competitor_results.csv)

``` r
#Load pre-computed competitor scores 
comp_scores = read.csv("bmmc_competitor_results.csv", row.names = 1)
head(comp_scores,2)
```

    ##                            CytoTRACE      ccat CytoTRACE_invert ccat_invert
    ## X10X_P7_3_AAACCTGAGCATCATC      2645 0.3818031        0.2281879   0.2388558
    ## X10X_P7_3_AAACCTGCAGAGTGTG      1520 0.2712764        0.5564634   0.4591965

``` r
adata@meta.data = cbind(adata@meta.data, comp_scores) #add to metadata

#Quantify performance
list_all_withcomp = compute_performance_single(adata, competitor = T, comp_id = 'CytoTRACE') 
```

    ## [1] "Single-cell Spearman Correlation, stemFinder: 0.74"
    ## [1] "AUC, stemFinder: 0.97"
    ## [1] "Phenotypic Spearman correlation, stemFinder: 0.89"

``` r
print(list_all_withcomp)
```

    ## $`stemFinder results`
    ## Spearman_SingleCell      Spearman_Pheno                 AUC 
    ##           0.7362893           0.8866655           0.9698956 
    ## 
    ## $`Competitor results`
    ## Spearman_SingleCell      Spearman_Pheno                 AUC 
    ##           0.5712992           0.6307109           0.8847005

## Visualize stemFinder and competitor results

##### UMAP embedding

``` r
p2 <- FeaturePlot(adata, features = c('Ground_truth','stemFinder_invert','CytoTRACE_invert','ccat_invert'), cols = c('blue','red'), ncol = 2)
```

<img src="figure/stemFinder-unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

##### Box plot of inverted stemFinder score

``` r
p3 <- ggplot(adata@meta.data, aes(x = Ground_truth, y = stemFinder_invert)) + geom_point() + geom_boxplot(aes(group = Ground_truth, color = Ground_truth)) + theme_bw() + ggtitle("Inverted stemFinder score vs. Ground truth") + ylab("Inverted stemFinder score") + xlab("Ground truth")
```

<img src="figure/stemFinder-unnamed-chunk-14-1.png" style="display: block; margin: auto;" />
