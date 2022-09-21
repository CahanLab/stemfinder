stemFinder vignette
================
Kathleen Noller
9/13/2022

## Setup

``` r
library(Seurat)
```

    ## Attaching SeuratObject

    ## Attaching sp

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
load("./data/s_genes_mouse.rda")
load("./data/g2m_genes_mouse.rda")
source("./R/run_stemFinder.R")
source("./R/compute_performance.R")
```

## Load example data - Mouse BMMC 10X dataset

#### Note: example data has already been filtered, normalized, and scaled

#### Query data must have two metadata columns: Phenotype (string of cell type annotations) and Ground\_truth (numeric, ground truth potency values)

[Query data available on S3](https://s3.console.aws.amazon.com/s3/object/cahanlab/kathleen.noller/stemFinder/Mouse10X_BoneMarrow_GSE109774_stemFinder.rds?region=us-east-1&tab=overview )

``` r
DimPlot(adata, group.by = 'Phenotype')
```

![](Vignette_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
head(adata)
```

    ##                            orig.ident nCount_RNA nFeature_RNA
    ## X10X_P7_3_AAACCTGAGCATCATC       X10X      12994         3468
    ## X10X_P7_3_AAACCTGCAGAGTGTG       X10X       5437         1764
    ## X10X_P7_3_AAACCTGGTCGAACAG       X10X       4466         1526
    ## X10X_P7_3_AAACCTGTCACTTCAT       X10X      23852         4043
    ## X10X_P7_3_AAACGGGAGAAGGTTT       X10X       4375          977
    ## X10X_P7_3_AAACGGGAGTAGTGCG       X10X       5451         1644
    ## X10X_P7_3_AAACGGGCACGACTCG       X10X      30906         3836
    ## X10X_P7_3_AAACGGGGTAGCACGA       X10X       6076         1940
    ## X10X_P7_3_AAAGATGCACTCGACG       X10X       8027         2981
    ## X10X_P7_3_AAAGATGCATGGATGG       X10X       6458         1487
    ##                                          Phenotype Ground_truth CytoTRACE
    ## X10X_P7_3_AAACCTGAGCATCATC    Monocyte_progenitors            2      2645
    ## X10X_P7_3_AAACCTGCAGAGTGTG               Monocytes            3      1520
    ## X10X_P7_3_AAACCTGGTCGAACAG    Monocyte_progenitors            2      2957
    ## X10X_P7_3_AAACCTGTCACTTCAT        Stem_Progenitors            1      3086
    ## X10X_P7_3_AAACGGGAGAAGGTTT            Granulocytes            3       390
    ## X10X_P7_3_AAACGGGAGTAGTGCG               Monocytes            3      1707
    ## X10X_P7_3_AAACGGGCACGACTCG Granulocyte_progenitors            2      1389
    ## X10X_P7_3_AAACGGGGTAGCACGA               Monocytes            3      1620
    ## X10X_P7_3_AAAGATGCACTCGACG              Immature_B            2      3053
    ## X10X_P7_3_AAAGATGCATGGATGG Granulocyte_progenitors            2      1022
    ##                            percent.mt percent.ribo     S.Score  G2M.Score Phase
    ## X10X_P7_3_AAACCTGAGCATCATC          0    22.625827  0.33601275  0.3821197   G2M
    ## X10X_P7_3_AAACCTGCAGAGTGTG          0    24.038992 -0.18945969 -0.3641231    G1
    ## X10X_P7_3_AAACCTGGTCGAACAG          0    33.631885  0.30172632 -0.1413534     S
    ## X10X_P7_3_AAACCTGTCACTTCAT          0    33.104142 -0.01163238 -0.3062905    G1
    ## X10X_P7_3_AAACGGGAGAAGGTTT          0     2.537143 -0.15402552 -0.1239491    G1
    ## X10X_P7_3_AAACGGGAGTAGTGCG          0    30.544854 -0.21395371 -0.3532216    G1
    ## X10X_P7_3_AAACGGGCACGACTCG          0    12.738627  0.02354409  0.4392959   G2M
    ## X10X_P7_3_AAACGGGGTAGCACGA          0    22.728769 -0.22305128 -0.3560751    G1
    ## X10X_P7_3_AAAGATGCACTCGACG          0    11.847515  0.44817941  0.9243623   G2M
    ## X10X_P7_3_AAAGATGCATGGATGG          0     7.742335  0.27757475  0.1743533     S
    ##                                 ccat CytoTRACE_invert ccat_invert stemFinder
    ## X10X_P7_3_AAACCTGAGCATCATC 0.3818031       0.22818792   0.2388558  17.193335
    ## X10X_P7_3_AAACCTGCAGAGTGTG 0.2712764       0.55646338   0.4591965   5.174950
    ## X10X_P7_3_AAACCTGGTCGAACAG 0.3050967       0.13714619   0.3917739  15.202528
    ## X10X_P7_3_AAACCTGTCACTTCAT 0.3830006       0.09950394   0.2364684  18.691755
    ## X10X_P7_3_AAACGGGAGAAGGTTT 0.1473567       0.88619784   0.7062368   3.007756
    ## X10X_P7_3_AAACGGGAGTAGTGCG 0.2461174       0.50189670   0.5093523   3.617926
    ## X10X_P7_3_AAACGGGCACGACTCG 0.3436792       0.59468923   0.3148577  14.930192
    ## X10X_P7_3_AAACGGGGTAGCACGA 0.2639813       0.52728334   0.4737396   4.753806
    ## X10X_P7_3_AAAGATGCACTCGACG 0.3629157       0.10913335   0.2765088  12.912956
    ## X10X_P7_3_AAAGATGCATGGATGG 0.2559489       0.70177998   0.4897526  17.345016
    ##                            stemFinder_invert stemFinder_comp
    ## X10X_P7_3_AAACCTGAGCATCATC        0.15611517      0.18688408
    ## X10X_P7_3_AAACCTGCAGAGTGTG        0.74600265      0.05624945
    ## X10X_P7_3_AAACCTGGTCGAACAG        0.25382815      0.16524487
    ## X10X_P7_3_AAACCTGTCACTTCAT        0.08256958      0.20317125
    ## X10X_P7_3_AAACGGGAGAAGGTTT        0.85237303      0.03269300
    ## X10X_P7_3_AAACGGGAGTAGTGCG        0.82242464      0.03932528
    ## X10X_P7_3_AAACGGGCACGACTCG        0.26719495      0.16228470
    ## X10X_P7_3_AAACGGGGTAGCACGA        0.76667325      0.05167181
    ## X10X_P7_3_AAAGATGCACTCGACG        0.36620513      0.14035822
    ## X10X_P7_3_AAAGATGCATGGATGG        0.14867037      0.18853278

## Run stemFinder

``` r
#PCA
adata <- RunPCA(adata)
```

    ## PC_ 1 
    ## Positive:  Ptma, Npm1, Rps2, Snrpf, Hmgb1, Erh, Hsp90aa1, Rpl3, Rpl36a, Serbp1 
    ##     Rps17, Ncl, Rpl12, Rpl14, Rpl15, Hsp90ab1, Set, Ran, Rpl10a, Snrpg 
    ##     Ranbp1, Rps8, Eef1g, Gnb2l1, Rps18, Cdk4, Eif4a1, Rplp1, Rps12, Slc25a5 
    ## Negative:  Ifitm6, Lcn2, Ngp, Cd177, Ltf, Camp, 1100001G20Rik, Pglyrp1, Mmp9, AA467197 
    ##     Anxa1, Lrg1, Retnlg, Fpr2, Mmp8, Chi3l1, St3gal5, Adpgk, Chi3l3, Itgb2l 
    ##     Gadd45a, Mapk13, G0s2, Actn1, Cd9, Trem3, Fam101b, Cxcr2, Ceacam10, Thbs1 
    ## PC_ 2 
    ## Positive:  Car2, Rhd, Ctse, Slc4a1, Gypa, Hemgn, Hmbs, Cldn13, Ermap, Fam55b 
    ##     Aqp1, Slc25a37, Tspan33, Cpox, Tspo2, Ubac1, Klf1, Rhag, Fech, Epb4.2 
    ##     Urod, Blvrb, Trim10, Dhrs11, Slc43a1, Spna1, Pdzk1ip1, Mgst3, Hebp1, Alas2 
    ## Negative:  Cst3, Plac8, Ly6c2, Ctss, Lgals3, Pld4, Ly86, Clec12a, Emb, Ahnak 
    ##     Nme2, Actb, Cstb, Zyx, Ms4a6c, Ccr2, Ybx1, Rap1b, Rnase6, Cd68 
    ##     H2-DMa, Ramp1, Pgam1, Ccl9, S100a4, Ctsc, Cmtm7, Ifi27l2a, Igsf6, Tgfbi 
    ## PC_ 3 
    ## Positive:  Ctss, Rnase6, Pld4, Ly86, Bst2, Ppfia4, Ahnak, Irf8, H2-DMa, Cd74 
    ##     Tcf4, Ctsl, Rgs10, Cd68, Bmyc, Siglech, Irf7, Cd7, Upb1, H2-DMb1 
    ##     Itgb7, Tsc22d1, Clec10a, Lgals1, Cox6a2, Ccnd1, Pltp, Ifi203, Ifi30, Ccr9 
    ## Negative:  Hmgn2, Trem3, 1100001G20Rik, Pglyrp1, Serpinb1a, Cebpe, Dstn, Ngp, Camp, Lcn2 
    ##     Anxa1, Cd63, Lta4h, Ltf, Lrg1, Mgst2, Fam101b, Adpgk, Cd9, Gca 
    ##     Itgb2l, Cd177, Ap3s1, Fcnb, Chi3l1, Gsr, 1700020L24Rik, Anxa3, Chi3l3, Tmem216 
    ## PC_ 4 
    ## Positive:  Prtn3, Ctsg, Mpo, Elane, Gstm1, Nkg7, Cst7, Mt1, Igfbp4, Ms4a3 
    ##     Tuba4a, Gsto1, Gatm, S100a1, Ap3s1, Gsr, Prss57, Dmkn, Lmo2, Cebpa 
    ##     Emb, Svip, Srm, Igsf6, Alas1, Atp8b4, Tfec, Ccl9, Ly6c2, Tspan4 
    ## Negative:  Pafah1b3, Blnk, Vpreb3, Cnp, Fcrla, Il7r, Ptprcap, Pou2af1, Cd79a, Spib 
    ##     Cd79b, Akap12, Cd72, Myl4, Bach2, 2010001M09Rik, Dnajc7, Atp1b1, Fam129c, Chchd10 
    ##     Bcl7a, Snn, Ly6d, Cecr2, Sdc4, Stambpl1, Cplx2, Cd19, Plekha2, Arl5c 
    ## PC_ 5 
    ## Positive:  Siglech, Spint2, Cox6a2, Lair1, Ccr9, Tsc22d1, Lag3, Runx2, Pacsin1, Cd7 
    ##     Prkca, Cd300c, Klk1, Ccnd1, Cybasc3, Tcf4, Ly6d, Fyn, Ctsl, Lefty1 
    ##     Atp2a1, Bcl11a, Paqr5, Nedd4, Il21r, Cmah, Nucb2, Klrd1, Pltp, P2ry14 
    ## Negative:  S100a4, Tgfbi, Ms4a6c, F13a1, Fn1, Ccr2, Ms4a4c, Ccdc109b, Ms4a6b, Lgals3 
    ##     Gm11428, Ctsc, Ccl9, Ifi30, Csf1r, Naaa, Clec4a3, Clec7a, Klf4, Ccl6 
    ##     Cd302, Mnda, Crip1, Clec4a1, Ifi204, Cx3cr1, Emb, AI607873, Stxbp6, Ifi27l2a

``` r
ElbowPlot(adata, ndims=50)
```

![](Vignette_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
pcs = 32

#K nearest neighbors
k = round(sqrt(ncol(adata)))
adata = FindNeighbors(adata, dims = 1:pcs, k.param = k)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
knn = adata@graphs$RNA_nn
expDat = as.matrix(adata@assays$RNA@scale.data)
cell_cycle_genes = c(s_genes_mouse, g2m_genes_mouse)[c(s_genes_mouse, g2m_genes_mouse) %in% rownames(expDat)]

#Compute single-cell potency
adata = run_stemFinder(adata, expDat, k = k, nn = knn, thresh = 0, markers = cell_cycle_genes)
head(adata) 
```

    ##                            orig.ident nCount_RNA nFeature_RNA
    ## X10X_P7_3_AAACCTGAGCATCATC       X10X      12994         3468
    ## X10X_P7_3_AAACCTGCAGAGTGTG       X10X       5437         1764
    ## X10X_P7_3_AAACCTGGTCGAACAG       X10X       4466         1526
    ## X10X_P7_3_AAACCTGTCACTTCAT       X10X      23852         4043
    ## X10X_P7_3_AAACGGGAGAAGGTTT       X10X       4375          977
    ## X10X_P7_3_AAACGGGAGTAGTGCG       X10X       5451         1644
    ## X10X_P7_3_AAACGGGCACGACTCG       X10X      30906         3836
    ## X10X_P7_3_AAACGGGGTAGCACGA       X10X       6076         1940
    ## X10X_P7_3_AAAGATGCACTCGACG       X10X       8027         2981
    ## X10X_P7_3_AAAGATGCATGGATGG       X10X       6458         1487
    ##                                          Phenotype Ground_truth CytoTRACE
    ## X10X_P7_3_AAACCTGAGCATCATC    Monocyte_progenitors            2      2645
    ## X10X_P7_3_AAACCTGCAGAGTGTG               Monocytes            3      1520
    ## X10X_P7_3_AAACCTGGTCGAACAG    Monocyte_progenitors            2      2957
    ## X10X_P7_3_AAACCTGTCACTTCAT        Stem_Progenitors            1      3086
    ## X10X_P7_3_AAACGGGAGAAGGTTT            Granulocytes            3       390
    ## X10X_P7_3_AAACGGGAGTAGTGCG               Monocytes            3      1707
    ## X10X_P7_3_AAACGGGCACGACTCG Granulocyte_progenitors            2      1389
    ## X10X_P7_3_AAACGGGGTAGCACGA               Monocytes            3      1620
    ## X10X_P7_3_AAAGATGCACTCGACG              Immature_B            2      3053
    ## X10X_P7_3_AAAGATGCATGGATGG Granulocyte_progenitors            2      1022
    ##                            percent.mt percent.ribo     S.Score  G2M.Score Phase
    ## X10X_P7_3_AAACCTGAGCATCATC          0    22.625827  0.33601275  0.3821197   G2M
    ## X10X_P7_3_AAACCTGCAGAGTGTG          0    24.038992 -0.18945969 -0.3641231    G1
    ## X10X_P7_3_AAACCTGGTCGAACAG          0    33.631885  0.30172632 -0.1413534     S
    ## X10X_P7_3_AAACCTGTCACTTCAT          0    33.104142 -0.01163238 -0.3062905    G1
    ## X10X_P7_3_AAACGGGAGAAGGTTT          0     2.537143 -0.15402552 -0.1239491    G1
    ## X10X_P7_3_AAACGGGAGTAGTGCG          0    30.544854 -0.21395371 -0.3532216    G1
    ## X10X_P7_3_AAACGGGCACGACTCG          0    12.738627  0.02354409  0.4392959   G2M
    ## X10X_P7_3_AAACGGGGTAGCACGA          0    22.728769 -0.22305128 -0.3560751    G1
    ## X10X_P7_3_AAAGATGCACTCGACG          0    11.847515  0.44817941  0.9243623   G2M
    ## X10X_P7_3_AAAGATGCATGGATGG          0     7.742335  0.27757475  0.1743533     S
    ##                                 ccat CytoTRACE_invert ccat_invert stemFinder
    ## X10X_P7_3_AAACCTGAGCATCATC 0.3818031       0.22818792   0.2388558  17.847745
    ## X10X_P7_3_AAACCTGCAGAGTGTG 0.2712764       0.55646338   0.4591965   6.359092
    ## X10X_P7_3_AAACCTGGTCGAACAG 0.3050967       0.13714619   0.3917739  15.843149
    ## X10X_P7_3_AAACCTGTCACTTCAT 0.3830006       0.09950394   0.2364684  19.008905
    ## X10X_P7_3_AAACGGGAGAAGGTTT 0.1473567       0.88619784   0.7062368   4.313703
    ## X10X_P7_3_AAACGGGAGTAGTGCG 0.2461174       0.50189670   0.5093523   5.033036
    ## X10X_P7_3_AAACGGGCACGACTCG 0.3436792       0.59468923   0.3148577  15.667912
    ## X10X_P7_3_AAACGGGGTAGCACGA 0.2639813       0.52728334   0.4737396   6.022407
    ## X10X_P7_3_AAAGATGCACTCGACG 0.3629157       0.10913335   0.2765088  13.696064
    ## X10X_P7_3_AAAGATGCATGGATGG 0.2559489       0.70177998   0.4897526  17.760414
    ##                            stemFinder_invert stemFinder_comp
    ## X10X_P7_3_AAACCTGAGCATCATC        0.12888390      0.19399723
    ## X10X_P7_3_AAACCTGCAGAGTGTG        0.68962423      0.06912057
    ## X10X_P7_3_AAACCTGGTCGAACAG        0.22672462      0.17220814
    ## X10X_P7_3_AAACCTGTCACTTCAT        0.07220976      0.20661854
    ## X10X_P7_3_AAACGGGAGAAGGTTT        0.78945597      0.04688808
    ## X10X_P7_3_AAACGGGAGTAGTGCG        0.75434661      0.05470692
    ## X10X_P7_3_AAACGGGCACGACTCG        0.23527762      0.17030339
    ## X10X_P7_3_AAACGGGGTAGCACGA        0.70605721      0.06546095
    ## X10X_P7_3_AAAGATGCACTCGACG        0.33151991      0.14887026
    ## X10X_P7_3_AAAGATGCATGGATGG        0.13314638      0.19304797

``` r
    ##3 columns are added to metadata: 
      ###1. raw potency score ("stemFinder")
      ###2. comparable potency score across datasets ("stemFinder_comp")
      ###3. inverted potency score to correspond with pseudotime / ground truth ("stemFinder_invert")
```

## Quantify performance

``` r
list_all = compute_performance(adata, d = 1, id = "Mouse BMMC", list_all = c())
```

    ## [1] "Single-cell Spearman Correlation, stemFinder: 0.736689385728266"
    ## [1] "Single-cell Spearman Correlation, CCAT: 0.675567206638321"
    ## [1] "Single-cell Spearman Correlation, CytoTRACE: 0.571299186770523"

``` r
pct.recov = pct_recover(adata)
```

    ## [1] "Percentage highly potent cells recovered by stemFinder: 83.2080200501253"
    ## [1] "Relative abundance of highly potent cells: 11.6428362999708"

## Plot single-cell potency scores

``` r
FeaturePlot(adata, features = c('Ground_truth','stemFinder_invert','CytoTRACE_invert','ccat_invert'), cols = c('blue','red'))
```

![](Vignette_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
