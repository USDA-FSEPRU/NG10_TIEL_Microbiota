Batch Integration
================
Jayne Wiarda
2022April28

# Overview

Epithelial cell fractions were collected from jejunum, ileum, and cecum
of eight pigs 2 and 4 weeks post-weaning. Pigs were weaned at \~3 weeks
of age, so timepoints are equivalent to \~5 or \~7 weeks of age as well.
Epithelial cell fractions were stained via flow cytometry to detect
intraepithelial T cells (CD3e+ lymphocytes) and associated phenotype
markers (CD4, CD8a, CD8b, gdTCR, CD2, CD16, CD27). Samples were
collected and stained across multiple batches.

We have already gated live CD3e+ lymphocytes from flow cytometry samples
imported into FlowJo (FlowJo, LLC). We exported channel values in .csv
format for all CD3e+ cells of each flow sample. What we want to do now
is merge together all the .csv files with channel values from our
different samples while also incorporating pertinent meta data from our
experiment.

Now we want to integrate data to minimize any batch effects, and we can
visualize both integrated and unintegrated data.

## Load required packages

See session information at bottom for further information.

``` r
library(Spectre)
library(data.table)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     between, first, last

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(harmony)
```

    ## Loading required package: Rcpp

``` r
library(Rcpp)
library(Rtsne)
library(scales)
library(colorRamps)
library(ggthemes)
library(RColorBrewer)
library(ggpointdensity)
```

## Set directory:

``` r
setwd('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/') # set working directory
dir.create('BatchIntegration')
setwd('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/BatchIntegration')
```

## Import merged flow cytometry data:

``` r
cell.dat <- fread("/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/AllCellsMerged_InputForCorrection.csv")
cell.dat
```

    ##          FSC_A FSC_H FSC_W SSC_A SSC_H SSC_W CD8A_APC Viability_eFluor480
    ##       1:   360   335   592   574   385     2      846                 540
    ##       2:   484   438   636   656   436   108      858                 549
    ##       3:   359   332   605   589   393     2      584                 556
    ##       4:   323   304   597   506   350     4      749                 468
    ##       5:   387   363   600   593   402    63      467                 402
    ##      ---                                                                 
    ## 1346771:   330   316   576   483   338     1      737                 470
    ## 1346772:   305   278   583   543   369     1      891                 591
    ## 1346773:   288   277   562   464   330     2      694                 445
    ## 1346774:   306   286   578   473   338     0      685                 493
    ## 1346775:   253   217   568   618   408    76      796                 517
    ##          CD27_FITC CD4_PerCPCy5.5 CD16_BUV496 CD2_BV650 CD8B_PE gdTCR_iFluor594
    ##       1:       319              0          99       761     527             341
    ##       2:       334              0         200       830     610             332
    ##       3:       260              0         584       814     165             554
    ##       4:       149            173         321       761     520             550
    ##       5:       431            251         746       777     207             400
    ##      ---                                                                       
    ## 1346771:       341            427         107       719     694             329
    ## 1346772:       599            223         206       899     745             317
    ## 1346773:       151            334         105       719     420             353
    ## 1346774:        95              0         137       788     534             342
    ## 1346775:       124            485         196       815     758             322
    ##          CD3E_PECy7                                                 FileName
    ##       1:        296 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209
    ##       2:        350 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209
    ##       3:        750 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209
    ##       4:        659 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209
    ##       5:        452 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209
    ##      ---                                                                    
    ## 1346771:        388 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209
    ## 1346772:        421 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209
    ## 1346773:        418 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209
    ## 1346774:        389 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209
    ## 1346775:        469 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209
    ##          FileNo Batch WeeksOfAge  Tissue AnimalID
    ##       1:      1     B          5   cecum       13
    ##       2:      1     B          5   cecum       13
    ##       3:      1     B          5   cecum       13
    ##       4:      1     B          5   cecum       13
    ##       5:      1     B          5   cecum       13
    ##      ---                                         
    ## 1346771:     48     B          5 jejunum        9
    ## 1346772:     48     B          5 jejunum        9
    ## 1346773:     48     B          5 jejunum        9
    ## 1346774:     48     B          5 jejunum        9
    ## 1346775:     48     B          5 jejunum        9

## Define the columns:

``` r
colnames(cell.dat)
```

    ##  [1] "FSC_A"               "FSC_H"               "FSC_W"              
    ##  [4] "SSC_A"               "SSC_H"               "SSC_W"              
    ##  [7] "CD8A_APC"            "Viability_eFluor480" "CD27_FITC"          
    ## [10] "CD4_PerCPCy5.5"      "CD16_BUV496"         "CD2_BV650"          
    ## [13] "CD8B_PE"             "gdTCR_iFluor594"     "CD3E_PECy7"         
    ## [16] "FileName"            "FileNo"              "Batch"              
    ## [19] "WeeksOfAge"          "Tissue"              "AnimalID"

``` r
cell.dat$SampleID <- paste0(cell.dat$AnimalID, cell.dat$Tissue) # create a column for sample IDs (animal ID + tissue)
cluster.cols <- names(cell.dat)[c(7, 9:14)] # define columns to use for dimensionality reduction & all clustering/integration
cluster.cols
```

    ## [1] "CD8A_APC"        "CD27_FITC"       "CD4_PerCPCy5.5"  "CD16_BUV496"    
    ## [5] "CD2_BV650"       "CD8B_PE"         "gdTCR_iFluor594"

## Subset data

``` r
nrow(cell.dat) # how many total cells do we have?
```

    ## [1] 1346775

``` r
min(table(cell.dat$SampleID)) # how many cells were in our smallest sample?
```

    ## [1] 4471

``` r
max(table(cell.dat$SampleID)) # how many cells were in our smallest sample?
```

    ## [1] 106951

``` r
cell.dat <- do.subsample(cell.dat, 
                         targets = rep(min(table(cell.dat$SampleID)), length(unique(cell.dat$SampleID))), # subset to a more manageable and more equal number of cells from each sample, specified as the number of cells in our smallest sample
                         divide.by = 'SampleID') # subset so that every sample has the same number of cells as was present in our smallest sample
nrow(cell.dat) # how many total cells now?
```

    ## [1] 214608

## Perform integration

Perform integration across batches with Harmony data alignment algorithm

``` r
# We opt to use all samples for harmonizing data since we don't have a true 'control' treatment group...
colnames(cell.dat)
```

    ##  [1] "FSC_A"               "FSC_H"               "FSC_W"              
    ##  [4] "SSC_A"               "SSC_H"               "SSC_W"              
    ##  [7] "CD8A_APC"            "Viability_eFluor480" "CD27_FITC"          
    ## [10] "CD4_PerCPCy5.5"      "CD16_BUV496"         "CD2_BV650"          
    ## [13] "CD8B_PE"             "gdTCR_iFluor594"     "CD3E_PECy7"         
    ## [16] "FileName"            "FileNo"              "Batch"              
    ## [19] "WeeksOfAge"          "Tissue"              "AnimalID"           
    ## [22] "SampleID"

``` r
cell.dat <- run.harmony(cell.dat, # this step can take a few minutes...
                        batch.col = 'Batch', 
                        align.cols = cluster.cols)
```

    ## run.harmony - preparing data (1/3)

    ## run.harmony - running harmony (2/3)

    ## run.harmony - harmony complete, returning data (3/3)

``` r
colnames(cell.dat) # note new columns with _aligned cell markers...this is batch-corrected data values we will use from now on
```

    ##  [1] "FSC_A"                   "FSC_H"                  
    ##  [3] "FSC_W"                   "SSC_A"                  
    ##  [5] "SSC_H"                   "SSC_W"                  
    ##  [7] "CD8A_APC"                "Viability_eFluor480"    
    ##  [9] "CD27_FITC"               "CD4_PerCPCy5.5"         
    ## [11] "CD16_BUV496"             "CD2_BV650"              
    ## [13] "CD8B_PE"                 "gdTCR_iFluor594"        
    ## [15] "CD3E_PECy7"              "FileName"               
    ## [17] "FileNo"                  "Batch"                  
    ## [19] "WeeksOfAge"              "Tissue"                 
    ## [21] "AnimalID"                "SampleID"               
    ## [23] "CD8A_APC_aligned"        "CD27_FITC_aligned"      
    ## [25] "CD4_PerCPCy5.5_aligned"  "CD16_BUV496_aligned"    
    ## [27] "CD2_BV650_aligned"       "CD8B_PE_aligned"        
    ## [29] "gdTCR_iFluor594_aligned"

## Save data

``` r
fwrite(cell.dat, "/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/BatchIntegration/AllCellsMerged_SubsetCells_HarmonyIntegratedData.csv") 
```

## Create t-SNE plots to visualize data before and after integration

First subset data again to make the process faster:

``` r
sub.all <- do.subsample(cell.dat, 
                        targets = rep(250, length(unique(cell.dat$SampleID))), 
                        divide.by = 'SampleID') # subset 250 random cells from each sample
rm(cell.dat) # free up space
```

Generate t-SNE coordinates for data before integration:

``` r
colnames(sub.all)
```

    ##  [1] "FSC_A"                   "FSC_H"                  
    ##  [3] "FSC_W"                   "SSC_A"                  
    ##  [5] "SSC_H"                   "SSC_W"                  
    ##  [7] "CD8A_APC"                "Viability_eFluor480"    
    ##  [9] "CD27_FITC"               "CD4_PerCPCy5.5"         
    ## [11] "CD16_BUV496"             "CD2_BV650"              
    ## [13] "CD8B_PE"                 "gdTCR_iFluor594"        
    ## [15] "CD3E_PECy7"              "FileName"               
    ## [17] "FileNo"                  "Batch"                  
    ## [19] "WeeksOfAge"              "Tissue"                 
    ## [21] "AnimalID"                "SampleID"               
    ## [23] "CD8A_APC_aligned"        "CD27_FITC_aligned"      
    ## [25] "CD4_PerCPCy5.5_aligned"  "CD16_BUV496_aligned"    
    ## [27] "CD2_BV650_aligned"       "CD8B_PE_aligned"        
    ## [29] "gdTCR_iFluor594_aligned"

``` r
set.seed(123)
tsne.out <- Rtsne(sub.all[,c(7,9:14)], # use only these columns to generate t-SNE coordinates
                  perplexity = 15) # toggle with perplexity to modify the spread of cells in t-SNE plot
sub.all$tsne1_Uncorrected <- tsne.out$Y[,1]
sub.all$tsne2_Uncorrected <- tsne.out$Y[,2]
rm(tsne.out) # free up space
```

Generate t-SNE coordinates for data after integration:

``` r
colnames(sub.all)
```

    ##  [1] "FSC_A"                   "FSC_H"                  
    ##  [3] "FSC_W"                   "SSC_A"                  
    ##  [5] "SSC_H"                   "SSC_W"                  
    ##  [7] "CD8A_APC"                "Viability_eFluor480"    
    ##  [9] "CD27_FITC"               "CD4_PerCPCy5.5"         
    ## [11] "CD16_BUV496"             "CD2_BV650"              
    ## [13] "CD8B_PE"                 "gdTCR_iFluor594"        
    ## [15] "CD3E_PECy7"              "FileName"               
    ## [17] "FileNo"                  "Batch"                  
    ## [19] "WeeksOfAge"              "Tissue"                 
    ## [21] "AnimalID"                "SampleID"               
    ## [23] "CD8A_APC_aligned"        "CD27_FITC_aligned"      
    ## [25] "CD4_PerCPCy5.5_aligned"  "CD16_BUV496_aligned"    
    ## [27] "CD2_BV650_aligned"       "CD8B_PE_aligned"        
    ## [29] "gdTCR_iFluor594_aligned" "tsne1_Uncorrected"      
    ## [31] "tsne2_Uncorrected"

``` r
set.seed(123)
tsne.out <- Rtsne(sub.all[,c(23:29)], # use only these columns to generate t-SNE coordinates
                  perplexity = 15) # toggle with perplexity to modify the spread of cells in t-SNE plot
sub.all$tsne1_Corrected <- tsne.out$Y[,1]
sub.all$tsne2_Corrected <- tsne.out$Y[,2]
rm(tsne.out) # free up space
```

Plot cell markers before integration:

``` r
colnames(sub.all)
```

    ##  [1] "FSC_A"                   "FSC_H"                  
    ##  [3] "FSC_W"                   "SSC_A"                  
    ##  [5] "SSC_H"                   "SSC_W"                  
    ##  [7] "CD8A_APC"                "Viability_eFluor480"    
    ##  [9] "CD27_FITC"               "CD4_PerCPCy5.5"         
    ## [11] "CD16_BUV496"             "CD2_BV650"              
    ## [13] "CD8B_PE"                 "gdTCR_iFluor594"        
    ## [15] "CD3E_PECy7"              "FileName"               
    ## [17] "FileNo"                  "Batch"                  
    ## [19] "WeeksOfAge"              "Tissue"                 
    ## [21] "AnimalID"                "SampleID"               
    ## [23] "CD8A_APC_aligned"        "CD27_FITC_aligned"      
    ## [25] "CD4_PerCPCy5.5_aligned"  "CD16_BUV496_aligned"    
    ## [27] "CD2_BV650_aligned"       "CD8B_PE_aligned"        
    ## [29] "gdTCR_iFluor594_aligned" "tsne1_Uncorrected"      
    ## [31] "tsne2_Uncorrected"       "tsne1_Corrected"        
    ## [33] "tsne2_Corrected"

``` r
make.multi.plot(sub.all, 
                'tsne1_Uncorrected', 
                'tsne2_Uncorrected', 
                colnames(sub.all[,c(1:15)]), # note: some of these parameters were not used to perform t-SNE dimensionality reduction, but we can still see where they show up in the plot. 
                colours = 'inferno',
                figure.title = 'AllCellMarkers_BeforeIntegration_250downsample',
                path = '/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/BatchIntegration')
```

    ## Loading required package: ggplot2

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-7.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-8.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-9.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-10.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-11.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-12.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-13.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-14.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-15.png)<!-- -->

    ## Check your working directory for a new .png called 'AllCellMarkers_BeforeIntegration_250downsample.png'

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-11-16.png)<!-- -->

Plot cell markers after integration:

``` r
#colnames(sub.all)
make.multi.plot(sub.all, 
                'tsne1_Corrected', 
                'tsne2_Corrected', 
                colnames(sub.all[,c(1:6, 8, 15, 23:29)]), # note: some of these parameters were not used to perform t-SNE dimensionality reduction, but we can still see where they show up in the plot. Some parameters also don't have aligned values, but we can still view the unaligned data for them.
                colours = 'inferno',
                figure.title = 'AllCellMarkers_AfterIntegration_250downsample',
                path = '/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/BatchIntegration')
```

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-7.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-8.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-9.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-10.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-11.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-12.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-13.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-14.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-15.png)<!-- -->

    ## Check your working directory for a new .png called 'AllCellMarkers_AfterIntegration_250downsample.png'

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-12-16.png)<!-- -->

Plot some variables from meta data before integration:

``` r
#colnames(sub.all)
sub.all$combo <- paste(sub.all$Tissue, sub.all$WeeksOfAge, sep = '_')
make.multi.plot(sub.all, 
                'tsne1_Uncorrected', 
                'tsne2_Uncorrected', 
                c('Batch', 'WeeksOfAge', 'Tissue', 'combo'), 
                col.type = 'factor', 
                dot.size = 0.2,
                figure.title = 'MetaData_BeforeIntegration_250downsample',
                path = '/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/BatchIntegration')
```

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-13-4.png)<!-- -->

    ## Check your working directory for a new .png called 'MetaData_BeforeIntegration_250downsample.png'

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-13-5.png)<!-- -->

Plot some variables from meta data after integration:

``` r
#colnames(sub.all)
#sub.all$combo <- paste(sub.all$Tissue, sub.all$WeeksOfAge, sep = '_')
make.multi.plot(sub.all, 
                'tsne1_Corrected', 
                'tsne2_Corrected', 
                c('Batch', 'WeeksOfAge', 'Tissue', 'combo'), 
                col.type = 'factor', 
                dot.size = 0.2,
                figure.title = 'MetaData_AfterIntegration_250downsample',
                path = '/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/BatchIntegration')
```

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->

    ## Check your working directory for a new .png called 'MetaData_AfterIntegration_250downsample.png'

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-14-5.png)<!-- -->

Plot cell density before integration:

``` r
make.colour.plot(sub.all, 
                 'tsne1_Uncorrected', 
                 'tsne2_Uncorrected', 
                 filename = 'CellDensity_BeforeIntegration_250downsample.png',
                 colours = 'inferno',
                path = '/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/BatchIntegration')
```

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Plot cell density after integration:

``` r
make.colour.plot(sub.all, 
                 'tsne1_Corrected', 
                 'tsne2_Corrected', 
                 filename = 'CellDensity_AfterIntegration_250downsample.png',
                 colours = 'inferno',
                path = '/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/BatchIntegration')
```

![](02_BatchIntegration_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## Session information

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggplot2_3.3.5        ggpointdensity_0.1.0 RColorBrewer_1.1-2  
    ##  [4] ggthemes_4.2.4       colorRamps_2.3       scales_1.1.1        
    ##  [7] Rtsne_0.15           harmony_0.1.0        Rcpp_1.0.7          
    ## [10] dplyr_1.0.7          data.table_1.14.0    Spectre_0.4.1       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.1  xfun_0.26         purrr_0.3.4       colorspace_2.0-2 
    ##  [5] vctrs_0.3.8       generics_0.1.0    htmltools_0.5.2   viridisLite_0.4.0
    ##  [9] yaml_2.2.1        utf8_1.2.2        rlang_0.4.11      pillar_1.6.2     
    ## [13] glue_1.4.2        withr_2.4.2       DBI_1.1.1         lifecycle_1.0.0  
    ## [17] stringr_1.4.0     munsell_0.5.0     gtable_0.3.0      ragg_1.2.0       
    ## [21] codetools_0.2-18  evaluate_0.14     labeling_0.4.2    knitr_1.34       
    ## [25] fastmap_1.1.0     fansi_0.5.0       highr_0.9         farver_2.1.0     
    ## [29] systemfonts_1.0.3 textshaping_0.3.6 gridExtra_2.3     digest_0.6.27    
    ## [33] stringi_1.7.4     cowplot_1.1.1     grid_4.1.3        tools_4.1.3      
    ## [37] magrittr_2.0.1    tibble_3.1.4      crayon_1.4.1      tidyr_1.1.3      
    ## [41] pkgconfig_2.0.3   ellipsis_0.3.2    viridis_0.6.2     assertthat_0.2.1 
    ## [45] rmarkdown_2.11    R6_2.5.1          compiler_4.1.3
