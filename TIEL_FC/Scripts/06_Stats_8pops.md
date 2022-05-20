Statistical analysis of 8 T-IEL populations
================
Jayne Wiarda
2022April30

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

We next integrated data to lessen batch effects. While visualizing
integrated data, we noted some CD4 T cells and what appears to be some
cellular debri that we elected to filter out. We are now left with an
equal number of integrated, filtered cells from each sample. At this
point, we elected to use manual gating to define biologically-relevant T
cell populations for further analysis.

We annotated cells into eight major cell populations as follows:

-   CD8a+CD8b- ab T
-   CD8a+CD8b+ ab T
-   CD8a-CD8b+ ab T
-   CD8a-CD8b- ab T
-   CD2-CD8a+ gd T
-   CD2+CD8a+ gd T
-   CD2+CD8a- gd T
-   CD2-CD8a- gd T

The objective of the below analyses is to further investigate abundances
and statistical comparisons for different T-IEL populations across age
and tissue.

## Load required packages

See session information at bottom for further information.

``` r
library(readxl)
library(rstatix)
library(Spectre)
library(dplyr)
library(writexl)
```

## Run stats

Read in the counts data for the 8 subsets of T-IELs:

``` r
dat <- as.data.frame(read_excel('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/CellGatingAnnotation/TotalCellCountsBySample_8pops.xlsx'))
rownames(dat) <- dat$Var2
```

Read in meta data & organize further:

``` r
meta.dat <- read_excel("/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/NG10_MetaData_FC.xlsx") # read in meta data
meta.dat$Var2 <- paste0(meta.dat$AnimalID, meta.dat$Tissue)
```

Convert counts data to percentages of total T-IELs:

``` r
colnames(dat)
```

    ## [1] "Var2"              "AB_CD8AnegCD8Bneg" "AB_CD8AnegCD8Bpos"
    ## [4] "AB_CD8AposCD8Bneg" "AB_CD8AposCD8Bpos" "GD_CD2negCD8Aneg" 
    ## [7] "GD_CD2negCD8Apos"  "GD_CD2posCD8Aneg"  "GD_CD2posCD8Apos"

``` r
rowSums(dat[2:9])
```

    ##   13cecum   13ileum 13jejunum   17cecum   17ileum 17jejunum    1cecum    1ileum 
    ##      3451      3451      3451      3451      3451      3451      3451      3451 
    ##  1jejunum   21cecum   21ileum 21jejunum   25cecum   25ileum 25jejunum   29cecum 
    ##      3451      3451      3451      3451      3451      3451      3451      3451 
    ##   29ileum 29jejunum   33cecum   33ileum 33jejunum   37cecum   37ileum 37jejunum 
    ##      3451      3451      3451      3451      3451      3451      3451      3451 
    ##   41cecum   41ileum 41jejunum   45cecum   45ileum 45jejunum   49cecum   49ileum 
    ##      3451      3451      3451      3451      3451      3451      3451      3451 
    ## 49jejunum   53cecum   53ileum 53jejunum   57cecum   57ileum 57jejunum    5cecum 
    ##      3451      3451      3451      3451      3451      3451      3451      3451 
    ##    5ileum  5jejunum   61cecum   61ileum 61jejunum    9cecum    9ileum  9jejunum 
    ##      3451      3451      3451      3451      3451      3451      3451      3451

``` r
dat <- (dat[2:9] / rowSums(dat[2:9])) * 100
rowSums(dat) # should all be equal to 100
```

    ##   13cecum   13ileum 13jejunum   17cecum   17ileum 17jejunum    1cecum    1ileum 
    ##       100       100       100       100       100       100       100       100 
    ##  1jejunum   21cecum   21ileum 21jejunum   25cecum   25ileum 25jejunum   29cecum 
    ##       100       100       100       100       100       100       100       100 
    ##   29ileum 29jejunum   33cecum   33ileum 33jejunum   37cecum   37ileum 37jejunum 
    ##       100       100       100       100       100       100       100       100 
    ##   41cecum   41ileum 41jejunum   45cecum   45ileum 45jejunum   49cecum   49ileum 
    ##       100       100       100       100       100       100       100       100 
    ## 49jejunum   53cecum   53ileum 53jejunum   57cecum   57ileum 57jejunum    5cecum 
    ##       100       100       100       100       100       100       100       100 
    ##    5ileum  5jejunum   61cecum   61ileum 61jejunum    9cecum    9ileum  9jejunum 
    ##       100       100       100       100       100       100       100       100

Calculate % of total gd T-IELs:

``` r
colnames(dat)
```

    ## [1] "AB_CD8AnegCD8Bneg" "AB_CD8AnegCD8Bpos" "AB_CD8AposCD8Bneg"
    ## [4] "AB_CD8AposCD8Bpos" "GD_CD2negCD8Aneg"  "GD_CD2negCD8Apos" 
    ## [7] "GD_CD2posCD8Aneg"  "GD_CD2posCD8Apos"

``` r
dat$GD <- rowSums(dat[5:8])
```

Calculate % of each subset within total gd T-IELs:

``` r
colnames(dat)
```

    ## [1] "AB_CD8AnegCD8Bneg" "AB_CD8AnegCD8Bpos" "AB_CD8AposCD8Bneg"
    ## [4] "AB_CD8AposCD8Bpos" "GD_CD2negCD8Aneg"  "GD_CD2negCD8Apos" 
    ## [7] "GD_CD2posCD8Aneg"  "GD_CD2posCD8Apos"  "GD"

``` r
dat$GD_CD2negCD8Aneg_GD <- dat$GD_CD2negCD8Aneg/dat$GD *100
dat$GD_CD2posCD8Aneg_GD <- dat$GD_CD2posCD8Aneg/dat$GD *100
dat$GD_CD2negCD8Apos_GD <- dat$GD_CD2negCD8Apos/dat$GD *100
dat$GD_CD2posCD8Apos_GD <- dat$GD_CD2posCD8Apos/dat$GD *100
```

Calculate % of each subset within total ab T-IELs:

``` r
colnames(dat)
```

    ##  [1] "AB_CD8AnegCD8Bneg"   "AB_CD8AnegCD8Bpos"   "AB_CD8AposCD8Bneg"  
    ##  [4] "AB_CD8AposCD8Bpos"   "GD_CD2negCD8Aneg"    "GD_CD2negCD8Apos"   
    ##  [7] "GD_CD2posCD8Aneg"    "GD_CD2posCD8Apos"    "GD"                 
    ## [10] "GD_CD2negCD8Aneg_GD" "GD_CD2posCD8Aneg_GD" "GD_CD2negCD8Apos_GD"
    ## [13] "GD_CD2posCD8Apos_GD"

``` r
dat$AB_CD8AnegCD8Bneg_AB <- dat$AB_CD8AnegCD8Bneg/(100-dat$GD) *100
dat$AB_CD8AposCD8Bneg_AB <- dat$AB_CD8AposCD8Bneg/(100-dat$GD) *100
dat$AB_CD8AnegCD8Bpos_AB <- dat$AB_CD8AnegCD8Bpos/(100-dat$GD) *100
dat$AB_CD8AposCD8Bpos_AB <- dat$AB_CD8AposCD8Bpos/(100-dat$GD) *100
```

Merge counts and meta data:

``` r
dat$SampleID <- rownames(dat)
dat <- do.add.cols(dat = dat, base.col = "SampleID", add.dat = meta.dat, add.by = "Var2")
```

    ## Loading required package: data.table

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## Step 1/3. Mapping data

    ## Step 2/3. Merging data

    ## Step 3/3. Returning data

``` r
head(dat)
```

    ##    AB_CD8AnegCD8Bneg AB_CD8AnegCD8Bpos AB_CD8AposCD8Bneg AB_CD8AposCD8Bpos
    ## 1:        0.08693132        0.05795422          7.215300          66.93712
    ## 2:        1.39090119        0.02897711          5.070994          56.85309
    ## 3:        3.94088670        0.11590843          2.897711          70.50130
    ## 4:        0.20283976        0.00000000          7.476094          67.74848
    ## 5:        1.27499276        0.02897711          4.926108          67.98030
    ## 6:        1.97044335        0.02897711          3.071573          76.00695
    ##    GD_CD2negCD8Aneg GD_CD2negCD8Apos GD_CD2posCD8Aneg GD_CD2posCD8Apos       GD
    ## 1:       0.05795422        1.0721530       0.05795422         24.51463 25.70269
    ## 2:       0.02897711        1.3619241       5.04201681         30.22312 36.65604
    ## 3:       0.08693132        1.4488554      10.92436975         10.08403 22.54419
    ## 4:       0.11590843        1.6516952       0.60851927         22.19646 24.57259
    ## 5:       0.20283976        0.7534048       4.28861200         20.54477 25.78963
    ## 6:       0.00000000        0.3767024       4.66531440         13.88003 18.92205
    ##    GD_CD2negCD8Aneg_GD GD_CD2posCD8Aneg_GD GD_CD2negCD8Apos_GD
    ## 1:          0.22547914           0.2254791            4.171364
    ## 2:          0.07905138          13.7549407            3.715415
    ## 3:          0.38560411          48.4575835            6.426735
    ## 4:          0.47169811           2.4764151            6.721698
    ## 5:          0.78651685          16.6292135            2.921348
    ## 6:          0.00000000          24.6554364            1.990812
    ##    GD_CD2posCD8Apos_GD AB_CD8AnegCD8Bneg_AB AB_CD8AposCD8Bneg_AB
    ## 1:            95.37768            0.1170047             9.711388
    ## 2:            82.45059            2.1957914             8.005489
    ## 3:            44.73008            5.0879162             3.741115
    ## 4:            90.33019            0.2689205             9.911640
    ## 5:            79.66292            1.7180789             6.638032
    ## 6:            73.35375            2.4303074             3.788420
    ##    AB_CD8AnegCD8Bpos_AB AB_CD8AposCD8Bpos_AB  SampleID
    ## 1:           0.07800312             90.09360   13cecum
    ## 2:           0.04574565             89.75297   13ileum
    ## 3:           0.14964459             91.02132 13jejunum
    ## 4:           0.00000000             89.81944   17cecum
    ## 5:           0.03904725             91.60484   17ileum
    ## 6:           0.03573981             93.74553 17jejunum
    ##                                                         FileName Batch
    ## 1:  13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209.csv     B
    ## 2:  13il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209.csv     B
    ## 3: 13jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209.csv     B
    ## 4:  17ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209.csv     C
    ## 5:  17il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209.csv     C
    ## 6: 17jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209.csv     C
    ##    WeeksOfAge  Tissue AnimalID
    ## 1:          5   cecum       13
    ## 2:          5   ileum       13
    ## 3:          5 jejunum       13
    ## 4:          5   cecum       17
    ## 5:          5   ileum       17
    ## 6:          5 jejunum       17

Make directory to store results:

``` r
setwd('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/')
dir.create('StatResults')
setwd('StatResults')
```

Compare population frequencies within total T-IELs

Across tissues:

``` r
colnames(dat)
```

    ##  [1] "AB_CD8AnegCD8Bneg"    "AB_CD8AnegCD8Bpos"    "AB_CD8AposCD8Bneg"   
    ##  [4] "AB_CD8AposCD8Bpos"    "GD_CD2negCD8Aneg"     "GD_CD2negCD8Apos"    
    ##  [7] "GD_CD2posCD8Aneg"     "GD_CD2posCD8Apos"     "GD"                  
    ## [10] "GD_CD2negCD8Aneg_GD"  "GD_CD2posCD8Aneg_GD"  "GD_CD2negCD8Apos_GD" 
    ## [13] "GD_CD2posCD8Apos_GD"  "AB_CD8AnegCD8Bneg_AB" "AB_CD8AposCD8Bneg_AB"
    ## [16] "AB_CD8AnegCD8Bpos_AB" "AB_CD8AposCD8Bpos_AB" "SampleID"            
    ## [19] "FileName"             "Batch"                "WeeksOfAge"          
    ## [22] "Tissue"               "AnimalID"

``` r
for (i in 1:9) {
  formula = as.formula( paste(colnames(dat[,..i]), 'Tissue', sep = '~'))
 stat <- dat %>%
        group_by(WeeksOfAge) %>%
        wilcox_test(formula = formula, paired = TRUE) %>% # use paired since tissues derived from same animal
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentTotalTIELs_PairedTissueWilcoxon_',
                                colnames(dat[,..i])), '.xlsx'))
}
```

Across ages:

``` r
colnames(dat)
```

    ##  [1] "AB_CD8AnegCD8Bneg"    "AB_CD8AnegCD8Bpos"    "AB_CD8AposCD8Bneg"   
    ##  [4] "AB_CD8AposCD8Bpos"    "GD_CD2negCD8Aneg"     "GD_CD2negCD8Apos"    
    ##  [7] "GD_CD2posCD8Aneg"     "GD_CD2posCD8Apos"     "GD"                  
    ## [10] "GD_CD2negCD8Aneg_GD"  "GD_CD2posCD8Aneg_GD"  "GD_CD2negCD8Apos_GD" 
    ## [13] "GD_CD2posCD8Apos_GD"  "AB_CD8AnegCD8Bneg_AB" "AB_CD8AposCD8Bneg_AB"
    ## [16] "AB_CD8AnegCD8Bpos_AB" "AB_CD8AposCD8Bpos_AB" "SampleID"            
    ## [19] "FileName"             "Batch"                "WeeksOfAge"          
    ## [22] "Tissue"               "AnimalID"

``` r
for (i in 1:9) {
formula = as.formula( paste(colnames(dat[,..i]), 'WeeksOfAge', sep = '~'))
 stat <- dat %>%
        group_by(Tissue) %>%
        wilcox_test(formula = formula, paired = FALSE) %>%
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentTotalTIELs_UnpairedAgeWilcoxon_', colnames(dat[,..i])), '.xlsx'))
}
```

Compare population frequencies within gd T-IELs:

Across tissues:

``` r
colnames(dat)
```

    ##  [1] "AB_CD8AnegCD8Bneg"    "AB_CD8AnegCD8Bpos"    "AB_CD8AposCD8Bneg"   
    ##  [4] "AB_CD8AposCD8Bpos"    "GD_CD2negCD8Aneg"     "GD_CD2negCD8Apos"    
    ##  [7] "GD_CD2posCD8Aneg"     "GD_CD2posCD8Apos"     "GD"                  
    ## [10] "GD_CD2negCD8Aneg_GD"  "GD_CD2posCD8Aneg_GD"  "GD_CD2negCD8Apos_GD" 
    ## [13] "GD_CD2posCD8Apos_GD"  "AB_CD8AnegCD8Bneg_AB" "AB_CD8AposCD8Bneg_AB"
    ## [16] "AB_CD8AnegCD8Bpos_AB" "AB_CD8AposCD8Bpos_AB" "SampleID"            
    ## [19] "FileName"             "Batch"                "WeeksOfAge"          
    ## [22] "Tissue"               "AnimalID"

``` r
for (i in 10:13) {
  formula = as.formula( paste(colnames(dat[,..i]), 'Tissue', sep = '~'))
 stat <- dat %>%
        group_by(WeeksOfAge) %>%
        wilcox_test(formula = formula, paired = TRUE) %>% # use paired since tissues derived from same animal
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentGDTIELs_PairedTissueWilcoxon_',
                                colnames(dat[,..i])), '.xlsx'))
}
```

Across ages:

``` r
colnames(dat)
```

    ##  [1] "AB_CD8AnegCD8Bneg"    "AB_CD8AnegCD8Bpos"    "AB_CD8AposCD8Bneg"   
    ##  [4] "AB_CD8AposCD8Bpos"    "GD_CD2negCD8Aneg"     "GD_CD2negCD8Apos"    
    ##  [7] "GD_CD2posCD8Aneg"     "GD_CD2posCD8Apos"     "GD"                  
    ## [10] "GD_CD2negCD8Aneg_GD"  "GD_CD2posCD8Aneg_GD"  "GD_CD2negCD8Apos_GD" 
    ## [13] "GD_CD2posCD8Apos_GD"  "AB_CD8AnegCD8Bneg_AB" "AB_CD8AposCD8Bneg_AB"
    ## [16] "AB_CD8AnegCD8Bpos_AB" "AB_CD8AposCD8Bpos_AB" "SampleID"            
    ## [19] "FileName"             "Batch"                "WeeksOfAge"          
    ## [22] "Tissue"               "AnimalID"

``` r
for (i in 10:13) {
formula = as.formula( paste(colnames(dat[,..i]), 'WeeksOfAge', sep = '~'))
 stat <- dat %>%
        group_by(Tissue) %>%
        wilcox_test(formula = formula, paired = FALSE) %>%
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentGDTIELs_UnpairedAgeWilcoxon_', colnames(dat[,..i])), '.xlsx'))
}
```

Compare population frequencies within ab T-IELs:

Across tissues:

``` r
colnames(dat)
```

    ##  [1] "AB_CD8AnegCD8Bneg"    "AB_CD8AnegCD8Bpos"    "AB_CD8AposCD8Bneg"   
    ##  [4] "AB_CD8AposCD8Bpos"    "GD_CD2negCD8Aneg"     "GD_CD2negCD8Apos"    
    ##  [7] "GD_CD2posCD8Aneg"     "GD_CD2posCD8Apos"     "GD"                  
    ## [10] "GD_CD2negCD8Aneg_GD"  "GD_CD2posCD8Aneg_GD"  "GD_CD2negCD8Apos_GD" 
    ## [13] "GD_CD2posCD8Apos_GD"  "AB_CD8AnegCD8Bneg_AB" "AB_CD8AposCD8Bneg_AB"
    ## [16] "AB_CD8AnegCD8Bpos_AB" "AB_CD8AposCD8Bpos_AB" "SampleID"            
    ## [19] "FileName"             "Batch"                "WeeksOfAge"          
    ## [22] "Tissue"               "AnimalID"

``` r
for (i in 14:17) {
  formula = as.formula( paste(colnames(dat[,..i]), 'Tissue', sep = '~'))
 stat <- dat %>%
        group_by(WeeksOfAge) %>%
        wilcox_test(formula = formula, paired = TRUE) %>% # use paired since tissues derived from same animal
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentABTIELs_PairedTissueWilcoxon_',
                                colnames(dat[,..i])), '.xlsx'))
}
```

Across ages:

``` r
colnames(dat)
```

    ##  [1] "AB_CD8AnegCD8Bneg"    "AB_CD8AnegCD8Bpos"    "AB_CD8AposCD8Bneg"   
    ##  [4] "AB_CD8AposCD8Bpos"    "GD_CD2negCD8Aneg"     "GD_CD2negCD8Apos"    
    ##  [7] "GD_CD2posCD8Aneg"     "GD_CD2posCD8Apos"     "GD"                  
    ## [10] "GD_CD2negCD8Aneg_GD"  "GD_CD2posCD8Aneg_GD"  "GD_CD2negCD8Apos_GD" 
    ## [13] "GD_CD2posCD8Apos_GD"  "AB_CD8AnegCD8Bneg_AB" "AB_CD8AposCD8Bneg_AB"
    ## [16] "AB_CD8AnegCD8Bpos_AB" "AB_CD8AposCD8Bpos_AB" "SampleID"            
    ## [19] "FileName"             "Batch"                "WeeksOfAge"          
    ## [22] "Tissue"               "AnimalID"

``` r
for (i in 14:17) {
formula = as.formula( paste(colnames(dat[,..i]), 'WeeksOfAge', sep = '~'))
 stat <- dat %>%
        group_by(Tissue) %>%
        wilcox_test(formula = formula, paired = FALSE) %>%
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentABTIELs_UnpairedAgeWilcoxon_', colnames(dat[,..i])), '.xlsx'))
}
```

### View session information

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
    ## [1] data.table_1.14.0 writexl_1.4.0     dplyr_1.0.7       Spectre_0.4.1    
    ## [5] rstatix_0.7.0     readxl_1.3.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] zip_2.2.0        Rcpp_1.0.7       cellranger_1.1.0 pillar_1.6.2    
    ##  [5] compiler_4.1.3   forcats_0.5.1    tools_4.1.3      digest_0.6.27   
    ##  [9] evaluate_0.14    lifecycle_1.0.0  tibble_3.1.4     pkgconfig_2.0.3 
    ## [13] rlang_0.4.11     openxlsx_4.2.4   DBI_1.1.1        curl_4.3.2      
    ## [17] yaml_2.2.1       haven_2.4.3      xfun_0.26        rio_0.5.27      
    ## [21] fastmap_1.1.0    stringr_1.4.0    knitr_1.34       hms_1.1.0       
    ## [25] generics_0.1.0   vctrs_0.3.8      tidyselect_1.1.1 glue_1.4.2      
    ## [29] R6_2.5.1         fansi_0.5.0      foreign_0.8-82   rmarkdown_2.11  
    ## [33] carData_3.0-4    purrr_0.3.4      tidyr_1.1.3      car_3.0-11      
    ## [37] magrittr_2.0.1   backports_1.2.1  ellipsis_0.3.2   htmltools_0.5.2 
    ## [41] abind_1.4-5      assertthat_0.2.1 utf8_1.2.2       stringi_1.7.4   
    ## [45] broom_0.7.9      crayon_1.4.1
