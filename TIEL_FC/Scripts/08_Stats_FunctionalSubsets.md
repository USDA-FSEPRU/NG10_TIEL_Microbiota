Statistical analysis of T-IEL functional subsets
================
Jayne Wiarda
2022May04

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

We further classified cells within eight major T-IEL populations as:

-   CD16+CD27+
-   CD16+CD27-
-   CD16-CD27-
-   CD16-CD27+

The objective of the below analyses is to further investigate abundances
and statistical comparisons for different T-IEL populations across age,
tissue, and specific T-IEL populations.

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
dat <- as.data.frame(read_excel('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/CellGatingAnnotation/TotalCellCountsBySample_32pops.xlsx'))
rownames(dat) <- dat$Var2
```

Read in meta data & organize further:

``` r
meta.dat <- read_excel("/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/NG10_MetaData_FC.xlsx") # read in meta data
meta.dat$Var2 <- paste0(meta.dat$AnimalID, meta.dat$Tissue)
```

Convert counts data to percentages of total T-IELs:

Calculate as %CD16/CD27 subsets within parent T-IEL population

``` r
colnames(dat)
```

    ##  [1] "Var2"                             "AB_CD8AnegCD8Bneg_CD16negCD27neg"
    ##  [3] "AB_CD8AnegCD8Bneg_CD16negCD27pos" "AB_CD8AnegCD8Bneg_CD16posCD27neg"
    ##  [5] "AB_CD8AnegCD8Bneg_CD16posCD27pos" "AB_CD8AnegCD8Bpos_CD16negCD27neg"
    ##  [7] "AB_CD8AnegCD8Bpos_CD16negCD27pos" "AB_CD8AnegCD8Bpos_CD16posCD27neg"
    ##  [9] "AB_CD8AnegCD8Bpos_CD16posCD27pos" "AB_CD8AposCD8Bneg_CD16negCD27neg"
    ## [11] "AB_CD8AposCD8Bneg_CD16negCD27pos" "AB_CD8AposCD8Bneg_CD16posCD27neg"
    ## [13] "AB_CD8AposCD8Bneg_CD16posCD27pos" "AB_CD8AposCD8Bpos_CD16negCD27neg"
    ## [15] "AB_CD8AposCD8Bpos_CD16negCD27pos" "AB_CD8AposCD8Bpos_CD16posCD27neg"
    ## [17] "AB_CD8AposCD8Bpos_CD16posCD27pos" "GD_CD2negCD8Aneg_CD16negCD27neg" 
    ## [19] "GD_CD2negCD8Aneg_CD16negCD27pos"  "GD_CD2negCD8Aneg_CD16posCD27neg" 
    ## [21] "GD_CD2negCD8Aneg_CD16posCD27pos"  "GD_CD2negCD8Apos_CD16negCD27neg" 
    ## [23] "GD_CD2negCD8Apos_CD16negCD27pos"  "GD_CD2negCD8Apos_CD16posCD27neg" 
    ## [25] "GD_CD2negCD8Apos_CD16posCD27pos"  "GD_CD2posCD8Aneg_CD16negCD27neg" 
    ## [27] "GD_CD2posCD8Aneg_CD16negCD27pos"  "GD_CD2posCD8Aneg_CD16posCD27neg" 
    ## [29] "GD_CD2posCD8Aneg_CD16posCD27pos"  "GD_CD2posCD8Apos_CD16negCD27neg" 
    ## [31] "GD_CD2posCD8Apos_CD16negCD27pos"  "GD_CD2posCD8Apos_CD16posCD27neg" 
    ## [33] "GD_CD2posCD8Apos_CD16posCD27pos"

``` r
dat <- dat %>% select(1, 10:17, 26:33) # select only counts data for four major T-IEL primary populations
dat <- do.add.cols(dat = dat, base.col = "Var2", add.dat = meta.dat, add.by = "Var2") # add meta data
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
dat$SampleID <- dat$Var2
colnames(dat)
```

    ##  [1] "Var2"                             "AB_CD8AposCD8Bneg_CD16negCD27neg"
    ##  [3] "AB_CD8AposCD8Bneg_CD16negCD27pos" "AB_CD8AposCD8Bneg_CD16posCD27neg"
    ##  [5] "AB_CD8AposCD8Bneg_CD16posCD27pos" "AB_CD8AposCD8Bpos_CD16negCD27neg"
    ##  [7] "AB_CD8AposCD8Bpos_CD16negCD27pos" "AB_CD8AposCD8Bpos_CD16posCD27neg"
    ##  [9] "AB_CD8AposCD8Bpos_CD16posCD27pos" "GD_CD2posCD8Aneg_CD16negCD27neg" 
    ## [11] "GD_CD2posCD8Aneg_CD16negCD27pos"  "GD_CD2posCD8Aneg_CD16posCD27neg" 
    ## [13] "GD_CD2posCD8Aneg_CD16posCD27pos"  "GD_CD2posCD8Apos_CD16negCD27neg" 
    ## [15] "GD_CD2posCD8Apos_CD16negCD27pos"  "GD_CD2posCD8Apos_CD16posCD27neg" 
    ## [17] "GD_CD2posCD8Apos_CD16posCD27pos"  "FileName"                        
    ## [19] "Batch"                            "WeeksOfAge"                      
    ## [21] "Tissue"                           "AnimalID"                        
    ## [23] "SampleID"

``` r
subdat1 <- dat %>% select(2:5) 
subdat1 <- (subdat1 / rowSums(subdat1)) * 100 # calculate percentages of each CD16/CD27 subset within total population
rowSums(subdat1) # should all be equal to 100
```

    ##  [1] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
    ## [20] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
    ## [39] 100 100 100 100 100 100 100 100 100 100

``` r
subdat2 <- dat %>% select(6:9) 
subdat2 <- (subdat2 / rowSums(subdat2)) * 100 # calculate percentages of each CD16/CD27 subset within total population
rowSums(subdat2) # should all be equal to 100
```

    ##  [1] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
    ## [20] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
    ## [39] 100 100 100 100 100 100 100 100 100 100

``` r
subdat3 <- dat %>% select(10:13) 
subdat3 <- (subdat3 / rowSums(subdat3)) * 100 # calculate percentages of each CD16/CD27 subset within total population
rowSums(subdat3) # should all be equal to 100
```

    ##  [1] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
    ## [20] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
    ## [39] 100 100 100 100 100 100 100 100 100 100

``` r
subdat4 <- dat %>% select(14:17) 
subdat4 <- (subdat4 / rowSums(subdat4)) * 100 # calculate percentages of each CD16/CD27 subset within total population
rowSums(subdat4) # should all be equal to 100
```

    ##  [1] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
    ## [20] 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
    ## [39] 100 100 100 100 100 100 100 100 100 100

``` r
colnames(dat)
```

    ##  [1] "Var2"                             "AB_CD8AposCD8Bneg_CD16negCD27neg"
    ##  [3] "AB_CD8AposCD8Bneg_CD16negCD27pos" "AB_CD8AposCD8Bneg_CD16posCD27neg"
    ##  [5] "AB_CD8AposCD8Bneg_CD16posCD27pos" "AB_CD8AposCD8Bpos_CD16negCD27neg"
    ##  [7] "AB_CD8AposCD8Bpos_CD16negCD27pos" "AB_CD8AposCD8Bpos_CD16posCD27neg"
    ##  [9] "AB_CD8AposCD8Bpos_CD16posCD27pos" "GD_CD2posCD8Aneg_CD16negCD27neg" 
    ## [11] "GD_CD2posCD8Aneg_CD16negCD27pos"  "GD_CD2posCD8Aneg_CD16posCD27neg" 
    ## [13] "GD_CD2posCD8Aneg_CD16posCD27pos"  "GD_CD2posCD8Apos_CD16negCD27neg" 
    ## [15] "GD_CD2posCD8Apos_CD16negCD27pos"  "GD_CD2posCD8Apos_CD16posCD27neg" 
    ## [17] "GD_CD2posCD8Apos_CD16posCD27pos"  "FileName"                        
    ## [19] "Batch"                            "WeeksOfAge"                      
    ## [21] "Tissue"                           "AnimalID"                        
    ## [23] "SampleID"

``` r
dat <- cbind(dat[,20:23], subdat1, subdat2, subdat3, subdat4) # combine T-IEL percentages back with meta data
head(dat)
```

    ##    WeeksOfAge  Tissue AnimalID  SampleID AB_CD8AposCD8Bneg_CD16negCD27neg
    ## 1:          5   cecum       13   13cecum                         20.08032
    ## 2:          5   ileum       13   13ileum                         25.14286
    ## 3:          5 jejunum       13 13jejunum                         27.00000
    ## 4:          5   cecum       17   17cecum                         23.64341
    ## 5:          5   ileum       17   17ileum                         17.05882
    ## 6:          5 jejunum       17 17jejunum                         42.45283
    ##    AB_CD8AposCD8Bneg_CD16negCD27pos AB_CD8AposCD8Bneg_CD16posCD27neg
    ## 1:                         7.228916                         55.82329
    ## 2:                        18.285714                         32.57143
    ## 3:                        32.000000                         23.00000
    ## 4:                         7.364341                         58.91473
    ## 5:                        18.823529                         32.94118
    ## 6:                        16.037736                         26.41509
    ##    AB_CD8AposCD8Bneg_CD16posCD27pos AB_CD8AposCD8Bpos_CD16negCD27neg
    ## 1:                         16.86747                         35.10823
    ## 2:                         24.00000                         27.01325
    ## 3:                         18.00000                         10.97411
    ## 4:                         10.07752                         49.61506
    ## 5:                         31.17647                         26.00171
    ## 6:                         15.09434                         25.84827
    ##    AB_CD8AposCD8Bpos_CD16negCD27pos AB_CD8AposCD8Bpos_CD16posCD27neg
    ## 1:                         21.25541                        35.194805
    ## 2:                         44.95413                        17.278287
    ## 3:                         83.76490                         1.972873
    ## 4:                         14.41403                        31.822070
    ## 5:                         46.67519                        15.899403
    ## 6:                         67.21311                         3.240564
    ##    AB_CD8AposCD8Bpos_CD16posCD27pos GD_CD2posCD8Aneg_CD16negCD27neg
    ## 1:                         8.441558                        50.00000
    ## 2:                        10.754332                        62.64368
    ## 3:                         3.288122                        88.85942
    ## 4:                         4.148845                        52.38095
    ## 5:                        11.423700                        34.45946
    ## 6:                         3.698056                        65.21739
    ##    GD_CD2posCD8Aneg_CD16negCD27pos GD_CD2posCD8Aneg_CD16posCD27neg
    ## 1:                        0.000000                       50.000000
    ## 2:                        5.172414                       21.264368
    ## 3:                        6.631300                        3.713528
    ## 4:                       19.047619                       28.571429
    ## 5:                       20.270270                       25.675676
    ## 6:                       11.801242                       13.664596
    ##    GD_CD2posCD8Aneg_CD16posCD27pos GD_CD2posCD8Apos_CD16negCD27neg
    ## 1:                        0.000000                       19.976359
    ## 2:                       10.919540                       15.436242
    ## 3:                        0.795756                       18.390805
    ## 4:                        0.000000                       32.898172
    ## 5:                       19.594595                        5.500705
    ## 6:                        9.316770                       18.580376
    ##    GD_CD2posCD8Apos_CD16negCD27pos GD_CD2posCD8Apos_CD16posCD27neg
    ## 1:                        8.865248                        52.95508
    ## 2:                       31.160115                        26.74976
    ## 3:                       37.068966                        15.51724
    ## 4:                       13.446475                        40.86162
    ## 5:                       36.671368                        20.73343
    ## 6:                       35.699374                        17.32777
    ##    GD_CD2posCD8Apos_CD16posCD27pos
    ## 1:                        18.20331
    ## 2:                        26.65388
    ## 3:                        29.02299
    ## 4:                        12.79373
    ## 5:                        37.09450
    ## 6:                        28.39248

``` r
dat$AB_CD8AposCD8Bneg_CD16pos <- dat$AB_CD8AposCD8Bneg_CD16posCD27neg + dat$AB_CD8AposCD8Bneg_CD16posCD27pos
dat$AB_CD8AposCD8Bneg_CD27 <- dat$AB_CD8AposCD8Bneg_CD16negCD27pos + dat$AB_CD8AposCD8Bneg_CD16posCD27pos

dat$AB_CD8AposCD8Bpos_CD16pos <- dat$AB_CD8AposCD8Bpos_CD16posCD27neg + dat$AB_CD8AposCD8Bpos_CD16posCD27pos
dat$AB_CD8AposCD8Bpos_CD27 <- dat$AB_CD8AposCD8Bpos_CD16negCD27pos + dat$AB_CD8AposCD8Bpos_CD16posCD27pos

dat$GD_CD2posCD8Apos_CD16pos <- dat$GD_CD2posCD8Apos_CD16posCD27neg + dat$GD_CD2posCD8Apos_CD16posCD27pos
dat$GD_CD2posCD8Apos_CD27 <- dat$GD_CD2posCD8Apos_CD16negCD27pos + dat$GD_CD2posCD8Apos_CD16posCD27pos

dat$GD_CD2posCD8Aneg_CD16pos <- dat$GD_CD2posCD8Aneg_CD16posCD27neg + dat$GD_CD2posCD8Aneg_CD16posCD27pos
dat$GD_CD2posCD8Aneg_CD27 <- dat$GD_CD2posCD8Aneg_CD16negCD27pos + dat$GD_CD2posCD8Aneg_CD16posCD27pos
```

Make directory to store results:

THis directory already exists if you ran script 06\_Stats\_8pops

``` r
setwd('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/')
dir.create('StatResults')
setwd('StatResults')
```

Compare population frequencies within each parent population:

Across tissues:

``` r
colnames(dat)
```

    ##  [1] "WeeksOfAge"                       "Tissue"                          
    ##  [3] "AnimalID"                         "SampleID"                        
    ##  [5] "AB_CD8AposCD8Bneg_CD16negCD27neg" "AB_CD8AposCD8Bneg_CD16negCD27pos"
    ##  [7] "AB_CD8AposCD8Bneg_CD16posCD27neg" "AB_CD8AposCD8Bneg_CD16posCD27pos"
    ##  [9] "AB_CD8AposCD8Bpos_CD16negCD27neg" "AB_CD8AposCD8Bpos_CD16negCD27pos"
    ## [11] "AB_CD8AposCD8Bpos_CD16posCD27neg" "AB_CD8AposCD8Bpos_CD16posCD27pos"
    ## [13] "GD_CD2posCD8Aneg_CD16negCD27neg"  "GD_CD2posCD8Aneg_CD16negCD27pos" 
    ## [15] "GD_CD2posCD8Aneg_CD16posCD27neg"  "GD_CD2posCD8Aneg_CD16posCD27pos" 
    ## [17] "GD_CD2posCD8Apos_CD16negCD27neg"  "GD_CD2posCD8Apos_CD16negCD27pos" 
    ## [19] "GD_CD2posCD8Apos_CD16posCD27neg"  "GD_CD2posCD8Apos_CD16posCD27pos" 
    ## [21] "AB_CD8AposCD8Bneg_CD16pos"        "AB_CD8AposCD8Bneg_CD27"          
    ## [23] "AB_CD8AposCD8Bpos_CD16pos"        "AB_CD8AposCD8Bpos_CD27"          
    ## [25] "GD_CD2posCD8Apos_CD16pos"         "GD_CD2posCD8Apos_CD27"           
    ## [27] "GD_CD2posCD8Aneg_CD16pos"         "GD_CD2posCD8Aneg_CD27"

``` r
for (i in 5:28) {
  formula = as.formula( paste(colnames(dat[,..i]), 'Tissue', sep = '~'))
 stat <- dat %>%
        group_by(WeeksOfAge) %>%
        wilcox_test(formula = formula, paired = TRUE) %>% # use paired since tissues derived from same animal
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentParentPopulation_PairedTissueWilcoxon_',
                                colnames(dat[,..i])), '.xlsx'))
}
```

Across ages:

``` r
colnames(dat)
```

    ##  [1] "WeeksOfAge"                       "Tissue"                          
    ##  [3] "AnimalID"                         "SampleID"                        
    ##  [5] "AB_CD8AposCD8Bneg_CD16negCD27neg" "AB_CD8AposCD8Bneg_CD16negCD27pos"
    ##  [7] "AB_CD8AposCD8Bneg_CD16posCD27neg" "AB_CD8AposCD8Bneg_CD16posCD27pos"
    ##  [9] "AB_CD8AposCD8Bpos_CD16negCD27neg" "AB_CD8AposCD8Bpos_CD16negCD27pos"
    ## [11] "AB_CD8AposCD8Bpos_CD16posCD27neg" "AB_CD8AposCD8Bpos_CD16posCD27pos"
    ## [13] "GD_CD2posCD8Aneg_CD16negCD27neg"  "GD_CD2posCD8Aneg_CD16negCD27pos" 
    ## [15] "GD_CD2posCD8Aneg_CD16posCD27neg"  "GD_CD2posCD8Aneg_CD16posCD27pos" 
    ## [17] "GD_CD2posCD8Apos_CD16negCD27neg"  "GD_CD2posCD8Apos_CD16negCD27pos" 
    ## [19] "GD_CD2posCD8Apos_CD16posCD27neg"  "GD_CD2posCD8Apos_CD16posCD27pos" 
    ## [21] "AB_CD8AposCD8Bneg_CD16pos"        "AB_CD8AposCD8Bneg_CD27"          
    ## [23] "AB_CD8AposCD8Bpos_CD16pos"        "AB_CD8AposCD8Bpos_CD27"          
    ## [25] "GD_CD2posCD8Apos_CD16pos"         "GD_CD2posCD8Apos_CD27"           
    ## [27] "GD_CD2posCD8Aneg_CD16pos"         "GD_CD2posCD8Aneg_CD27"

``` r
for (i in 5:28) {
formula = as.formula( paste(colnames(dat[,..i]), 'WeeksOfAge', sep = '~'))
 stat <- dat %>%
        group_by(Tissue) %>%
        wilcox_test(formula = formula, paired = FALSE) %>%
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentParentPopulation_UnpairedAgeWilcoxon_', colnames(dat[,..i])), '.xlsx'))
}
```

Across primary T-IEL populations:

``` r
colnames(dat)
```

    ##  [1] "WeeksOfAge"                       "Tissue"                          
    ##  [3] "AnimalID"                         "SampleID"                        
    ##  [5] "AB_CD8AposCD8Bneg_CD16negCD27neg" "AB_CD8AposCD8Bneg_CD16negCD27pos"
    ##  [7] "AB_CD8AposCD8Bneg_CD16posCD27neg" "AB_CD8AposCD8Bneg_CD16posCD27pos"
    ##  [9] "AB_CD8AposCD8Bpos_CD16negCD27neg" "AB_CD8AposCD8Bpos_CD16negCD27pos"
    ## [11] "AB_CD8AposCD8Bpos_CD16posCD27neg" "AB_CD8AposCD8Bpos_CD16posCD27pos"
    ## [13] "GD_CD2posCD8Aneg_CD16negCD27neg"  "GD_CD2posCD8Aneg_CD16negCD27pos" 
    ## [15] "GD_CD2posCD8Aneg_CD16posCD27neg"  "GD_CD2posCD8Aneg_CD16posCD27pos" 
    ## [17] "GD_CD2posCD8Apos_CD16negCD27neg"  "GD_CD2posCD8Apos_CD16negCD27pos" 
    ## [19] "GD_CD2posCD8Apos_CD16posCD27neg"  "GD_CD2posCD8Apos_CD16posCD27pos" 
    ## [21] "AB_CD8AposCD8Bneg_CD16pos"        "AB_CD8AposCD8Bneg_CD27"          
    ## [23] "AB_CD8AposCD8Bpos_CD16pos"        "AB_CD8AposCD8Bpos_CD27"          
    ## [25] "GD_CD2posCD8Apos_CD16pos"         "GD_CD2posCD8Apos_CD27"           
    ## [27] "GD_CD2posCD8Aneg_CD16pos"         "GD_CD2posCD8Aneg_CD27"

``` r
dat <- dat %>% gather(key = "subset", value = "Percentages", 5:28)
dat$celltype <- sub("^([^_]*_[^_]*).*", "\\1", dat$subset)
dat$cellsubset <- sub("^.*\\_","",dat$subset)
head(dat)
```

    ##   WeeksOfAge  Tissue AnimalID  SampleID                           subset
    ## 1          5   cecum       13   13cecum AB_CD8AposCD8Bneg_CD16negCD27neg
    ## 2          5   ileum       13   13ileum AB_CD8AposCD8Bneg_CD16negCD27neg
    ## 3          5 jejunum       13 13jejunum AB_CD8AposCD8Bneg_CD16negCD27neg
    ## 4          5   cecum       17   17cecum AB_CD8AposCD8Bneg_CD16negCD27neg
    ## 5          5   ileum       17   17ileum AB_CD8AposCD8Bneg_CD16negCD27neg
    ## 6          5 jejunum       17 17jejunum AB_CD8AposCD8Bneg_CD16negCD27neg
    ##   Percentages          celltype     cellsubset
    ## 1    20.08032 AB_CD8AposCD8Bneg CD16negCD27neg
    ## 2    25.14286 AB_CD8AposCD8Bneg CD16negCD27neg
    ## 3    27.00000 AB_CD8AposCD8Bneg CD16negCD27neg
    ## 4    23.64341 AB_CD8AposCD8Bneg CD16negCD27neg
    ## 5    17.05882 AB_CD8AposCD8Bneg CD16negCD27neg
    ## 6    42.45283 AB_CD8AposCD8Bneg CD16negCD27neg

``` r
dat$set <- paste(dat$Tissue, dat$WeeksOfAge, sep = '_')
dat <- select(dat, c(celltype, Percentages, cellsubset, set, SampleID))
dat <- spread(dat, key=cellsubset, value=Percentages)

colnames(dat)
```

    ## [1] "celltype"       "set"            "SampleID"       "CD16negCD27neg"
    ## [5] "CD16negCD27pos" "CD16pos"        "CD16posCD27neg" "CD16posCD27pos"
    ## [9] "CD27"

``` r
for (i in 4:9) {
  formula = as.formula( paste(colnames(dat[i]), 'celltype', sep = '~'))
 stat <- dat %>%
        group_by(set) %>%
        wilcox_test(formula = formula, paired = TRUE) %>% # use paired since tissues derived from same animal
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentParentPopulation_PairedCellTypeWilcoxon_',
                                colnames(dat[i])), '.xlsx'))
}
```

Across only ab primary T-IEL populations:

``` r
subdat <- subset(dat, celltype == 'AB_CD8AposCD8Bpos' | celltype == 'AB_CD8AposCD8Bneg')
colnames(subdat)
```

    ## [1] "celltype"       "set"            "SampleID"       "CD16negCD27neg"
    ## [5] "CD16negCD27pos" "CD16pos"        "CD16posCD27neg" "CD16posCD27pos"
    ## [9] "CD27"

``` r
for (i in 4:9) {
  formula = as.formula( paste(colnames(subdat[i]), 'celltype', sep = '~'))
 stat <- subdat %>%
        group_by(set) %>%
        wilcox_test(formula = formula, paired = TRUE) %>% # use paired since tissues derived from same animal
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentParentPopulation_CD8AposABonly_PairedCellTypeWilcoxon_',
                                colnames(subdat[i])), '.xlsx'))
}
```

Across only gd primary T-IEL populations:

``` r
subdat <- subset(dat, celltype == 'GD_CD2posCD8Apos' | celltype == 'GD_CD2posCD8Aneg')
colnames(subdat)
```

    ## [1] "celltype"       "set"            "SampleID"       "CD16negCD27neg"
    ## [5] "CD16negCD27pos" "CD16pos"        "CD16posCD27neg" "CD16posCD27pos"
    ## [9] "CD27"

``` r
for (i in 4:9) {
  formula = as.formula( paste(colnames(subdat[i]), 'celltype', sep = '~'))
 stat <- subdat %>%
        group_by(set) %>%
        wilcox_test(formula = formula, paired = TRUE) %>% # use paired since tissues derived from same animal
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") 
 write_xlsx(stat, paste0(paste0('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/StatResults/PercentParentPopulation_CD2posGDonly_PairedCellTypeWilcoxon_',
                                colnames(subdat[i])), '.xlsx'))
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
