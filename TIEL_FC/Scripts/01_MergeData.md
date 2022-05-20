Merging Data
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

## Load required packages

See session information at bottom for further information.

``` r
library(Spectre)
library(readxl)
library(data.table)
```

## Set directory:

``` r
dirData <- '/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/CD3EposOnly_FlowJo_ChannelValues/'
setwd('/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC') # set working directory
```

## Import flow data files:

``` r
list.files(dirData, ".csv") # See the list of files present in the input directory
```

    ##  [1] "13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ##  [2] "13il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ##  [3] "13jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ##  [4] "17ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ##  [5] "17il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ##  [6] "17jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ##  [7] "1ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"  
    ##  [8] "1il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"  
    ##  [9] "1jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [10] "21ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [11] "21il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [12] "21jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [13] "25ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [14] "25il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [15] "25jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [16] "29ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [17] "29il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [18] "29jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [19] "33ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [20] "33il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [21] "33jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [22] "37ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [23] "37il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [24] "37jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [25] "41ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [26] "41il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [27] "41jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [28] "45ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [29] "45il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [30] "45jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [31] "49ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [32] "49il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [33] "49jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [34] "53ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [35] "53il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [36] "53jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [37] "57ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [38] "57il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [39] "57jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [40] "5ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"  
    ## [41] "5il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"  
    ## [42] "5jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [43] "61ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [44] "61il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209.csv" 
    ## [45] "61jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"
    ## [46] "9ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"  
    ## [47] "9il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"  
    ## [48] "9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209.csv"

``` r
data.list <- Spectre::read.files(file.loc = dirData,
                                 file.type = ".csv",
                                 do.embed.file.names = TRUE) # read in all .csv data files as a list

names(data.list) # names in data list should correspond to listed files in dirData, minus the .csv suffix
```

    ##  [1] "13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ##  [2] "13il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ##  [3] "13jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ##  [4] "17ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ##  [5] "17il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ##  [6] "17jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ##  [7] "1ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209"  
    ##  [8] "1il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209"  
    ##  [9] "1jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [10] "21ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [11] "21il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [12] "21jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [13] "25ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [14] "25il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [15] "25jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [16] "29ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [17] "29il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [18] "29jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [19] "33ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [20] "33il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [21] "33jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [22] "37ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [23] "37il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [24] "37jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [25] "41ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [26] "41il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [27] "41jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [28] "45ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [29] "45il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [30] "45jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [31] "49ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [32] "49il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [33] "49jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [34] "53ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [35] "53il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [36] "53jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [37] "57ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [38] "57il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [39] "57jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [40] "5ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209"  
    ## [41] "5il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209"  
    ## [42] "5jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [43] "61ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [44] "61il_E8_E08_008_CD3e+_CD3EposOnly_ChannelValues_20220209" 
    ## [45] "61jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelValues_20220209"
    ## [46] "9ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelValues_20220209"  
    ## [47] "9il_E2_E02_002_CD3e+_CD3EposOnly_ChannelValues_20220209"  
    ## [48] "9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209"

## QC check flow data:

``` r
check <- do.list.summary(data.list)
check$name.table # Review column names and their subsequent values
```

    ##       X1    X2    X3    X4    X5    X6                 X7                    X8
    ## 1  FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 2  FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 3  FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 4  FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 5  FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 6  FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 7  FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 8  FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 9  FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 10 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 11 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 12 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 13 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 14 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 15 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 16 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 17 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 18 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 19 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 20 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 21 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 22 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 23 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 24 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 25 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 26 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 27 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 28 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 29 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 30 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 31 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 32 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 33 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 34 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 35 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 36 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 37 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 38 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 39 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 40 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 41 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 42 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 43 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 44 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 45 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 46 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 47 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ## 48 FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a Comp-APC-Cy7-A :: L_D
    ##                      X9                 X10                   X11
    ## 1  Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 2  Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 3  Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 4  Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 5  Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 6  Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 7  Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 8  Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 9  Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 10 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 11 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 12 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 13 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 14 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 15 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 16 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 17 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 18 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 19 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 20 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 21 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 22 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 23 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 24 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 25 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 26 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 27 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 28 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 29 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 30 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 31 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 32 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 33 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 34 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 35 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 36 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 37 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 38 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 39 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 40 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 41 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 42 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 43 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 44 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 45 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 46 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 47 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ## 48 Comp-BB515-A :: CD27 Comp-BB700-A :: CD4 Comp-BUV496-A :: CD16
    ##                    X12               X13                      X14
    ## 1  Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 2  Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 3  Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 4  Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 5  Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 6  Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 7  Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 8  Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 9  Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 10 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 11 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 12 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 13 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 14 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 15 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 16 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 17 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 18 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 19 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 20 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 21 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 22 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 23 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 24 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 25 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 26 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 27 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 28 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 29 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 30 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 31 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 32 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 33 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 34 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 35 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 36 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 37 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 38 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 39 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 40 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 41 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 42 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 43 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 44 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 45 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 46 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 47 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ## 48 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b Comp-PE-CF594-A :: gdTCR
    ##                      X15      X16    X17
    ## 1  Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 2  Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 3  Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 4  Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 5  Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 6  Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 7  Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 8  Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 9  Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 10 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 11 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 12 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 13 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 14 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 15 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 16 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 17 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 18 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 19 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 20 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 21 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 22 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 23 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 24 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 25 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 26 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 27 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 28 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 29 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 30 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 31 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 32 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 33 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 34 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 35 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 36 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 37 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 38 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 39 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 40 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 41 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 42 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 43 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 44 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 45 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 46 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 47 Comp-PE-Cy7-A :: CD3e FileName FileNo
    ## 48 Comp-PE-Cy7-A :: CD3e FileName FileNo

``` r
#check$ncol.check # Review number of columns (features/markers) in each sample
#check$nrow.check # Review number of rows (cells) in each sample
#data.list[[1]] # look at top 5 rows & bottom 5 rows of first sample in list
rm(check) # free up space
```

## Merge data

``` r
cell.dat <- Spectre::do.merge.files(dat = data.list) # if all looks good from QC, merge data into a single table
cell.dat # look at top 5 rows & bottom 5 rows of merged table
```

    ##          FSC-A FSC-H FSC-W SSC-A SSC-H SSC-W Comp-APC-A :: CD8a
    ##       1:   360   335   592   574   385     2                846
    ##       2:   484   438   636   656   436   108                858
    ##       3:   359   332   605   589   393     2                584
    ##       4:   323   304   597   506   350     4                749
    ##       5:   387   363   600   593   402    63                467
    ##      ---                                                       
    ## 1346771:   330   316   576   483   338     1                737
    ## 1346772:   305   278   583   543   369     1                891
    ## 1346773:   288   277   562   464   330     2                694
    ## 1346774:   306   286   578   473   338     0                685
    ## 1346775:   253   217   568   618   408    76                796
    ##          Comp-APC-Cy7-A :: L_D Comp-BB515-A :: CD27 Comp-BB700-A :: CD4
    ##       1:                   540                  319                   0
    ##       2:                   549                  334                   0
    ##       3:                   556                  260                   0
    ##       4:                   468                  149                 173
    ##       5:                   402                  431                 251
    ##      ---                                                               
    ## 1346771:                   470                  341                 427
    ## 1346772:                   591                  599                 223
    ## 1346773:                   445                  151                 334
    ## 1346774:                   493                   95                   0
    ## 1346775:                   517                  124                 485
    ##          Comp-BUV496-A :: CD16 Comp-BV650-A :: CD2 Comp-PE-A :: CD8b
    ##       1:                    99                 761               527
    ##       2:                   200                 830               610
    ##       3:                   584                 814               165
    ##       4:                   321                 761               520
    ##       5:                   746                 777               207
    ##      ---                                                            
    ## 1346771:                   107                 719               694
    ## 1346772:                   206                 899               745
    ## 1346773:                   105                 719               420
    ## 1346774:                   137                 788               534
    ## 1346775:                   196                 815               758
    ##          Comp-PE-CF594-A :: gdTCR Comp-PE-Cy7-A :: CD3e
    ##       1:                      341                   296
    ##       2:                      332                   350
    ##       3:                      554                   750
    ##       4:                      550                   659
    ##       5:                      400                   452
    ##      ---                                               
    ## 1346771:                      329                   388
    ## 1346772:                      317                   421
    ## 1346773:                      353                   418
    ## 1346774:                      342                   389
    ## 1346775:                      322                   469
    ##                                                          FileName FileNo
    ##       1: 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209      1
    ##       2: 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209      1
    ##       3: 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209      1
    ##       4: 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209      1
    ##       5: 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelValues_20220209      1
    ##      ---                                                                
    ## 1346771: 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209     48
    ## 1346772: 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209     48
    ## 1346773: 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209     48
    ## 1346774: 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209     48
    ## 1346775: 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelValues_20220209     48

``` r
rm(data.list)
```

## Rename data columns

``` r
colnames(cell.dat) # list current column names
```

    ##  [1] "FSC-A"                    "FSC-H"                   
    ##  [3] "FSC-W"                    "SSC-A"                   
    ##  [5] "SSC-H"                    "SSC-W"                   
    ##  [7] "Comp-APC-A :: CD8a"       "Comp-APC-Cy7-A :: L_D"   
    ##  [9] "Comp-BB515-A :: CD27"     "Comp-BB700-A :: CD4"     
    ## [11] "Comp-BUV496-A :: CD16"    "Comp-BV650-A :: CD2"     
    ## [13] "Comp-PE-A :: CD8b"        "Comp-PE-CF594-A :: gdTCR"
    ## [15] "Comp-PE-Cy7-A :: CD3e"    "FileName"                
    ## [17] "FileNo"

``` r
colnames(cell.dat) <- c('FSC_A', 'FSC_H', 'FSC_W', 'SSC_A', 'SSC_H', 'SSC_W', # new column names
                        'CD8A_APC', 'Viability_eFluor480', 'CD27_FITC', 
                        'CD4_PerCPCy5.5', 'CD16_BUV496', 'CD2_BV650', 
                        'CD8B_PE', 'gdTCR_iFluor594', 'CD3E_PECy7',
                        'FileName', 'FileNo')
colnames(cell.dat)
```

    ##  [1] "FSC_A"               "FSC_H"               "FSC_W"              
    ##  [4] "SSC_A"               "SSC_H"               "SSC_W"              
    ##  [7] "CD8A_APC"            "Viability_eFluor480" "CD27_FITC"          
    ## [10] "CD4_PerCPCy5.5"      "CD16_BUV496"         "CD2_BV650"          
    ## [13] "CD8B_PE"             "gdTCR_iFluor594"     "CD3E_PECy7"         
    ## [16] "FileName"            "FileNo"

## Data transformation

We skip this step since we used .csv files with channel values exported
from FlowJo, which incorporates axis scaling that we toggled with in
FlowJo. See this link as a great resource for explanation:

<https://wiki.centenary.org.au/display/SPECTRE/Data+transformation>
<https://wiki.centenary.org.au/display/SPECTRE/Exporting+data+from+FlowJo+for+analysis+in+Spectre>

## Add metadata to cell dataframe

``` r
meta.dat <- read_excel("/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/NG10_MetaData_FC.xlsx") # read in meta data
meta.dat # look at meta data... if you don't want to include all metadata columns, subset before next step
```

    ## # A tibble: 48 × 5
    ##    FileName                                     Batch WeeksOfAge Tissue AnimalID
    ##    <chr>                                        <chr>      <dbl> <chr>  <chr>   
    ##  1 1ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelVal… A              5 cecum  1       
    ##  2 1il_E2_E02_002_CD3e+_CD3EposOnly_ChannelVal… A              5 ileum  1       
    ##  3 1jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelVa… A              5 jejun… 1       
    ##  4 5ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelVal… A              5 cecum  5       
    ##  5 5il_E8_E08_008_CD3e+_CD3EposOnly_ChannelVal… A              5 ileum  5       
    ##  6 5jej_E7_E07_007_CD3e+_CD3EposOnly_ChannelVa… A              5 jejun… 5       
    ##  7 9ce_E3_E03_003_CD3e+_CD3EposOnly_ChannelVal… B              5 cecum  9       
    ##  8 9il_E2_E02_002_CD3e+_CD3EposOnly_ChannelVal… B              5 ileum  9       
    ##  9 9jej_E1_E01_001_CD3e+_CD3EposOnly_ChannelVa… B              5 jejun… 9       
    ## 10 13ce_E9_E09_009_CD3e+_CD3EposOnly_ChannelVa… B              5 cecum  13      
    ## # … with 38 more rows

``` r
meta.dat$FileName <- gsub(".csv","",as.character(meta.dat$FileName)) # remove .csv suffix so FileName in meta.dat matches FileName in cell.dat
cell.dat <- do.add.cols(dat = cell.dat, base.col = "FileName", add.dat = meta.dat, add.by = "FileName") # add metadata based on corresponding file names between cell.dat and meta.dat
```

    ## Step 1/3. Mapping data

    ## Step 2/3. Merging data

    ## Step 3/3. Returning data

``` r
cell.dat # check that new meta data is present
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

``` r
rm(meta.dat)
```

## Save merged data:

``` r
fwrite(cell.dat, "/home/Jayne.Wiarda/NG10/Dissertation/TIEL_FC/AllCellsMerged_InputForCorrection.csv") # save merged data into a single .csv file
```

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
    ## [1] data.table_1.14.0 readxl_1.3.1      Spectre_0.4.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7       rstudioapi_0.13  knitr_1.34       magrittr_2.0.1  
    ##  [5] rlang_0.4.11     fastmap_1.1.0    fansi_0.5.0      stringr_1.4.0   
    ##  [9] tools_4.1.3      xfun_0.26        utf8_1.2.2       cli_3.0.1       
    ## [13] htmltools_0.5.2  ellipsis_0.3.2   yaml_2.2.1       digest_0.6.27   
    ## [17] tibble_3.1.4     lifecycle_1.0.0  crayon_1.4.1     vctrs_0.3.8     
    ## [21] evaluate_0.14    rmarkdown_2.11   stringi_1.7.4    compiler_4.1.3  
    ## [25] pillar_1.6.2     cellranger_1.1.0 pkgconfig_2.0.3
