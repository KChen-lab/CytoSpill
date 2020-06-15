
<!-- README.md is generated from README.Rmd. Please edit that file -->
CytoSpill
=========

<!-- badges: start -->
<!-- badges: end -->
The goal of CytoSpill is to compensate the spillover effects in CyTOF data which caused by technical effects without relying on control experiment.

Installation
------------

You can install the development version of CytoSpill from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("qmiao19/CytoSpill")
```

Example
-------

This is a example which shows you how to compensate a CyTOF dataset, we used a sample dataset that included in our package here:

``` r
library(CytoSpill)
getwd()
#> [1] "/Users/qmiao/CytoSpill"
load(file="./data/Levine32.RData")
str(data_Levine32)
#>  num [1:10000, 1:37] 459501 88341 129149 59689 267864 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:37] "Time" "Cell_length" "Ir191Di" "Ir193Di" ...
```

We sampled 5,000 cells from the healthy human bone marrow data used in Levine, 2015. We can check the colunmns of this CyTOF dataset:

``` r
colnames(data_Levine32)
#>  [1] "Time"        "Cell_length" "Ir191Di"     "Ir193Di"     "La139Di"    
#>  [6] "Pr141Di"     "Nd142Di"     "Nd143Di"     "Nd144Di"     "Nd145Di"    
#> [11] "Nd146Di"     "Nd148Di"     "Nd150Di"     "Sm147Di"     "Sm149Di"    
#> [16] "Sm152Di"     "Sm154Di"     "Eu151Di"     "Eu153Di"     "Gd156Di"    
#> [21] "Gd158Di"     "Gd160Di"     "Tb159Di"     "Dy162Di"     "Dy164Di"    
#> [26] "Ho165Di"     "Er166Di"     "Er167Di"     "Er168Di"     "Er170Di"    
#> [31] "Tm169Di"     "Yb171Di"     "Yb172Di"     "Yb174Di"     "Yb176Di"    
#> [36] "Lu175Di"     "Pt195Di"
```

And the corresponding markers:

``` r
markers
#>  [1] "Time"        "Cell_length" "DNA1"        "DNA2"        "CD45RA"     
#>  [6] "CD133"       "CD19"        "CD22"        "CD11b"       "CD4"        
#> [11] "CD8"         "CD34"        "Flt3"        "CD20"        "CXCR4"      
#> [16] "CD235ab"     "CD45"        "CD123"       "CD321"       "CD14"       
#> [21] "CD33"        "CD47"        "CD11c"       "CD7"         "CD15"       
#> [26] "CD16"        "CD44"        "CD38"        "CD13"        "CD3"        
#> [31] "CD61"        "CD117"       "CD49d"       "HLA-DR"      "CD64"       
#> [36] "CD41"        "Viability"
```

After checking the markers we decided to compensate for columns 5 to 36. We can use the function GetSpillMat to obtain the compensation matrix. We need to specify the data we are using, which column are used to calculate the spillover effects and using how many cells for calculation.

``` r
GetSpillMat(data = data_Levine32, cols = c(5:36), n = 5000)
```

Or we can directly obtain the compensated data and write it as a fcs format file using SpillComp function. In addition, we need to give the function an output file name.

``` r
SpillComp(data = data_Levine32, cols = c(5:36), n = 5000, output = 'Levine32_compensation.fcs')
```
