
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pyhloprocessing

<!-- badges: start -->
<!-- badges: end -->

The goal of pyhloprocessing is to …

## Installation

You can install the development version of pyhloprocessing from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("BJ-Chen-Eric/pyhloprocessing")
```

## Example with in-package dataset

This is a basic example which shows you how to solve a common problem:

``` r
library(pyhloprocessing)
## basic example code: sequence merge and filtering
m_path <- system.file('extdata/raw/gisaid_raw_processing_test.csv', package = 'pyhloprocessing')
gis_fa_path <- system.file('extdata/raw/gisaid_raw_processing_test.fasta', package = 'pyhloprocessing')
ird_fa_path <- system.file('extdata/raw/ird_raw_processing_test.fa', package = 'pyhloprocessing')

raw_processing(meta_path = m_path, 
               gisaid_seq_path = gis_fa_path, 
               ird_seq_path = ird_fa_path , 
               out_fasta_path = paste(getwd(), '/', sep = ''),
               out_fasta_prefix = 'processing_test')
               
## basic example code: Open reading frame identification
fa_path <- system.file('extdata/align/', package = 'pyhloprocessing')

ORF_ident(fasta_path = fa_path, coor_gap_thres = 86, 
          out_fasta_path = paste(getwd(), '/', sep = ''), out_fasta_prefix = 'test_orf')
               
```

<img src="man/figures/test.jeg" width="100%" />
