
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

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(pyhloprocessing)
## basic example code: sequence merge and filtering
raw_processing(meta_path = here::here('test_data/raw/', 'gisaid_raw_processing_test.csv'), 
               gisaid_seq_path = here::here('test_data/raw/', 'gisaid_raw_processing_test.fasta'), 
               ird_seq_path = here::here('test_data/raw/', 'ird_raw_processing_test.fa') , 
               out_fasta_path = paste(getwd(), '/', sep = ''),
               out_fasta_prefix = 'processing_test')
               
## basic example code: Open reading frame identification
ORF_ident(fasta_path = here::here('test_data/align/'), coor_gap_thres = 86, 
          out_fasta_path = '~/', out_fasta_prefix = 'test_orf')
               
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r

```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/test.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
