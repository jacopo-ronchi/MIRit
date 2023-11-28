
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MIRit

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/jacopo-ronchi/MIRit)](https://github.com/jacopo-ronchi/MIRit/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/jacopo-ronchi/MIRit)](https://github.com/jacopo-ronchi/MIRit/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/MIRit.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MIRit)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/MIRit.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/MIRit)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/MIRit.svg)](http://bioconductor.org/packages/stats/bioc/MIRit/)
[![Bioc
support](https://bioconductor.org/shields/posts/MIRit.svg)](https://support.bioconductor.org/tag/MIRit)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/MIRit.svg)](https://bioconductor.org/packages/release/bioc/html/MIRit.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/MIRit.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/MIRit/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/MIRit.svg)](https://bioconductor.org/packages/release/bioc/html/MIRit.html#since)
[![R-CMD-check-bioc](https://github.com/jacopo-ronchi/MIRit/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/jacopo-ronchi/MIRit/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/jacopo-ronchi/MIRit/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/jacopo-ronchi/MIRit?branch=devel)
<!-- badges: end -->

## Overview

<img src="man/figures/logo.png" align="right" width="330" style="float:right; width:330px;">

`MIRit` (miRNA integration tool) is an open-source R package that aims
to facilitate the comprehension of microRNA (miRNA) biology through the
integrative analysis of gene and miRNA expression data deriving from
different platforms, including microarrays, RNA-Seq, miRNA-Seq,
proteomics and single-cell transcriptomics. Given their regulatory
importance, a complete characterization of miRNA dysregulations results
crucial to explore the molecular networks that may lead to the
insurgence of complex diseases. To this purpose, we developed MIRit, an
all-in-one framework that provides flexible and powerful methods for
performing integrative miRNA-mRNA multi-omic analyses from start to
finish.

## Author

**Dr. Jacopo Ronchi**
<a itemprop="sameAs" content="https://orcid.org/0000-0001-5520-4631" href="https://orcid.org/0000-0001-5520-4631" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a><sup>1</sup>

<sup>1</sup>School of Medicine and Surgery, University of
Milano-Bicocca, Italy

## Installation

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `MIRit` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("MIRit")
```

Alternatively, the development version of `MIRit` can be installed from
[GitHub](https://github.com/jacopo-ronchi/MIRit) with:

``` r
BiocManager::install("jacopo-ronchi/MIRit")
```

## Usage

For detailed instructions on how to use `MIRit` for integrative
miRNA-mRNA analysis, please refer to the package vignette on
[Bioconductor](). Alternatively, you can refer to the [documentation
website](http://jacopo-ronchi.github.io/MIRit).

## Citation

If you use `MIRit` in published research, please cite the corresponding
paper:

> Ronchi J and Foti M. ‘MIRit: an integrative R framework for the
> identification of impaired miRNA-mRNA regulatory networks in complex
> diseases’. bioRxiv (2023). <doi:10.1101/2023.11.24.568528>

Please note that the `MIRit` package was made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `MIRit` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.17/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://jacopo-ronchi.github.io/MIRit) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
