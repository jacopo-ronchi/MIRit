
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MIRit <img src="man/figures/logo.svg" align="right" height="139" alt="" />

<!-- badges: start -->

[![Devel
version](https://img.shields.io/badge/devel%20version-0.99.4-blue.svg)](https://github.com/jacopo-ronchi/MIRit)
[![GitHub
issues](https://img.shields.io/github/issues/jacopo-ronchi/MIRit)](https://github.com/jacopo-ronchi/MIRit/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/jacopo-ronchi/MIRit)](https://github.com/jacopo-ronchi/MIRit/pulls)
[![Last
commit](https://img.shields.io/github/last-commit/jacopo-ronchi/MIRit.svg)](https://github.com/jacopo-ronchi/MIRit/commits/devel)
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL (\>=
3)](https://img.shields.io/badge/license-GPL%20(%3E=%203)-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![DOI
badge](https://img.shields.io/badge/doi-10.1101/2023.11.24.568528-yellow.svg)](https://doi.org/10.1101/2023.11.24.568528)
[![R-CMD-check-bioc](https://github.com/jacopo-ronchi/MIRit/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/jacopo-ronchi/MIRit/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/jacopo-ronchi/MIRit/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/jacopo-ronchi/MIRit?branch=devel)
<!-- badges: end -->

## Overview

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

## Authors

**Dr. Jacopo Ronchi**
<a itemprop="sameAs" content="https://orcid.org/0000-0001-5520-4631" href="https://orcid.org/0000-0001-5520-4631" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a><sup>1</sup>
(author and maintainer)

**Dr. Maria Foti**
<a itemprop="sameAs" content="https://orcid.org/0000-0002-4481-1900" href="https://orcid.org/0000-0002-4481-1900" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a><sup>1</sup>

<sup>1</sup>School of Medicine and Surgery, University of
Milano-Bicocca, Italy

## Installation

`MIRit` is currently undergoing Bioconductor submission. In the
meantime, you can install it from
[GitHub](https://github.com/jacopo-ronchi/MIRit) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

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

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
