# MIRit

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

**Dr. Jacopo Ronchi** [![ORCID iD
icon](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-5520-4631)¹
(author and maintainer)

**Dr. Maria Foti** [![ORCID iD
icon](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-4481-1900)¹

¹School of Medicine and Surgery, University of Milano-Bicocca, Italy

## Installation

`MIRit` is available on Bioconductor and can be installed using:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("MIRit")
```

Alternatively, the development version of MIRit can be installed from
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
[Bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/MIRit/inst/doc/MIRit.html).
Alternatively, you can refer to the [documentation
website](http://jacopo-ronchi.github.io/MIRit).

## Citation

If you use `MIRit` in published research, please cite the corresponding
paper:

> Ronchi, J., & Foti, M. (2026). MIRit: An integrative R framework for
> the identification of impaired miRNA–mRNA regulatory networks in
> complex diseases. Bioinformatics Advances, vbag042.
> <https://doi.org/10.1093/bioadv/vbag042>

Please note that the `MIRit` package was made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `MIRit` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
