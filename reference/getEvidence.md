# Get the scientific evidence for a particular disease-SNP association

This function returns the biomedical evidence that supports the
association between a particular SNP and a phenotypic trait.

## Usage

``` r
getEvidence(variant, diseaseEFO)
```

## Arguments

- variant:

  The SNP ID of a particular variant of interest (e.g. 'rs394581')

- diseaseEFO:

  The EFO identifier of a disease of interest. This can be identified
  with the function
  [`searchDisease()`](https://jacopo-ronchi.github.io/MIRit/reference/searchDisease.md)

## Value

A `tbl_df` dataframe containing information about literature evidences
for a disease-SNP association.

## Note

To retrieve evidences of disease-SNP association, this function makes
use of the `gwasrapidd` package.

## References

Ramiro Magno, Ana-Teresa Maia, gwasrapidd: an R package to query,
download and wrangle GWAS catalog data, Bioinformatics, Volume 36, Issue
2, January 2020, Pages 649–650,
<https://doi.org/10.1093/bioinformatics/btz605>.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# \donttest{
# searchDisease("Alzheimer disease")
evidence <- getEvidence("rs2075650", diseaseEFO = "Alzheimer disease")
#> Retrieving biomedical evidence for the association between Alzheimer disease and rs2075650 variant...
#> 137 studies reporting this association were found!
# }
```
