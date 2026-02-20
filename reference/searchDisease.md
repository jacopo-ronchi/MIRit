# Search for disease EFO identifiers

This function allows to retrieve the Experimental Factor Ontology (EFO)
identifier of a particular disease. This ID is then needed to use the
function
[`findMirnaSNPs()`](https://jacopo-ronchi.github.io/MIRit/reference/findMirnaSNPs.md).

## Usage

``` r
searchDisease(diseaseName)
```

## Arguments

- diseaseName:

  The name of a particular disease (ex. `Alzheimer disease`).

## Value

A `character` object containing EFO identifiers.

## Note

To retrieve EFO IDs for specific diseases, this function makes use of
the `gwasrapidd` package.

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
# search the EFO identifier of Alzheimer disease
searchDisease("Alzheimer disease")
#> Checking for cached EFO traits...
#> Reading EFO traits from cache...
#> Searching for disease: Alzheimer disease
#> [1] "Alzheimer's disease biomarker measurement" 
#> [2] "Alzheimer's disease neuropathologic change"
#> [3] "Alzheimer disease"                         
#> [4] "late-onset Alzheimer's disease"            
#> [5] "family history of Alzheimer’s disease"     
#> [6] "age of onset of Alzheimer disease"         
#> [7] "early-onset Alzheimers disease"            
# }
```
