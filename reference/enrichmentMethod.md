# Access the method used for functional enrichment analyses

This function accesses the `method` slot of a
[`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
object and returns a the name of the enrichment strategy used by the
[`enrichGenes()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichGenes.md)
function to perform the enrichment analysis.

## Usage

``` r
enrichmentMethod(object)
```

## Arguments

- object:

  An object of class
  [`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  containing enrichment results

## Value

A `character` containing the enrichment method, such as `GSEA`.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example FunctionalEnrichment object
obj <- loadExamples("FunctionalEnrichment")

# see the method
enrichmentMethod(obj)
#> [1] "Gene-Set Enrichment Analysis (GSEA)"
```
