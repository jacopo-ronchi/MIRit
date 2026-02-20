# Access the results of functional enrichment analyses

This function accesses the `data` slot of a
[`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
object and returns a `data.frame` with enrichment results.

## Usage

``` r
enrichmentResults(object)
```

## Arguments

- object:

  An object of class
  [`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  containing enrichment results

## Value

A `data.frame` object containing the results of functional enrichment
analyses, as returned by the
[`enrichGenes()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichGenes.md)
function.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example FunctionalEnrichment object
obj <- loadExamples("FunctionalEnrichment")

# extract results
de_df <- enrichmentResults(obj)
```
