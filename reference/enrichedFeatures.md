# Extract the names of the pre-ranked features in a GSEA experiment

This function accesses the `features` slot of a
[`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
object and returns a `character` vector with the names of the features
considered in GSEA in the same order as the ranking metric.

## Usage

``` r
enrichedFeatures(object)
```

## Arguments

- object:

  An object of class
  [`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  containing enrichment results

## Value

A `character` vector with the names of the genes ordered based on the
ranking metric.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example FunctionalEnrichment object
obj <- loadExamples("FunctionalEnrichment")

# extract the ranking metric
rmet <- enrichmentMetric(obj)

## extract the corresponding names
rnames <- enrichedFeatures(obj)
```
