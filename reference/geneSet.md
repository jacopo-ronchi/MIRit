# Extract the gene-sets used for functional enrichment analyses

This function accesses the `geneSet` slot of a
[`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
object and returns a `list` with the collection of genes used for the
enrichment.

## Usage

``` r
geneSet(object)
```

## Arguments

- object:

  An object of class
  [`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  containing enrichment results

## Value

A `list` containing the gene-sets.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example FunctionalEnrichment object
obj <- loadExamples("FunctionalEnrichment")

# extract the gene-sets
gs <- geneSet(obj)
```
