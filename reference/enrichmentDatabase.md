# Access the database used for functional enrichment analyses

This function accesses the `database` slot of a
[`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
object and returns a the name of the database used by the
[`enrichGenes()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichGenes.md)
function to perform the enrichment analysis.

## Usage

``` r
enrichmentDatabase(object)
```

## Arguments

- object:

  An object of class
  [`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  containing enrichment results

## Value

A `character` containing the name of the database, such as `KEGG`.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example FunctionalEnrichment object
obj <- loadExamples("FunctionalEnrichment")

# see the database
enrichmentDatabase(obj)
#> [1] "KEGG (category: pathway)"
```
