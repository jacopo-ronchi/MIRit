# Extract the database used for integrative miRNA-mRNA pathway analyses

This function accesses the `database` slot of a
[`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
object and returns the name of the database used by the
[`topologicalAnalysis()`](https://jacopo-ronchi.github.io/MIRit/reference/topologicalAnalysis.md)
function to perform the integrative topological analysis.

## Usage

``` r
integrationDatabase(object)
```

## Arguments

- object:

  An object of class
  [`IntegrativePathwayAnalysis`](https://jacopo-ronchi.github.io/MIRit/reference/IntegrativePathwayAnalysis-class.md)
  containing the results of a miRNA-mRNA pathway analysis

## Value

A `character` object with the name of the database used by the
[`topologicalAnalysis()`](https://jacopo-ronchi.github.io/MIRit/reference/topologicalAnalysis.md)
function, such as `KEGG`.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load the example IntegrativePathwayAnalysis object
obj <- loadExamples("IntegrativePathwayAnalysis")

# see the database
integrationDatabase(obj)
#> [1] "KEGG"
```
