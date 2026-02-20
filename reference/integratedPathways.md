# Access the results of integrative miRNA-mRNA pathway analyses

This function accesses the `data` slot of a
[`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
object and returns a `data.frame` with the results of an integrative
topological analysis carried out through the
[`topologicalAnalysis()`](https://jacopo-ronchi.github.io/MIRit/reference/topologicalAnalysis.md)
function.

## Usage

``` r
integratedPathways(object)
```

## Arguments

- object:

  An object of class
  [`IntegrativePathwayAnalysis`](https://jacopo-ronchi.github.io/MIRit/reference/IntegrativePathwayAnalysis-class.md)
  containing the results of a miRNA-mRNA pathway analysis

## Value

A `data.frame` object containing the results of the topological
analysis, as returned by the
[`topologicalAnalysis()`](https://jacopo-ronchi.github.io/MIRit/reference/topologicalAnalysis.md)
function.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load the example IntegrativePathwayAnalysis object
obj <- loadExamples("IntegrativePathwayAnalysis")

# extract results
taipaRes <- integratedPathways(obj)
```
