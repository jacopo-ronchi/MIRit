# Access the miRNA-augmented pathways that were used during TAIPA

This function accesses the `pathways` slot of a
[`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
object and returns a `list` object with the augmented pathways that were
considered by the
[`topologicalAnalysis()`](https://jacopo-ronchi.github.io/MIRit/reference/topologicalAnalysis.md)
function to perform the integrative analysis.

## Usage

``` r
augmentedPathways(object)
```

## Arguments

- object:

  An object of class
  [`IntegrativePathwayAnalysis`](https://jacopo-ronchi.github.io/MIRit/reference/IntegrativePathwayAnalysis-class.md)
  containing the results of a miRNA-mRNA pathway analysis

## Value

A `list` object with the miRNA-augmented biological pathways.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load the example IntegrativePathwayAnalysis object
obj <- loadExamples("IntegrativePathwayAnalysis")

# extract the pathways
ps <- augmentedPathways(obj)
```
