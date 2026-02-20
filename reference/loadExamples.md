# Load example MIRit objects

This helper function allows to create a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object containing miRNA and gene expression data deriving from
Riesco-Eizaguirre et al (2015), an
[`IntegrativePathwayAnalysis`](https://jacopo-ronchi.github.io/MIRit/reference/IntegrativePathwayAnalysis-class.md)
object containing TAIPA results for the same dataset, or a
[`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
with example GSEA enrichment results.

## Usage

``` r
loadExamples(class = "MirnaExperiment")
```

## Arguments

- class:

  It must be `MirnaExperiment` (default) to load an example object of
  class
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md),
  `IntegrativePathwayAnalysis`, to load an example object of class
  [`IntegrativePathwayAnalysis`](https://jacopo-ronchi.github.io/MIRit/reference/IntegrativePathwayAnalysis-class.md),
  or `FunctionalEnrichment`, to load an example object of class
  [`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md).

## Value

An example `MirnaExperiment` object, an `IntegrativePathwayAnalysis`
object, or a `FunctionalEnrichment` object.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example MirnaExperiment object
obj <- loadExamples()

# load example IntegrativePathwayAnalysis object
obj <- loadExamples("IntegrativePathwayAnalysis")

# load example FunctionalEnrichment object
obj <- loadExamples("FunctionalEnrichment")
```
