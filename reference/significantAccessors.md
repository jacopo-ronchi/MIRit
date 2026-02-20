# Get the IDs of statistically differentially expressed miRNAs/genes

The `significantMirnas()` and `significantGenes()` functions access the
`significant` features contained in the `mirnaDE` or `geneDE` slots of a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object, and can be used to obtain the IDs of statistically
differentially expressed miRNAs and genes.

## Usage

``` r
significantMirnas(object)

significantGenes(object)
```

## Arguments

- object:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

## Value

A `character` vector of miRNA IDs (e.g. 'hsa-miR-16-5p',
hsa-miR-29a-3p'...), or a`character` vector of gene symbols (e.g.
'TP53', 'FOXP2', 'TIGAR', CASP1'...).

## Functions

- `significantMirnas()`: Get the IDs of differentially expressed miRNAs

- `significantGenes()`: Get the IDs of differentially expressed genes

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example MirnaExperiment object
obj <- loadExamples()

# extract significant DE-miRNAs
sigMirnas <- significantMirnas(obj)

# extract significant DEGs
sigGenes <- significantGenes(obj)
```
