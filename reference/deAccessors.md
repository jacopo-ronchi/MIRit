# Extract differentially expressed miRNAs and genes from a [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md) object

The `mirnaDE()` and `geneDE()` are two accessor functions for the
`mirnaDE` and `geneDE` slots of
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
class, respectively. Thus, they can be used to explore the results of
miRNA and gene differential expression analysis stored in a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object.

## Usage

``` r
mirnaDE(object, onlySignificant = TRUE, param = FALSE, returnObject = FALSE)

geneDE(object, onlySignificant = TRUE, param = FALSE, returnObject = FALSE)
```

## Arguments

- object:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

- onlySignificant:

  Logical, if `TRUE` differential expression results will be returned
  just for statistically significant miRNAs/genes, if `FALSE` the full
  table of miRNA/gene differential expression will be provided. Default
  is `TRUE` to only report significant miRNAs/genes

- param:

  Logical, whether to return the complete `list` object with the
  parameters used, or just the results stored in `data`. Default is
  FALSE

- returnObject:

  Logical, if `TRUE` this function will return the limma/edgeR/DESeq2
  object used for differential expression analysis

## Value

A `data.frame` with miRNA/gene differential expression, or a `list`
object with the parameters used if `param = TRUE`.

## Functions

- `mirnaDE()`: Extract miRNA differential expression results

- `geneDE()`: Extract gene differential expression results

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example MirnaExperiment object
obj <- loadExamples()

# access miRNA differential expression of a MirnaExperiment object
sig <- mirnaDE(obj)
all <- mirnaDE(obj, onlySignificant = FALSE)

# access gene differential expression of a MirnaExperiment object
sig <- geneDE(obj)
all <- geneDE(obj, onlySignificant = FALSE)
```
