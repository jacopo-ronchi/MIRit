# The constructor function for [MirnaExperiment](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)

This is the constructor function that allows to easily create objects of
class
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md).
This function requires as inputs miRNA and gene expression matrices, as
well as sample metadata.

## Usage

``` r
MirnaExperiment(mirnaExpr, geneExpr, samplesMetadata, pairedSamples = TRUE)
```

## Arguments

- mirnaExpr:

  A `matrix` object containing microRNA expression levels. Other objects
  coercible to `matrix` are also accepted (e.g. `data.frame`). This
  object must be structured as specified in the *details* section

- geneExpr:

  A `matrix` object containing gene expression levels. Other objects
  coercible to `matrix` are also accepted (e.g. `data.frame`). This
  object must be structured as specified in the *details* section

- samplesMetadata:

  A `data.frame` object containing information about samples used for
  microRNA and gene expression profiling. For further information see
  the *details* section

- pairedSamples:

  Logical, whether miRNA and gene expression levels derive from the same
  subjects or not. Check the *details* section for additional
  instructions. Default is `TRUE`

## Value

A valid
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object containing information about miRNA and gene expression.

## Details

This function requires data to be prepared as described below.

### mirnaExpr and geneExpr

`mirnaExpr` and `geneExpr` must be `matrix` objects (or objects
coercible to one) that contain miRNA and gene expression values,
respectively. Rows must represent the different miRNAs/genes analyzed
while columns must represent the different samples in study. For
`mirnaExpr`, row names must contain miRNA names according to miRBase
nomenclature, whereas for `geneExpr`, row names must contain gene
symbols according to hgnc nomenclature. The values contained in these
objects can derive from both microarray and RNA-Seq experiments.

For NGS experiments, `mirnaExpr` and `geneExpr` should just be
un-normalized count matrices. Instead, for microarray experiments, data
should be normalized and log2 transformed, for example with the RMA
algorithm.

### samplesMetadata

`samplesMetadata` must be a `data.frame` object containing information
about samples used for miRNA profiling and for gene expression analysis.
Specifically, this data.frame must contain:

- A column named `primary`, specifying an identifier for each sample;

- A column named `mirnaCol`, containing the column names used for each
  sample in the `mirnaExpr` object;

- A column named `geneCol`, containing the column names used for each
  sample in the `geneExpr` object;

- Other eventual columns that define specific sample metadata, such as
  disease condition, age, sex and so on...

For unpaired samples, NAs can be used for missing entries in
`mirnaCol`/`geneCol`.

### pairedSamples

MicroRNA and gene expression measurements may derive from the same
subjects (i.e. samples used to generate both miRNA and gene expression
data) or from different individuals (i.e. miRNA expression assayed on a
group of samples and gene expression retrieved from a different group of
samples). `pairedSamples` is a `logical` parameter that defines the
relationship between miRNA and gene expression measurements. It must be
`TRUE` if data derive from the same individuals, while it must be
`FALSE` when data derive from different subjects.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example data
data(geneCounts, package = "MIRit")
data(mirnaCounts, package = "MIRit")

# create samples metadata
meta <- data.frame(
    "primary" = colnames(geneCounts),
    "mirnaCol" = colnames(mirnaCounts), "geneCol" = colnames(geneCounts),
    "disease" = c(rep("PTC", 8), rep("NTH", 8)),
    "patient" = c(rep(paste("Sample_", seq(8), sep = ""), 2))
)

# create a 'MirnaExperiment' object
obj <- MirnaExperiment(
    mirnaExpr = mirnaCounts, geneExpr = geneCounts,
    samplesMetadata = meta, pairedSamples = TRUE
)
```
