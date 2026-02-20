# Manually add differential expression results to a MirnaExperiment object

This function allows to add miRNA and gene differential expression
results to a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object. Instead of running
[`performMirnaDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAnalysis.md)
and
[`performGeneDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAnalysis.md)
functions, this one allows to use differential expression analyses
carried out in other ways. Note that it is possible to manually add
differential expression results just for miRNAs or just for genes. This
is particularly useful in order to use the pipeline implemented in MIRit
for proteomic data and for expression data deriving from different
technologies.

## Usage

``` r
addDifferentialExpression(
  mirnaObj,
  mirnaDE = NULL,
  geneDE = NULL,
  mirna.logFC = 0,
  mirna.pCutoff = 0.05,
  mirna.pAdjustment = "fdr",
  gene.logFC = 0,
  gene.pCutoff = 0.05,
  gene.pAdjustment = "fdr"
)
```

## Arguments

- mirnaObj:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

- mirnaDE:

  A `data.frame` containing the output of miRNA differential expression
  analysis. Check the *details* section to see the required format.
  Default is NULL not to add miRNA differential expression results

- geneDE:

  A `data.frame` containing the output of gene differential expression
  analysis. Check the *details* section to see the required format.
  Default is NULL not to add gene differential expression results

- mirna.logFC:

  The minimum log2 fold change required to consider a miRNA as
  differentially expressed. Optional, default is 0

- mirna.pCutoff:

  The adjusted p-value cutoff to use for miRNA statistical significance.
  The default value is `0.05`

- mirna.pAdjustment:

  The p-value correction method for miRNA multiple testing. It must be
  one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
  `bonferroni`, `BY`

- gene.logFC:

  The minimum log2 fold change required to consider a gene as
  differentially expressed. Optional, default is 0

- gene.pCutoff:

  The adjusted p-value cutoff to use for gene statistical significance.
  The default value is `0.05`

- gene.pAdjustment:

  The p-value correction method for gene multiple testing. It must be
  one of: `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
  `bonferroni`, `BY`

## Value

A
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object containing differential expression results. To access these
results, the user may run the
[`mirnaDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md)
and
[`geneDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md)
functions for miRNAs and genes, respectively.

## Details

The following paragraphs briefly explain the formats needed for mirnaDE,
geneDE, and differential expression parameters.

### mirnaDE and geneDE

`mirnaDE` and `geneDE` are two objects of class `data.frame` containing
the results of miRNA and gene differential expression analysis
respectively. These tables should contain the differential expression
results for all miRNAs/genes analyzed, not just for statistically
significant species.

Note that you can individually add differential expression results for
miRNAs and genes. For instance, it is possible to manually add gene
differential expression through this function, while performing miRNA
differential expression through the
[`performMirnaDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAnalysis.md)
function, and vice versa. In order to only add miRNA or gene
differential expression results, you must leave `mirnaDE` or `geneDE`
slots to NULL.

All `data.frame` objects can be used, as long as they have:

- One column containing miRNA/gene names (according to miRBase/hgnc
  nomenclature). Accepted column names are: `ID`, `Symbol`,
  `Gene_Symbol`, `Mirna`, `mir`, `Gene`, `gene.symbol`, `Gene.symbol`;

- One column with log2 fold changes. Accepted column names are: `logFC`,
  `log2FoldChange`, `FC`, `lFC`;

- One column with average expression. Accepted column names are:
  `AveExpr`, `baseMean`, `logCPM`;

- One column with the p-values resulting from the differential
  expression analysis. Accepted column names are: `P.Value`, `pvalue`,
  `PValue`, `Pvalue`;

- One column containing p-values adjusted for multiple testing. Accepted
  column names are: `adj.P.Val`, `padj`, `FDR`, `fdr`, `adj`, `adj.p`,
  `adjp`.

### Differential expression cutoffs

`mirna.logFC`, `mirna.pCutoff`, `mirna.pAdjustment`, and `gene.logFC`,
`gene.pCutoff`, `gene.pAdjustment` represent the parameters used to
define the significance of differential expression results. These are
needed in order to inform MIRit about the features that are considered
as differentially expressed.

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

# load example tables with differential expression results
de_m <- mirnaDE(object = loadExamples(), onlySignificant = FALSE)
de_g <- geneDE(object = loadExamples(), onlySignificant = FALSE)

# add DE results to MirnaExperiment object
obj <- addDifferentialExpression(obj, de_m, de_g,
    mirna.pCutoff = 0.05,
    gene.pCutoff = 0.05
)
```
