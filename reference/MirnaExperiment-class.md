# The 'MirnaExperiment' class

This class extends the
[`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment-class.html)
from the homonym package to provide the flexibility of handling genomic
data of both microRNAs and their targets, allowing to store information
about microRNA and gene expression, differential expression results,
microRNA targets and miRNA-gene integration analysis.

This class can be used to manage genomic data deriving from different
sources, like RNA-Seq, microarrays and mass spectrometry. Moreover,
microRNA and gene expression levels may originate from the same
individuals (paired samples) or from different subjects (unpaired
samples).

## Usage

``` r
# S4 method for class 'MirnaExperiment'
mirnaDE(object, onlySignificant = TRUE, param = FALSE, returnObject = FALSE)

# S4 method for class 'MirnaExperiment'
geneDE(object, onlySignificant = TRUE, param = FALSE, returnObject = FALSE)

# S4 method for class 'MirnaExperiment'
significantMirnas(object)

# S4 method for class 'MirnaExperiment'
significantGenes(object)

# S4 method for class 'MirnaExperiment'
pairedSamples(object)

# S4 method for class 'MirnaExperiment'
mirnaTargets(object)

# S4 method for class 'MirnaExperiment'
integration(object, param = FALSE)

# S4 method for class 'MirnaExperiment'
show(object)
```

## Arguments

- object:

  An object of class `MirnaExperiment`

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

- the
  [`mirnaDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md)
  function returns a `data.frame` object with the results of miRNA
  differential expression analysis.

&nbsp;

- the
  [`geneDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md)
  function returns a `data.frame` object with the results of gene
  differential expression analysis.

&nbsp;

- the
  [`significantMirnas()`](https://jacopo-ronchi.github.io/MIRit/reference/significantAccessors.md)
  function gives back a `character` vector with the names of
  differentially expressed miRNAs.

&nbsp;

- the
  [`significantGenes()`](https://jacopo-ronchi.github.io/MIRit/reference/significantAccessors.md)
  function gives back a `character` vector with the names of
  differentially expressed genes.

&nbsp;

- the
  [`pairedSamples()`](https://jacopo-ronchi.github.io/MIRit/reference/pairedSamples.md)
  function returns a `logical` object that specifies if the
  MirnaExperiment object has paired measurements or not.

&nbsp;

- the
  [`mirnaTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/mirnaTargets.md)
  function provides a `data.frame` object with the list of target genes
  for the differentially expressed miRNAs.

&nbsp;

- the
  [`integration()`](https://jacopo-ronchi.github.io/MIRit/reference/integration.md)
  function gives back a `data.frame` object with the results of the
  integrative miRNA-mRNA analysis.

## Functions

- `mirnaDE(MirnaExperiment)`: Access the results of miRNA differential
  expression

- `geneDE(MirnaExperiment)`: Access the results of gene differential
  expression

- `significantMirnas(MirnaExperiment)`: Access the names of
  differentially expressed miRNAs

- `significantGenes(MirnaExperiment)`: Access the names of
  differentially expressed genes

- `pairedSamples(MirnaExperiment)`: Check if the object derives from
  sample-matched data

- `mirnaTargets(MirnaExperiment)`: Extract the miRNA-targets
  interactions retrieved for the differentially expressed miRNAs

- `integration(MirnaExperiment)`: Access the results of the integrative
  miRNA-mRNA analysis

- `show(MirnaExperiment)`: Show method for objects of class
  MirnaExperiment

## Slots

- `ExperimentList`:

  An `ExperimentList` class object for each assay dataset

- `colData`:

  A `DataFrame` of all clinical/specimen data available across
  experiments

- `sampleMap`:

  A `DataFrame` of translatable identifiers of samples and participants

- `metadata`:

  Additional data describing the object

- `drops`:

  A metadata `list` of dropped information

- `mirnaDE`:

  A `list` object containing the results of miRNA differential
  expression

- `geneDE`:

  A `list` object containing the results of gene differential expression

- `pairedSamples`:

  A `logical` parameter that specifies whether miRNA and gene expression
  measurements derive from the same individuals (`TRUE`) or from
  different subjects (`FALSE`)

- `targets`:

  A `data.frame` object containing miRNA-target pairs. This slot is
  commonly populated by the
  [`getTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/getTargets.md)
  function

- `integration`:

  A `list` object containing the results of the integration analysis
  between miRNA and gene expression values. This slot is commonly
  populated by the
  [`mirnaIntegration()`](https://jacopo-ronchi.github.io/MIRit/reference/mirnaIntegration.md)
  function

## Note

To create a `MirnaExperiment` object, you can use the
[`MirnaExperiment()`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment.md)
constructor function, which allows to easily build and verify a valid
object starting from miRNA and gene expression matrices.

## ExperimentList

The
[`ExperimentList`](https://github.com/waldronlab/MultiAssayExperiment/reference/ExperimentList.html)
slot is designed to contain results from each experiment/assay. In this
case, it holds miRNA and gene expression matrices. It contains a
[`SimpleList-class`](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html).

## colData

The `colData` slot is a collection of primary specimen data valid across
all experiments. This slot is strictly of class
[`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
but arguments for the constructor function allow arguments to be of
class `data.frame` and subsequently coerced.

## sampleMap

The `sampleMap` contains a
[`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
of translatable identifiers of samples and participants or biological
units. The standard column names of the sampleMap are "assay",
"primary", and "colname". Note that the "assay" column is a factor
corresponding to the names of each experiment in the
[`ExperimentList`](https://github.com/waldronlab/MultiAssayExperiment/reference/ExperimentList.html).
In the case where these names do not match between the sampleMap and the
experiments, the documented experiments in the `sampleMap` take
precedence and experiments are dropped by the harmonization procedure.
The constructor function will generate a `sampleMap` in the case where
it is not provided and this method may be a 'safer' alternative for
creating the `MultiAssayExperiment` (so long as the rownames are
identical in the `colData`, if provided). An empty `sampleMap` may
produce empty experiments if the levels of the "assay" factor in the
`sampleMap` do not match the names in the
[`ExperimentList`](https://github.com/waldronlab/MultiAssayExperiment/reference/ExperimentList.html).

## mirnaDE and geneDE

`mirnaDE` and `geneDE` consist of two `list` objects storing information
regarding miRNA and gene differential expression, including:

- `data`, which contains differential expression results in a
  `data.frame` with five columns:

  - `ID`: indicates the name of the miRNA/gene;

  - `logFC`: indicates the fold change of each feature in logarithmic
    scale;

  - `AveExpr`: represents the average expression of each miRNA/gene;

  - `P.Value`: indicates the resulting p-value;

  - `adj.P.Val`: contains the p-values adjusted for multiple testing.

- `significant`, which is a `character` vector containing the names of
  significantly differentially expressed miRNAs/genes that passed the
  thresholds;

- `method`, which specifies the procedure used to determine
  differentially expressed miRNAs/gens (eg. "limma-voom", "edgeR",
  "DESeq2", "limma");

- `group`, which is the column name of the variable (in colData) used
  for differential expression analysis;

- `contrast`, which represents the groups that are compared during
  differential expression analysis (e.g. 'disease-healthy');

- `design`, which outlines the R `formula` used for fitting the model.
  It includes the variable of interest (`group`) together with eventual
  covariates (e.g. '~ 0 + disease + sex');

- `pCutoff`, which indicates the p-value cutoff used for DE analysis;

- `pAdjustment`, the approach used for multiple testing correction;

- `logFC`, which states the log2 Fold Change cutoff used for DE
  analysis;

- `deObject`, an object deriving from limma/edgeR/DESeq2, that holds
  additional information regarding data processing.

MiRNA differential expression results can be accessed through the
[`mirnaDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md)
function, for additional details see
[`?mirnaDE`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md).
Similarly, gene differential expression results can be accessed through
the
[`geneDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md)
function, for additional details see
[`?geneDE`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md).

## pairedSamples

As already mentioned, `pairedSamples` must be `TRUE` when miRNA and gene
expression derive from the same subjects, while it must `FALSE` if this
is not the case.

## targets

`targets` is a `data.frame` with miRNA-target interactions, as retrieved
by
[`getTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/getTargets.md)
function.

## integration

Lastly, `integration` slot contains a `list` object that stores the
results and the options used for performing the integrative miRNA-gene
analysis. In particular, `integration` contains:

- `data`, which is a `data.frame` object with the results of the
  integrative analysis;

- `method`, which specifies the procedure used to perform the
  integrative analysis;

- `pCutoff`, which indicates the p-value cutoff used for the analysis;

- `pAdjustment`, the approach used for multiple testing correction.

Moreover, `data` differs on the basis of the integration strategy used.
For the one-sided association test integration, and for integration
based on rotation gene set tests, this `data.frame` has seven columns:

- `microRNA`: the miRNA ID;

- `mirna.direction`: the fold change direction of the DE-miRNA (`Up` or
  `Down`);

- `gene.direction`: the fold change direction of target genes (`Up` or
  `Down`);

- `DE`: represents the number of differentially expressed targets;

- `targets`: represents the total number of targets for this miRNA;

- `P.Val`: indicates the resulting p-value;

- `adj.P.Val`: contains test p-values corrected for multiple testing;

- `DE.targets`: contains the list of differentially expressed targets
  whose expression is negatively associated with miRNA expression.

Instead, when a correlation analysis is performed, `data` has six
columns:

- `microRNA`: the miRNA ID;

- `Target`: the correlated target gene;

- `microRNA.Direction`: the fold change direction of the DE-miRNA;

- `Corr.Coeff`: the value of the correlation coefficient used;

- `Corr.P.Value`: the p-value resulting from the correlation analysis;

- `Corr.Adjusted.P.Val`: contains the correlation p-values corrected for
  multiple testing.

To access the results of the integrative analysis, the `data` slot can
be accessed through the
[`integration()`](https://jacopo-ronchi.github.io/MIRit/reference/integration.md)
function.

## References

Marcel Ramos et al. Software For The Integration Of Multiomics
Experiments In Bioconductor. Cancer Research, 2017 November 1; 77(21);
e39-42. DOI:
[10.1158/0008-5472.CAN-17-0344](https://jacopo-ronchi.github.io/MIRit/reference/10.1158/0008-5472.CAN-17-0344)

Ronchi, J., & Foti, M. (2026). MIRit: An integrative R framework for the
identification of impaired miRNA–mRNA regulatory networks in complex
diseases. Bioinformatics Advances, vbag042.
<https://doi.org/10.1093/bioadv/vbag042>

## See also

See
[`MultiAssayExperiment-class`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment-class.html)
for additional information.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>
