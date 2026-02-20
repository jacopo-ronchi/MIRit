# Correct for batch effects in miRNA and gene expression measurements

This function allows to remove unwanted batch effects from miRNA and
gene expression matrices. In particular, this function fits a linear
model to miRNA/gene expression levels, and then removes the variability
caused by batch effects. Furthermore, a weighted surrogate variable
analysis (WSVA) can also be included to remove the effects due to
surrogate variables.

## Usage

``` r
batchCorrection(
  mirnaObj,
  assay,
  batch = NULL,
  batch2 = NULL,
  covariates = NULL,
  includeWsva = FALSE,
  n.sv = 1L,
  weight.by.sd = TRUE
)
```

## Arguments

- mirnaObj:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

- assay:

  The expression matrix to correct. It must be one of `genes` and
  `microRNA`

- batch:

  It must be the name of a variable present in the `colData` of a
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object (eg. "disease"), or, alternatively, it must be a
  `character`/`factor` object that defines batch memberships. See the
  **details** section for additional information

- batch2:

  It must be the name of a variable present in the `colData` of a
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object (eg. "disease"), or, alternatively, it must be a
  `character`/`factor` object that defines another series of batches
  that have additive effects to those specified in `batch`. See the
  **details** section for additional information

- covariates:

  Additional numeric covariates that we want to correct for. It must be
  a `character` vector containing the names of numeric variables present
  in the `colData` of a
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object (eg. `c("age", "RIN", "quantity")`), or, alternatively, it must
  be a simple `matrix` object. See the **details** section for
  additional information

- includeWsva:

  Logical, whether to correct for surrogate variables or not. Default is
  FALSE

- n.sv:

  The number of surrogate variables to estimate

- weight.by.sd:

  Logical, whether to specifically tune the surrogate variables to the
  more variable genes or not. Default is TRUE

## Value

A
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object containing batch effect-corrected expression matrices.

## Details

Batch effects consist in unwanted sources of technical variation that
confound expression variability and limit downstream analyses. Since the
reliability of biological conclusions of integrative miRNA-mRNA analyses
depends on the association between miRNA and gene expression levels, it
is pivotal to ensure that expression measurements are not affected by
technical variations. To account for effects that influence both gene
and miRNA expression (e.g. sex, smoking, hospitals), the most
appropriate way is to perform a partial correlation analysis by
including `partial = TRUE` when calling the
[`mirnaIntegration()`](https://jacopo-ronchi.github.io/MIRit/reference/mirnaIntegration.md)
function. Nevertheless, for effects that are impacting either gene or
miRNA expression matrix, the only option is to regress out their
influence using this function before using the
[`mirnaIntegration()`](https://jacopo-ronchi.github.io/MIRit/reference/mirnaIntegration.md)
function to perform a correlation analysis.

Usually, given a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object, the user should specify:

- the `assay` from which we want to remove batch effects (one between
  `genes` and `microRNA`);

- the `batch` variable, which is a variable that defines the different
  batches;

- the `batch2` variable, which can be included to correct for a second
  series of batches that have additive effects to those specified in
  `batch`;

- the `covariates` variables, which allows correction for one or more
  continuous numeric effects.

In particular, `batch` and `batch2` could be provided as the names of
covariates included in the `colData` of a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object. Alternatively, they can be `character`/`factor` objects that
declare batch memberships. Similarly, `covariates` can be supplied as a
vector containing the names of numeric variables listed in the `colData`
of
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
objects, or they can be provided as a simple `matrix`.

Additionally, the influence of unknown sources of technical variation
can be removed by including surrogate variables estimated through WSVA.
To do so, we can set `includeWsva` to TRUE, and then we can specify the
number of surrogate variables to use through the `n.sv` parameter.
Further, the surrogate variables can be tuned to the more variable genes
by setting `weight.by.sd` to TRUE.

Please note that we only recommend to remove batch effects directly from
expression measurements prior to correlation analysis. This function
can't be used to remove batch effects before differential expression
analysis, because for that purpose, it is better to include batch
variables in the linear model. In this way, we do not underestimate the
residual degrees of freedom, so that the calculated standard errors,
t-statistics and p-values are not overoptimistic.

## Note

To estimate surrogate variables and to remove batch effects from
expression data, MIRit uses the
[`limma::wsva()`](https://rdrr.io/pkg/limma/man/wsva.html) and
[`limma::removeBatchEffect()`](https://rdrr.io/pkg/limma/man/removeBatchEffect.html)
functions, respectively.

## References

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
“limma powers differential expression analyses for RNA-sequencing and
microarray studies.” Nucleic Acids Research, 43(7), e47.
<doi:10.1093/nar/gkv007>.

Ronchi, J., & Foti, M. (2026). MIRit: An integrative R framework for the
identification of impaired miRNA–mRNA regulatory networks in complex
diseases. Bioinformatics Advances, vbag042.
<https://doi.org/10.1093/bioadv/vbag042>

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example MirnaExperiment object
obj <- loadExamples()

# correct batch effects due to the patient from miRNA expression matrix
obj <- batchCorrection(obj, "microRNA", batch = "patient")
```
