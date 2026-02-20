# Perform differential expression analysis

`performMirnaDE()` and `performGeneDE()` are two functions provided by
MIRit to conduct miRNA and gene differential expression analysis,
respectively. In particular, these functions allow the user to compute
differential expression through different methods, namely `edgeR`
(Quasi-Likelihood framework), `DESeq2`, `limma-voom` and `limma`. Data
deriving from NGS experiments and microarray technology are all suitable
for these functions. For precise indications about how to use these
functions, please refer to the *details* section.

## Usage

``` r
performMirnaDE(
  mirnaObj,
  group,
  contrast,
  design,
  method = "edgeR",
  logFC = 0,
  pCutoff = 0.05,
  pAdjustment = "fdr",
  filterByExpr.args = list(),
  calcNormFactors.args = list(),
  estimateDisp.args = list(robust = TRUE),
  glmQLFit.args = list(),
  glmQLFTest.args = list(),
  DESeq.args = list(),
  useVoomWithQualityWeights = TRUE,
  voom.args = list(),
  lmFit.args = list(),
  eBayes.args = list(),
  useArrayWeights = TRUE,
  useWsva = FALSE,
  wsva.args = list(),
  arrayWeights.args = list(),
  useDuplicateCorrelation = FALSE,
  correlationBlockVariable = NULL,
  duplicateCorrelation.args = list()
)

performGeneDE(
  mirnaObj,
  group,
  contrast,
  design,
  method = "edgeR",
  logFC = 0,
  pCutoff = 0.05,
  pAdjustment = "fdr",
  filterByExpr.args = list(),
  calcNormFactors.args = list(),
  estimateDisp.args = list(robust = TRUE),
  glmQLFit.args = list(),
  glmQLFTest.args = list(),
  DESeq.args = list(),
  useVoomWithQualityWeights = TRUE,
  voom.args = list(),
  lmFit.args = list(),
  eBayes.args = list(),
  useArrayWeights = TRUE,
  useWsva = FALSE,
  wsva.args = list(),
  arrayWeights.args = list(),
  useDuplicateCorrelation = FALSE,
  correlationBlockVariable = NULL,
  duplicateCorrelation.args = list()
)
```

## Arguments

- mirnaObj:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

- group:

  The variable of interest for differential expression analysis. It must
  be the column name of a variable present in the metadata (colData) of
  a
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object. See the *details* section for additional information

- contrast:

  A `character` object that specifies the groups to be compared during
  differential expression analysis, separated by a dash (e.g.
  'disease-healthy'). Note that reference group must be the last one,
  for additional information see the *details* section

- design:

  An R `formula` that indicates the model to fit. It must include the
  variable of interest (`group`) together with eventual covariates (e.g.
  '~ 0 + disease + sex'). Please note that `group` variable must be the
  first one. See the *details* section for additional information

- method:

  The statistical package used to compute differential expression. For
  NGS experiments, it must be one of `edgeR` (default), `DESeq2`, and
  `voom` (for limma-voom). Instead, for microarray data, only `limma`
  can be used

- logFC:

  The minimum log2 fold change required to consider a gene as
  differentially expressed. Optional, default is 0

- pCutoff:

  The adjusted p-value cutoff to use for statistical significance. The
  default value is `0.05`

- pAdjustment:

  The p-value correction method for multiple testing. It must be one of:
  `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
  `bonferroni`, `BY`

- filterByExpr.args:

  A `list` object containing additional arguments passed to
  [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
  function. It is used when `method` is set to `edgeR` or `voom`

- calcNormFactors.args:

  A `list` object containing additional arguments passed to
  [`edgeR::calcNormFactors()`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)
  function. It is used when `method` is set to `edgeR` or `voom`

- estimateDisp.args:

  A `list` object containing additional arguments passed to
  [`edgeR::estimateDisp()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html)
  function. It is used when `method` is set to `edgeR`. Default is
  `list(robust = TRUE)` to use the robust parameter

- glmQLFit.args:

  A `list` object containing additional arguments passed to
  [`edgeR::glmQLFit()`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html)
  function. It is used when `method` is set to `edgeR`

- glmQLFTest.args:

  A `list` object containing additional arguments passed to
  [`edgeR::glmQLFTest()`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html)
  function. It is used when `method` is set to `edgeR`

- DESeq.args:

  A `list` object containing additional arguments passed to
  [`DESeq2::DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html)
  function. It is used when `method` is set to `DESeq`

- useVoomWithQualityWeights:

  Logical, whether to use the
  [`limma::voomWithQualityWeights()`](https://rdrr.io/pkg/limma/man/voomWithQualityWeights.html)
  function or just the
  [`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html) function.
  It is used when `method` is set to `voom`. Default is TRUE

- voom.args:

  A `list` object containing additional arguments passed to
  [`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html) function or
  [`limma::voomWithQualityWeights()`](https://rdrr.io/pkg/limma/man/voomWithQualityWeights.html)
  function. It is used when `method` is set to `voom`

- lmFit.args:

  A `list` object containing additional arguments passed to
  [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html) function.
  It is used when `method` is set to `voom` or `limma`

- eBayes.args:

  A `list` object containing additional arguments passed to
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html)
  function. It is used when `method` is set to `voom` or `limma`

- useArrayWeights:

  Logical, whether to use the
  [`limma::arrayWeights()`](https://rdrr.io/pkg/limma/man/arrayWeights.html)
  function or not. It is used when `method` is set to `limma`. Default
  is TRUE

- useWsva:

  Logical, whether to use the
  [`limma::wsva()`](https://rdrr.io/pkg/limma/man/wsva.html) function or
  not. It is used when `method` is set to `limma`. Default is FALSE

- wsva.args:

  A `list` object containing additional arguments passed to
  [`limma::wsva()`](https://rdrr.io/pkg/limma/man/wsva.html) function.
  It is used when `method` is set to `limma`

- arrayWeights.args:

  A `list` object containing additional arguments passed to
  [`limma::arrayWeights()`](https://rdrr.io/pkg/limma/man/arrayWeights.html)
  function. It is used when `method` is set to `limma`

- useDuplicateCorrelation:

  Logical, whether to use the
  [`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html)
  function or not. It is used when `method` is set to `limma`. Default
  is FALSE

- correlationBlockVariable:

  It is the blocking variable to use for
  [`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html).
  Default is NULL

- duplicateCorrelation.args:

  A `list` object containing additional arguments passed to
  [`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html)
  function. It is used when `method` is set to `limma`

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

When performing differential expression for NGS experiments, count
matrices are detected and `method` parameter must be one of `edgeR`,
`DESeq2`, and `voom`. On the other hand, when dealing with microarray
studies, only `limma` can be used.

To calculate differential expression, MIRit must be informed about the
variable of interest and the desired contrast. In particular, the
`group` parameter must be the name of a variable present in the metadata
(colData) of a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object, which specifies the variable used to compute differential
expression analysis, between the groups indicated in `contrast`.
Specifically, `contrast` must be a character vector that defines the
levels to compare separated by a dash. For example, if we have a
variable named 'condition', with two levels, namely 'disease' and
'healthy', we can identify differentially expressed genes in 'disease'
samples compared to 'healthy' subjects by specifying:
`group = 'condition'` and `contrast = 'disease-healthy'`. Furthermore,
the user needs to specify the model to fit expression values. To do so,
the user has to state the model formula in the `design` parameter.
Please note that for a correct inner working of these functions, the
`group` variable of interest must be the *first* variable in model
formula. Moreover, the user can include in the design any other sources
of variation by specifying covariates that will be taken into account.
For instance, if we want to compare 'disease' subjects against 'healthy'
individuals, without the influence of sex differences, we may specify
`design = ~ condition + sex`, where 'sex' is also a variable present in
the metadata (colData) of `mirnaObj`.

Notably, for all the methods available, the user can supply additional
arguments to the functions implemented in `edgeR`, `DESeq2` and `limma`.
Therefore, the user has finer control over how the differential
expression analysis is performed. In this regard, for microarray
studies, the user may opt to include weighted surrogate variable
analysis (WSVA) to correct for unknown sources of variation
(`useWsva = TRUE`). Moreover, for microarray data, the `arrayWeights()`
function in `limma` can be used to assess differential expression with
respect to array qualities. Also, the `duplicateCorrelation()` function
in `limma` may be included in the pipeline in order to block the effect
of correlated samples. To do this, the user must set
`useDuplicateCorrelation = TRUE`, and must specify the blocking variable
through the `correlationBlockVariable` parameter. Additionally, when
using `limma-voom`, the user may estimate voom transformation with or
without quality weights (by specifying
`useVoomWithQualityWeights = TRUE`).

## Functions

- `performMirnaDE()`: Perform differential expression analysis for
  miRNAs

- `performGeneDE()`: Perform differential expression analysis for genes

## References

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
“limma powers differential expression analyses for RNA-sequencing and
microarray studies.” Nucleic Acids Research, 43(7), e47.
<doi:10.1093/nar/gkv007>.

Law, CW, Chen, Y, Shi, W, and Smyth, GK (2014). "Voom: precision weights
unlock linear model analysis tools for RNA-seq read counts". Genome
Biology 15, R29

Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor
package for differential expression analysis of digital gene expression
data.” Bioinformatics, 26(1), 139-140.
<doi:10.1093/bioinformatics/btp616>.

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change
and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550.
<doi:10.1186/s13059-014-0550-8>.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# \donttest{
# load example MirnaExperiment object
obj <- loadExamples()

# perform miRNA DE with edgeR
obj <- performMirnaDE(obj,
    group = "disease", contrast = "PTC-NTH",
    design = ~ 0 + disease + patient, method = "edgeR"
)
#> Performing differential expression analysis with edgeR...
#> Differential expression analysis reported 61 significant miRNAs with p < 0.05 (correction: fdr). You can use the 'mirnaDE()' function to access results.

# perform miRNA DE with DESeq2
obj <- performMirnaDE(obj,
    group = "disease", contrast = "PTC-NTH",
    design = ~ 0 + disease + patient, method = "DESeq2"
)
#> Performing differential expression analysis with DESeq2...
#> Warning: some variables in design formula are characters, converting to factors
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> Differential expression analysis reported 87 significant miRNAs with p < 0.05 (correction: fdr). You can use the 'mirnaDE()' function to access results.

# perform miRNA DE with limma-voom
obj <- performMirnaDE(obj,
    group = "disease", contrast = "PTC-NTH",
    design = ~ 0 + disease + patient, method = "voom"
)
#> Performing differential expression analysis with voom...
#> Differential expression analysis reported 54 significant miRNAs with p < 0.05 (correction: fdr). You can use the 'mirnaDE()' function to access results.
# }
```
