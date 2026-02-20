# Integrate microRNA and gene expression

This function allows to identify microRNAs that are significantly
associated/correlated with their targets. The principle is that, since
the biological role of miRNAs is mainly to negatively regulate gene
expression post-transcriptionally, the expression of a microRNA should
be negatively correlated with the expression of its targets. To test
this assumption for matched-sample data, this function performs a
correlation analysis. On the other hand, for unpaired data, it offers
different one-sided association tests to estimate if targets of
down-regulated miRNAs are enriched in up-regulated genes and vice versa.
Additionally, for unpaired data, miRNA effects on target gene expression
can also be quantified through a fast approximation to rotation gene-set
testing ('fry' method). For correlation analyses, the default behavior
is to use Spearman's correlation analysis, whereas for association tests
the default option makes use of a one-sided Boschloo's exact test. See
the *details* section for further information.

## Usage

``` r
mirnaIntegration(
  mirnaObj,
  test = "auto",
  pCutoff = 0.05,
  pAdjustment = "fdr",
  corMethod = "spearman",
  corCutoff = 0.5,
  partial = FALSE,
  partialCovs = NULL,
  associationMethod = "boschloo",
  nuisanceParam = 100,
  BPPARAM = bpparam()
)
```

## Arguments

- mirnaObj:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

- test:

  The statistical test to evaluate the association between miRNAs and
  genes. It must be one of `auto` (default), to automatically determine
  the appropriate statistical test; `correlation`, to perform a
  correlation analysis; `association`, to perform a one-sided
  association test; `fry` to perform the integrative analysis through
  rotation gene-set testing

- pCutoff:

  The adjusted p-value cutoff to use for statistical significance. The
  default value is `0.05`. When a lot of interactions are considered, a
  p-value cutoff after multiple testing correction could result
  excessively restrictive. In such cases, it is wise to just consider a
  threshold on the correlation strength and ignore p-values by setting
  `pCutoff = 1`.

- pAdjustment:

  The p-value correction method for multiple testing. It must be one of:
  `fdr` (default), `BH`, `none`, `holm`, `hochberg`, `hommel`,
  `bonferroni`, `BY`

- corMethod:

  The correlation method to be used for correlation analysis. It must be
  one of: `spearman` (default), `pearson`, `kendall`. See the *details*
  section for further information

- corCutoff:

  The minimum (negative) value of correlation coefficient to consider
  meaningful a miRNA-target relationship. Default is `0.5`

- partial:

  Logical, whether a partial correlation analysis should be performed.
  Default is `FALSE`. See the **details section** for further
  information

- partialCovs:

  Additional covariates to be considered in partial correlation
  analysis. This parameter is only considered when `TRUE`. It is an
  optional parameter that allows to include other covariates in the
  analysis in addition to `group`

- associationMethod:

  The statistical test used for evaluating the association between
  miRNAs and their targets for unpaired data. It must be one of
  `boschloo` (default), to perform a one-sided Boschloo's exact test;
  `fisher-midp`, to compute a one-sided Fisher's exact test with
  Lancaster's mid-p correction; `fisher`, to perform a one-sided
  Fisher's exact test

- nuisanceParam:

  The number of nuisance parameter values considered for p-value
  calculation in `boschloo` method. The higher this value, the better
  the p-value estimation accuracy. Default is 100

- BPPARAM:

  The desired parallel computing behavior. This parameter defaults to
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html),
  but this can be edited. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for information on parallel computing in R

## Value

A
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object containing integration results. To access these results, the user
can make use of the
[`integration()`](https://jacopo-ronchi.github.io/MIRit/reference/integration.md)
function. For additional details on how to interpret the results of
miRNA-gene integrative analysis, please see
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md).

## Details

As already pointed out, if miRNA and gene expression data derive from
the same samples, a correlation analysis is used. For evaluating these
relationships, the default method used is Spearman's correlation
coefficient, as:

- it does not need normally distributed data;

- it does not assume linearity;

- it is much more resistant to outliers.

However, the user can also decide to use other correlation methods, such
as Pearson's and Kendall's correlation. Nevertheless, for NGS data it
may happen that a certain number of ties is present in the expression
values. This can be handled by `spearman` method as it computes a
tie-corrected version of Spearman's coefficients. However, another
correlation method that is suitable to perform rank correlation on tied
data is the Kendall's tau-b method, usable with `kendall`.

Regarding correlation direction, since miRNAs mainly act as negative
regulators, only negatively correlated miRNA-target pairs are evaluated,
and statistical significance is calculated through a one-tailed t-test.

Additionally, when enough observations are present, it is appropriate to
account for the group effect by performing a partial correlation
analysis. In particular, a partial correlation analysis evaluates the
strength and the direction of a relationship between two variables –
miRNA and gene expression in our case – while accounting for the effect
of other factors. In integrative miRNA-mRNA analyses, the group effect
considered for differential expression analysis may lead to the
identification of several spurious correlated pairs, which result
anti-correlated simply because they are dysregulated in opposing
directions (upregulated miRNA and downregulated gene). This phenomenon,
known as Simpson's paradox, may therefore inflate false positive
relationships. By accounting for the group variable using partial
correlation analysis, the association between miRNA and gene expression
is evaluated within each group, thereby leading to reliable
identification of influential miRNAs. To perform such analysis, the
`partial` argument must be set to `TRUE`. Furthermore, the effect of
other covariates can be considered by passing a `character` vector with
the names of variables to account for to the `partialCovs` parameter.
However, partial correlation analyses are only effective when a
medium-large number of samples are available in each group. Our
simulations show that partial correlation outperforms standard
correlation when there are at least 20–30 samples for each condition,
this is way the default is set to `partial = FALSE`. Furthermore, for
batch effects that individually affect either miRNA or gene expression
matrices, the only way is to remove them using the
[`batchCorrection()`](https://jacopo-ronchi.github.io/MIRit/reference/batchCorrection.md)
function implemented in MIRit.

Moreover, if gene expression data and miRNA expression data derive from
different samples (unpaired data), a correlation analysis can't be
performed. However, one-sided association tests can be applied in these
cases to evaluate if targets of down-regulated miRNAs are statistically
enriched in up-regulated genes, and, conversely, if targets of
up-regulated miRNAs are statistically enriched in down-regulated genes.
In this case, Fisher's exact test can be used to assess the statistical
significance of this inverse association. Moreover, Lancaster's mid-p
adjustment can be applied since it has been shown that it increases
statistical power while retaining Type I error rates. However, Fisher's
exact test is a conditional test that requires the sum of both rows and
columns of a contingency table to be fixed. Notably, this is not true
for genomic data because it is likely that different datasets may lead
to a different number of DEGs. Therefore, the default behavior in MIRit
is to use a variant of Barnard's exact test, named Boschloo's exact
test, that is suitable when group sizes of contingency tables are
variable. Moreover, it is possible to demonstrate that Boschloo's test
is uniformly more powerful compared to Fisher's exact test.

Finally, for unpaired data, the effect of DE-miRNAs on the expression of
target genes can be estimated through rotation gene-set tests. In
particular, a fast approximation to rotation gene-set testing called
`fry`, implemented in the `limma` package, can be used to statistically
quantify the influence of miRNAs on the expression changes of their
target genes.

To speed up the identification of anti-correlated or anti-associated
miRNA-target pairs, this function implements parallel computation via
[`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).
In this regard, the parallelization behavior can be specified via the
`BPPARAM` parameter.

## References

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
“limma powers differential expression analyses for RNA-sequencing and
microarray studies.” Nucleic Acids Research, 43(7), e47.
<doi:10.1093/nar/gkv007>.

Di Wu and others, ROAST: rotation gene set tests for complex microarray
experiments, Bioinformatics, Volume 26, Issue 17, September 2010, Pages
2176–2182, <https://doi.org/10.1093/bioinformatics/btq401>.

Routledge, R. D. (1994). Practicing Safe Statistics with the Mid-p. The
Canadian Journal of Statistics / La Revue Canadienne de Statistique,
22(1), 103–110, <https://doi.org/10.2307/3315826>.

Boschloo R.D. (1970). "Raised Conditional Level of Significance for the
2x2-table when Testing the Equality of Two Probabilities". Statistica
Neerlandica. 24: 1–35. <doi:10.1111/j.1467-9574.1970.tb00104.x>.

Simpson, E. H. (1951). The Interpretation of Interaction in Contingency
Tables. Journal of the Royal Statistical Society: Series B
(Methodological), 13(2), 238–241.
<https://doi.org/10.1111/j.2517-6161.1951.tb00088.x>

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

# perform integration analysis with default settings
obj <- mirnaIntegration(obj)
#> Since data derive from paired samples, a correlation test will be used.
#> Performing Spearman's correlation analysis...
#> A statistically significant correlation between 215 miRNA-target pairs was found!
```
