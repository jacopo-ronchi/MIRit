# Create a ridgeplot to display the results of GSEA analysis

This function creates a ridgeplot that is useful for showing the results
of GSEA analyses. The output of this function is a plot where enriched
terms/pathways found with
[`enrichGenes()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichGenes.md)
function are visualized on the basis of the ranking metric used for the
analysis. The resulting areas represent the density of signed p-values,
log2 fold changes, or log.p-values belonging to genes annotated to that
category.

## Usage

``` r
gseaRidgeplot(
  enrichment,
  showTerms = 10,
  showTermsParam = "padj",
  title = NULL
)
```

## Arguments

- enrichment:

  An object of class
  [`FunctionalEnrichment`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  containing enrichment results

- showTerms:

  It is the number of terms to be shown, based on the order determined
  by the parameter `showTermsParam`; or, alternatively, a character
  vector indicating the terms to plot. Default is `10`

- showTermsParam:

  The order in which the top terms are selected as specified by the
  `showTerms` parameter. It must be one of `ratio`, `padj` (default),
  `pval` and `overlap`

- title:

  The title of the plot. Default is `NULL` not to include a plot title

## Value

An object of class `ggplot` containing the ridgeplot of GSEA results.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example FunctionalEnrichment object
obj <- loadExamples("FunctionalEnrichment")

# plot results
gseaRidgeplot(obj)
#> Picking joint bandwidth of 0.561

```
