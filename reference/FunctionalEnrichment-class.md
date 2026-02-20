# The `FunctionalEnrichment` class

This class introduces the possibility to store the results of functional
enrichment analyses such as over-representation analysis (ORA), gene set
enrichment analysis (GSEA), and competitive gene set test accounting for
inter-gene correlation (CAMERA). The different slots contained in this
class are used to store enrichment results generated through the
[`enrichGenes()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichGenes.md)
function.

## Usage

``` r
# S4 method for class 'FunctionalEnrichment'
enrichmentResults(object)

# S4 method for class 'FunctionalEnrichment'
enrichmentDatabase(object)

# S4 method for class 'FunctionalEnrichment'
enrichmentMethod(object)

# S4 method for class 'FunctionalEnrichment'
geneSet(object)

# S4 method for class 'FunctionalEnrichment'
enrichmentMetric(object)

# S4 method for class 'FunctionalEnrichment'
enrichedFeatures(object)

# S4 method for class 'FunctionalEnrichment'
show(object)
```

## Arguments

- object:

  An object of class `FunctionalEnrichment` containing enrichment
  results

## Value

- the
  [`enrichmentResults()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentResults.md)
  function returns a `data.frame` object with the results of the
  functional enrichment analysis.

&nbsp;

- the
  [`enrichmentDatabase()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentDatabase.md)
  function returns a `character` that specifies the database that
  provided the gene-set collection.

&nbsp;

- the
  [`enrichmentMethod()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentMethod.md)
  function returns a `character` with the enrichment method used for the
  analysis.

&nbsp;

- the
  [`geneSet()`](https://jacopo-ronchi.github.io/MIRit/reference/geneSet.md)
  function returns a `list` with the gene-set collections retrieved.

&nbsp;

- the
  [`enrichmentMetric()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentMetric.md)
  function returns a `numeric` object with ranked gene statistic used
  for GSEA.

&nbsp;

- the
  [`enrichmentMetric()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentMetric.md)
  function returns a `character` object with the names of the ranked
  features used for GSEA.

## Functions

- `enrichmentResults(FunctionalEnrichment)`: Access the `data` slot to
  take a closer look at all the enriched terms of an enrichment analysis

- `enrichmentDatabase(FunctionalEnrichment)`: See the database used for
  the functional enrichment

- `enrichmentMethod(FunctionalEnrichment)`: Visualize the approach used
  for the functional enrichment analysis

- `geneSet(FunctionalEnrichment)`: Access the `geneSet` slot to see the
  collection of gene sets used for GSEA

- `enrichmentMetric(FunctionalEnrichment)`: View the ranking metric used
  for GSEA

- `enrichedFeatures(FunctionalEnrichment)`: View the names of the
  pre-ranked features used for GSEA

- `show(FunctionalEnrichment)`: Show method for objects of class
  FunctionalEnrichment

## Slots

- `data`:

  A `data.frame` object holding the output of enrichment analysis

- `method`:

  The method used to perform functional enrichment analysis (e.g.
  `Gene Set Enrichment Analysis (GSEA)`)

- `organism`:

  The name of the organism under consideration (e.g. `Homo sapiens`)

- `database`:

  The name of the database used for the enrichment analysis (e.g.
  `KEGG`)

- `pCutoff`:

  A `numeric` value defining the threshold used for statistical
  significance in the enrichment analysis (e.g. `0.05`)

- `pAdjustment`:

  A `character` indicating the method used to correct p-values for
  multiple testing (e.g. `fdr`)

- `features`:

  A `character` vector containing the list of features used for the
  enrichment

- `statistic`:

  A `numeric` vector containing the statistic used to run GSEA. This
  parameter is empty for ORA and CAMERA

- `universe`:

  The background universe of features. Typically, this is equal to the
  complete list of features assayed. This slot is NULL for GSEA

- `geneSet`:

  The gene set used for the functional enrichment analysis. It is a
  `list` object where each element contains the list of genes belonging
  to a specific pathway.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>
