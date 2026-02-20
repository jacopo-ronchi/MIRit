# The `IntegrativePathwayAnalysis` class

This class stores the output of integrative multi-omic pathway analyses.
In particular, the slots of this class are suitable to represent the
results of topologically-aware integrative pathway analysis (TAIPA)
returned from the
[`topologicalAnalysis()`](https://jacopo-ronchi.github.io/MIRit/reference/topologicalAnalysis.md)
function.

## Usage

``` r
# S4 method for class 'IntegrativePathwayAnalysis'
integratedPathways(object)

# S4 method for class 'IntegrativePathwayAnalysis'
integrationDatabase(object)

# S4 method for class 'IntegrativePathwayAnalysis'
augmentedPathways(object)

# S4 method for class 'IntegrativePathwayAnalysis'
show(object)
```

## Arguments

- object:

  An object of class `IntegrativePathwayAnalysis` containing the results
  of a miRNA-mRNA pathway analysis

## Value

- the
  [`integratedPathways()`](https://jacopo-ronchi.github.io/MIRit/reference/integratedPathways.md)
  function returns a `data.frame` with the results of the integrative
  pathway analysis.

&nbsp;

- the
  [`integrationDatabase()`](https://jacopo-ronchi.github.io/MIRit/reference/integrationDatabase.md)
  function returns a `character` object with the database used for the
  analysis.

&nbsp;

- the
  [`augmentedPathways()`](https://jacopo-ronchi.github.io/MIRit/reference/augmentedPathways.md)
  function returns a `list` object containing the miRNA-augmented
  pathways considered for TAIPA.

## Details

### Analysis results

The `data` slot of this class consists in a `data.frame` object with six
columns, namely:

- `pathway`, which indicates the name of the biological network;

- `coverage`, which specifies the fraction of nodes with expression
  measurement available;

- `score`, which expresses the score of each individual pathway;

- `normalized.score`, which indicates the pathway scores after
  standardizing the values for the null distribution computed through
  permutations;

- `P.Val`, the resulting p-value of each pathway;

- `adj.P.Val`, the p-value adjusted for multiple testing.

### Organisms and databases

The `organism` and `database` slots specify the organism in study and
the database used for retrieving biological interactions, respectively.
In particular, the
[`topologicalAnalysis()`](https://jacopo-ronchi.github.io/MIRit/reference/topologicalAnalysis.md)
function supports `KEGG`, `WikiPathways`, and `Reactome` databases.
Regarding organisms, the
[`supportedOrganisms()`](https://jacopo-ronchi.github.io/MIRit/reference/supportedOrganisms.md)
function can be used to retrieve the available species for each
database.

### Statistical significance of the permutation test

`pCutoff` and `pAdjustment` slots refer to the cutoff used for the
analysis. `pCutoff` is the threshold used for defining statistically
significant pathways, whereas `pAdjustment` refers to the multiple
testing correction method used.

Furthermore, since the statistical significance of each pathway is
defined on the basis of a permutation test, the number of permutations
is also specified in the `nPerm` slot.

### Augmented pathways

The `pathways` slot contains a `list` with weighted `graph` objects,
each representing a biological pathway. These networks have been
enlarged by adding the observed miRNA-mRNA interactions. Each network
has been processed so that the weight of each edge is +1 for activation
interactions, and -1 for repression interactions, such as those
occurring between miRNAs and mRNAs.

### Differential expression results for both miRNAs and genes

The expression variation of all miRNAs and genes measured in the study
is stored in the `expression` slot. In particular, this slot consists of
a `data.frame` object with different information, including log2 fold
changes, node weights and p-values.

### Minimum percentage of measured features

The `minPc` slot indicates the minimum percentage of miRNAs/mRNAs above
which pathways have been considered for the integrative analysis. This
is needed because often, when differential expression analysis is
performed, lowly expressed features are removed. Therefore, some
pathways might result significantly affected even if only 1% of nodes is
perturbed.

## Functions

- `integratedPathways(IntegrativePathwayAnalysis)`: Access the results
  of integrative miRNA-mRNA pathway analysis

- `integrationDatabase(IntegrativePathwayAnalysis)`: View the database
  used for the integrative pathway analysis

- `augmentedPathways(IntegrativePathwayAnalysis)`: Extract the list of
  biological networks augmented with miRNA-mRNA interactions

- `show(IntegrativePathwayAnalysis)`: Show method for objects of class
  IntegrativePathwayAnalysis

## Slots

- `data`:

  A `data.frame` object that contains the results of the integrative
  pathway analysis. See the *details* section for further details

- `method`:

  The method used for the analysis

- `organism`:

  The name of the organism under consideration (e.g. `Homo sapiens`)

- `database`:

  The name of the database used for retrieving biological pathways (e.g.
  `KEGG`)

- `pCutoff`:

  A `numeric` value defining the threshold used for statistical
  significance (e.g. `0.05`)

- `pAdjustment`:

  A `character` indicating the method used to correct p-values for
  multiple testing (e.g. `fdr`)

- `pathways`:

  A `list` of `graph` objects containing the biological networks
  retrieved from `database`, and augmented with miRNA-mRNA interactions

- `expression`:

  A `data.frame` object containing differential expression results for
  both miRNAs and genes

- `minPc`:

  The minimum percentage of measured features that a pathway must have
  for being considered in the analysis

- `nPerm`:

  The number of permutation used for assessing the statistical
  significance of each pathway

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>
