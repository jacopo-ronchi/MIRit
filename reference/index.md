# Package index

## The `MirnaExperiment` class

Methods and accessors for objects of class `MirnaExperiment`, the main
class in MIRit to work with miRNA and gene expression data.

- [`mirnaDE(`*`<MirnaExperiment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  [`geneDE(`*`<MirnaExperiment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  [`significantMirnas(`*`<MirnaExperiment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  [`significantGenes(`*`<MirnaExperiment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  [`pairedSamples(`*`<MirnaExperiment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  [`mirnaTargets(`*`<MirnaExperiment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  [`integration(`*`<MirnaExperiment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  [`show(`*`<MirnaExperiment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  : The 'MirnaExperiment' class
- [`MirnaExperiment()`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment.md)
  : The constructor function for MirnaExperiment
- [`pairedSamples()`](https://jacopo-ronchi.github.io/MIRit/reference/pairedSamples.md)
  : View the relationship between miRNA and gene samples

## Create example objects

Example datasets provided by MIRit for exploring the capabilities of the
software.

- [`mirnaCounts`](https://jacopo-ronchi.github.io/MIRit/reference/mirnaCounts.md)
  : Count matrix for microRNA expression in thyroid cancer
- [`geneCounts`](https://jacopo-ronchi.github.io/MIRit/reference/geneCounts.md)
  : Count matrix for gene expression in thyroid cancer
- [`loadExamples()`](https://jacopo-ronchi.github.io/MIRit/reference/loadExamples.md)
  : Load example MIRit objects

## Differential expression analysis

The functions used to perform miRNA and gene differential expression
analyses from start to finish.

- [`plotDimensions()`](https://jacopo-ronchi.github.io/MIRit/reference/plotDimensions.md)
  : Generate multidimensional scaling (MDS) plots to explore miRNA/gene
  expression distances

- [`performMirnaDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAnalysis.md)
  [`performGeneDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAnalysis.md)
  : Perform differential expression analysis

- [`addDifferentialExpression()`](https://jacopo-ronchi.github.io/MIRit/reference/addDifferentialExpression.md)
  : Manually add differential expression results to a MirnaExperiment
  object

- [`mirnaDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md)
  [`geneDE()`](https://jacopo-ronchi.github.io/MIRit/reference/deAccessors.md)
  :

  Extract differentially expressed miRNAs and genes from a
  `MirnaExperiment` object

- [`significantMirnas()`](https://jacopo-ronchi.github.io/MIRit/reference/significantAccessors.md)
  [`significantGenes()`](https://jacopo-ronchi.github.io/MIRit/reference/significantAccessors.md)
  : Get the IDs of statistically differentially expressed miRNAs/genes

- [`plotDE()`](https://jacopo-ronchi.github.io/MIRit/reference/plotDE.md)
  : Represent differentially expressed miRNAs/genes as boxplots,
  barplots or violinplots

- [`plotVolcano()`](https://jacopo-ronchi.github.io/MIRit/reference/plotVolcano.md)
  : Produce volcano plots to display miRNA/gene differential expression

## Functional enrichment

Methods and functions used to perform functional enrichment analysis of
genes.

- [`enrichmentResults(`*`<FunctionalEnrichment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  [`enrichmentDatabase(`*`<FunctionalEnrichment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  [`enrichmentMethod(`*`<FunctionalEnrichment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  [`geneSet(`*`<FunctionalEnrichment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  [`enrichmentMetric(`*`<FunctionalEnrichment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  [`enrichedFeatures(`*`<FunctionalEnrichment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  [`show(`*`<FunctionalEnrichment>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/FunctionalEnrichment-class.md)
  :

  The `FunctionalEnrichment` class

- [`enrichGenes()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichGenes.md)
  : Perform functional enrichment analysis of genes

- [`enrichTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichTargets.md)
  : Perform an enrichment analysis of integrated microRNA targets

- [`enrichedFeatures()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichedFeatures.md)
  : Extract the names of the pre-ranked features in a GSEA experiment

- [`enrichmentBarplot()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentBarplot.md)
  : Create a barplot for functional enrichment analysis

- [`enrichmentDatabase()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentDatabase.md)
  : Access the database used for functional enrichment analyses

- [`enrichmentDotplot()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentDotplot.md)
  : Create a dotplot for functional enrichment analysis

- [`enrichmentMethod()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentMethod.md)
  : Access the method used for functional enrichment analyses

- [`enrichmentMetric()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentMetric.md)
  : Extract the GSEA ranking metric used for functional enrichment
  analyses

- [`enrichmentResults()`](https://jacopo-ronchi.github.io/MIRit/reference/enrichmentResults.md)
  : Access the results of functional enrichment analyses

- [`geneSet()`](https://jacopo-ronchi.github.io/MIRit/reference/geneSet.md)
  : Extract the gene-sets used for functional enrichment analyses

- [`gseaPlot()`](https://jacopo-ronchi.github.io/MIRit/reference/gseaPlot.md)
  : Create a GSEA plot that displays the running enrichment score (ES)
  for a given pathway

- [`gseaRidgeplot()`](https://jacopo-ronchi.github.io/MIRit/reference/gseaRidgeplot.md)
  : Create a ridgeplot to display the results of GSEA analysis

- [`supportedOrganisms()`](https://jacopo-ronchi.github.io/MIRit/reference/supportedOrganisms.md)
  : Get the list of supported organisms for a given database

## Retrieval of miRNA targets

The functions for obtaining and visualizing the target genes of
differentially expressed miRNAs.

- [`getTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/getTargets.md)
  : Get microRNA targets
- [`mirnaTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/mirnaTargets.md)
  : Explore miRNA-target pairs
- [`setTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/setTargets.md)
  : Use custom miRNA-target interactions

## Disease-SNPs association

The functions used to associate differentially expressed miRNAs with
disease-related SNPs.

- [`searchDisease()`](https://jacopo-ronchi.github.io/MIRit/reference/searchDisease.md)
  : Search for disease EFO identifiers
- [`findMirnaSNPs()`](https://jacopo-ronchi.github.io/MIRit/reference/findMirnaSNPs.md)
  : Find disease-associated SNPs occurring at DE-miRNA loci
- [`mirVariantPlot()`](https://jacopo-ronchi.github.io/MIRit/reference/mirVariantPlot.md)
  : Create a trackplot to show the association between miRNAs and
  disease-SNPs
- [`getEvidence()`](https://jacopo-ronchi.github.io/MIRit/reference/getEvidence.md)
  : Get the scientific evidence for a particular disease-SNP association

## Integrate miRNA and gene expression

The functions for integrating miRNA and gene expression levels for both
paired and unpaired samples.

- [`batchCorrection()`](https://jacopo-ronchi.github.io/MIRit/reference/batchCorrection.md)
  : Correct for batch effects in miRNA and gene expression measurements
- [`mirnaIntegration()`](https://jacopo-ronchi.github.io/MIRit/reference/mirnaIntegration.md)
  : Integrate microRNA and gene expression
- [`integration()`](https://jacopo-ronchi.github.io/MIRit/reference/integration.md)
  : Explore the results of the integration analysis between miRNAs and
  genes
- [`plotCorrelation()`](https://jacopo-ronchi.github.io/MIRit/reference/plotCorrelation.md)
  : Plot correlation between miRNAs and genes within biological groups

## Topological pathway analysis

The functions needed to perform a comprehensive topology-aware
integrative pathway analysis (TAIPA) with MIRit.

- [`integratedPathways(`*`<IntegrativePathwayAnalysis>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/IntegrativePathwayAnalysis-class.md)
  [`integrationDatabase(`*`<IntegrativePathwayAnalysis>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/IntegrativePathwayAnalysis-class.md)
  [`augmentedPathways(`*`<IntegrativePathwayAnalysis>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/IntegrativePathwayAnalysis-class.md)
  [`show(`*`<IntegrativePathwayAnalysis>`*`)`](https://jacopo-ronchi.github.io/MIRit/reference/IntegrativePathwayAnalysis-class.md)
  :

  The `IntegrativePathwayAnalysis` class

- [`listPathways()`](https://jacopo-ronchi.github.io/MIRit/reference/listPathways.md)
  : List all the available biological pathways in KEGG, Reactome and
  WikiPathways

- [`preparePathways()`](https://jacopo-ronchi.github.io/MIRit/reference/preparePathways.md)
  : Prepare miRNA-augmented pathways for integrative miRNA-mRNA pathway
  analyses

- [`topologicalAnalysis()`](https://jacopo-ronchi.github.io/MIRit/reference/topologicalAnalysis.md)
  : Perform a topologically-aware integrative pathway analysis (TAIPA)

- [`integratedPathways()`](https://jacopo-ronchi.github.io/MIRit/reference/integratedPathways.md)
  : Access the results of integrative miRNA-mRNA pathway analyses

- [`visualizeNetwork()`](https://jacopo-ronchi.github.io/MIRit/reference/visualizeNetwork.md)
  : Visualize the relationships between miRNAs and genes in a biological
  pathway

- [`integrationDotplot()`](https://jacopo-ronchi.github.io/MIRit/reference/integrationDotplot.md)
  : Display integrated miRNA-mRNA augmented pathways in a dotplot

- [`augmentedPathways()`](https://jacopo-ronchi.github.io/MIRit/reference/augmentedPathways.md)
  : Access the miRNA-augmented pathways that were used during TAIPA

- [`integrationDatabase()`](https://jacopo-ronchi.github.io/MIRit/reference/integrationDatabase.md)
  : Extract the database used for integrative miRNA-mRNA pathway
  analyses
