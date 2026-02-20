# Find disease-associated SNPs occurring at DE-miRNA loci

This function allows to identify disease-associated genomic variants
affecting differentially expressed miRNA genes or their host genes. To
do so, this function uses `gwasrapidd` to retrieve SNPs-disease
associations, and then retains only SNPs that affect DE-miRNA genes or
their relative host genes (for intronic miRNAs).

## Usage

``` r
findMirnaSNPs(mirnaObj, diseaseEFO)
```

## Arguments

- mirnaObj:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

- diseaseEFO:

  The EFO identifier of a disease of interest. This can be identified
  with the
  [`searchDisease()`](https://jacopo-ronchi.github.io/MIRit/reference/searchDisease.md)
  function

## Value

A `data.frame` containing details about disease-SNPs and the associated
differentially expressed miRNAs:

- `variant` contains SNP identifiers;

- `gene` defines the gene affected by a disease-SNP (it may be the miRNA
  gene itself or the host gene of an intronic miRNA);

- `miRNA.gene` specifies the DE-miRNA gene present;

- `miRNA.precursor` specifies the name of the miRNA precursor affected
  by disease-SNPs;

- `chr` indicates the chromosome of SNPs;

- `position` shows the SNP position;

- `allele` displays possible alleles for this SNPs;

- `distance` specifies the distance between SNPs and miRNAs;

- `is_upstream` indicates whether SNP is upstream of miRNA gene;

- `is_downstream` indicates whether SNP is downstream of miRNA gene;

- `mirnaStrand` shows the strand;

- `mirnaStartPosition` displays the start position of DE-miRNA gene;

- `mirnaEndPosition` displays the end position of DE-miRNA gene.

## Details

SNPs occurring within miRNAs may have important effects on the
biological function of these transcripts. Indeed, a SNP present within a
miRNA gene might alter its expression or the spectrum of miRNA targets.

To retrieve disease-SNPs, this function uses `gwasrapidd` package, which
directly queries the NHGRI-EBI Catalog of published genome-wide
association studies. After running this function, the user can use the
[`mirVariantPlot()`](https://jacopo-ronchi.github.io/MIRit/reference/mirVariantPlot.md)
function to produce a trackplot for visualizing the genomic location of
SNPs within miRNA genes.

## Note

To retrieve disease-associated SNPs, this function makes use of the
`gwasrapidd` package.

## References

Ramiro Magno, Ana-Teresa Maia, gwasrapidd: an R package to query,
download and wrangle GWAS catalog data, Bioinformatics, Volume 36, Issue
2, January 2020, Pages 649–650,
<https://doi.org/10.1093/bioinformatics/btz605>

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example MirnaExperiment object
obj <- loadExamples()

# \donttest{
# search disease
searchDisease("response to antidepressant")
#> Checking for cached EFO traits...
#> Downloading EFO traits, this may take some minutes...
#> Searching for disease: response to antidepressant
#> [1] "response to antidepressant" "response to anticonvulsant"
disId <- "response to antidepressant"

# retrieve associated SNPs
association <- findMirnaSNPs(obj, disId)
#> Querying GWAS Catalog, this may take some time...
#> Finding genomic information of differentially expressed miRNAs...
#> Ensembl site unresponsive, trying useast mirror
#> After the analysis, 1 variants associated with response to antidepressant were found within differentially expressed miRNA genes
# }
```
