#' Search for disease CUI identifiers
#'
#' This function allows to retrieve the UMLS Concept Unique Identifier (CUI)
#' of a particular disease. This identifier is then needed to use the function
#' [findMirnaSNPs()].
#'
#' @param diseaseName The name of a particular disease
#' (ex. `Creutzfeldt-Jakob disease`).
#'
#' @returns
#' A `data.frame` object containing CUI identifiers in the first column, and
#' full disease names in the second column.
#'
#' @examples
#' # search the CUI identifier of Creutzfeldt-Jakob disease
#' cui <- searchDisease("Creutzfeldt-Jakob disease")
#' cui
#'
#' # search the CUI identifier of Alzheimer's disease
#' cui <- searchDisease("Alzheimer's disease")
#' cui
#'
#' @note
#' To retrieve CUIs for specific diseases, this function makes use of the
#' `disgenet2r` package.
#'
#' @references
#' Piñero, J., Ramírez-Anguita, J. M., Saüch-Pitarch, J., Ronzano, F., Centeno,
#' E., Sanz, F., & Furlong, L. I. (2020). The DisGeNET knowledge platform for
#' disease genomics: 2019 update. Nucleic Acids Research, 48(D1), D845-D855.
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
searchDisease <- function(diseaseName) {

  ## check input
  if (!is.character(diseaseName) | length(diseaseName) != 1) {
    stop(paste("'diseaseName' must be a string of length 1",
                 "containing the name of the disease"), call. = FALSE)
  }

  ## set the DisGeNet api key
  disgenetApiKey <- "4cc61ac212f71bfefcdf375e31e879a135ffbb4d"

  ## search for the disease specified by the user
  disList <- disgenet2r::get_umls_from_vocabulary(disease = diseaseName,
                                                  vocabulary = "NAME",
                                                  api_key = disgenetApiKey)

  ## return disease list
  return(disList)

}





#' Find disease-associated SNPs occurring at DE-miRNA loci
#'
#' This function allows to identify genomic variants affecting
#' differentially expressed miRNA genes, that are associated with a
#' particular disease of interest. To do so, this function uses `gwasrapidd`
#' to retrieve SNPs-disease associations, and then retains only
#' SNPs that affect DE-miRNA genes.
#'
#' @details
#' SNPs occurring within miRNAs may have important effects on the biological
#' function of these transcripts. Indeed, a SNP present within a miRNA gene
#' might alter its expression or the spectrum of miRNA targets.
#'
#' To retrieve disease-SNPs, this function uses `gwasrapidd` package, which
#' directly queries the NHGRI-EBI Catalog of published genome-wide association
#' studies. After running this function, the user can use the [mirVariantPlot()]
#' function to produce a trackplot for visualizing the genomic presence of
#' SNPs within miRNA genes.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param diseaseEFO The EFO identifier of a disease of interest. This can be
#' identified with the [searchDisease()] function
#'
#' @returns
#' A `data.frame` containing details about disease-SNPs and the associated
#' differentially expressed miRNAs:
#' * `variant` contains SNP identifiers
#' * `miRNA` specifies the DE-miRNA gene present
#' * `chr` indicates the chromosome of SNPs
#' * `position` shows the SNP position
#' * `allele` displays possible alleles for this SNPs
#' * `distance` specifies the distance between SNPs and miRNAs
#' * `is_upstream` indicates wheter SNP is upstream of miRNA gene
#' * `is_downstream` indicates wheter SNP is downstream of miRNA gene
#' * `mirnaStrand` shows the strand
#' * `mirnaStartPosition` displays the start position of DE-miRNA gene
#' * `mirnaEndPosition` displays the end position of DE-miRNA gene
#'
#' @examples
#' # search disease
#' efo <- searchDisease("Alzheimer disease")
#' disId <- efo[1, 1]
#'
#' # retrieve associated SNPs
#' association <- findMirnaSNPs(obj, disId)
#' association
#'
#' @note
#' To retrieve disease-associated SNPs, this function makes use of the
#' `gwasrapidd` package.
#'
#' @references
#' Ramiro Magno, Ana-Teresa Maia, gwasrapidd: an R package to query, download
#' and wrangle GWAS catalog data, Bioinformatics, Volume 36, Issue 2, January
#' 2020, Pages 649–650, \url{https://doi.org/10.1093/bioinformatics/btz605}
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @importFrom rlang .data
#' @export
findMirnaSNPs <- function(mirnaObj,
                          diseaseEFO) {
  
  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (!is.character(diseaseEFO) | length(diseaseEFO) != 1) {
    stop(paste("'diseaseEFO' must be a string of length 1",
               "containing the EFO ID of the disease. You can obtain the",
               "disease identifier with 'searchDisease()'"), call. = FALSE)
  }
  
  ## inform the user about database querying
  message("Querying GWAS Catalog, this may take some time...")
  
  ## retrieve all SNPs related to the disease
  varDf <- gwasrapidd::get_variants(efo_trait = diseaseEFO)
  
  ## extract genomic context
  varDf <- varDf@genomic_contexts
  varDf <- varDf[varDf$chromosome_name %in% c(seq(22), "X", "Y"), seq(10)]
  varDf <- unique(varDf)
  
  ## only retain mapped genes
  varDf <- varDf[varDf$is_mapped_gene == "TRUE", ]
  
  ## find genomic details of differentially expressed miRNAs
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl",
                                 dataset = "hsapiens_gene_ensembl",
                                 mirror = "useast")
  deMirnas <- mirnaDE(mirnaObj)
  deMirnas <- deMirnas$ID
  deMirnas <- gsub("-5p|-3p", "", deMirnas)
  
  mirGenes <- biomaRt::select(ensembl,
                              keys = deMirnas,
                              keytype = "mirbase_id",
                              columns = c("hgnc_symbol",
                                          "gene_biotype",
                                          "chromosome_name",
                                          "strand",
                                          "start_position",
                                          "end_position"))
  
  ## intersect differentially expressed miRNAs with genes mapped with SNPs
  names(mirGenes)[which(names(mirGenes) == "hgnc_symbol")] <- "gene_name"
  mirMatches <- paste(mirGenes$gene_name, collapse = "|")
  resDf <- varDf[grepl(mirMatches, varDf$gene_name), ]
  resDf <- merge(resDf, mirGenes, by = c("gene_name", "chromosome_name"))
  
  ## add allele information
  snpData <- biomaRt::useEnsembl(biomart = "snps", dataset = "hsapiens_snp")
  snpInfo <- AnnotationDbi::select(snpData,
                                   keys = resDf$variant,
                                   keytype = "snp_filter",
                                   columns = c("refsnp_id",
                                               "allele"))
  colnames(snpInfo)[1] <- "variant_id"
  resDf <- merge(resDf, snpInfo, by = "variant_id")
  
  ## prepare resulting data.frame
  resDf <- resDf[, c("variant_id", "gene_name", "chromosome_name",
                     "chromosome_position", "allele", "distance", "is_upstream",
                     "is_downstream", "strand", "start_position",
                     "end_position")]
  colnames(resDf) <- c("variant", "miRNA", "chr", "position", "allele",
                       "distance", "is_upstream", "is_downstream",
                       "mirnaStrand", "mirnaStartPosition", "mirnaEndPosition")
  
  ## stop the code if no SNP is associated with DE-miRNAs
  if (nrow(resDf) == 0) {
    stop("No disease-associated SNPs occur at DE-miRNA loci!", call. = FALSE)
  }
  
  ## display a message with miRNA-SNPs association results
  message(paste("After the analysis,", nrow(resDf), "variants associated with",
                diseaseEFO, "were found within differentially",
                "expressed miRNA genes"))
  
  ## return SNPs association data.frame
  return(resDf)
  
}





#' Get the scientific evidence for a particular disease-SNP association
#'
#' This function returns the biomedical evidence that supports the association
#' between a particular SNP and a phenotypic trait.
#'
#' To retrieve disease-SNPs, this function uses DisGeNet, which is one of the
#' largest public databases containing information about human gene-disease
#' associations and human variants-disease associations. For variants, it
#' integrates information from multiple sources, such as GWAS catalogues and
#' biomedical literature. The `database` parameter specifies the database from
#' which the function obtains disease-associated SNPs. The options include:
#' * `CLINVAR`, to use ClinVar;
#' * `UNIPROT`, to use the Universal Protein Resource;
#' * `GWASCAT`, to use the NHGRI-EBI GWAS Catalog;
#' * `GWASDB`, to use GWAS Database GWASdb;
#' * `BEFREE`, to use text mining data, generated using BeFree System;
#' * `CURATED`, to use expert curated, human databases,
#' * `ALL`, to use all the above mentioned databases together (default).
#'
#' @param variant The SNP ID of a particular variant of interest
#' (e.g. 'rs394581')
#' @param diseaseUML The CUI identifier of a disease of interest. This can be
#' identified with the function [searchDisease()]
#' @param database The name of the database of disease-SNPs associations to
#' draw from. The default value is `ALL` to use all databases included in
#' DisGeNet. For a detailed list of all the possible values, see the details
#' section
#' @param score A numeric vector of length two containing the minimum value of
#' score and the maximum value of score. Default is `c(0, 1)`. The score is a
#' DisGeNet metric that estimates the strength of the association between
#' the variant and the disease
#'
#' @returns
#' A `data.frame` containing information about literature evidences for a
#' disease-SNP association.
#'
#' @examples
#' dis <- searchDisease("Alzheimer's disease")
#' dis
#'
#' association <- findMirnaSNPs(obj, diseaseUML = "C0002395")
#' evidence <- getEvidence("rs2075650", diseaseUML = "C0002395")
#'
#' @note
#' To retrieve evidences of disease-SNP association, this function makes use of
#' the `disgenet2r` package.
#'
#' @references
#' Piñero, J., Ramírez-Anguita, J. M., Saüch-Pitarch, J., Ronzano, F., Centeno,
#' E., Sanz, F., & Furlong, L. I. (2020). The DisGeNET knowledge platform for
#' disease genomics: 2019 update. Nucleic Acids Research, 48(D1), D845-D855.
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
getEvidence <- function(variant,
                        diseaseUML,
                        database = "ALL",
                        score = c(0, 1)) {

  ## check inputs
  if (!is.character(variant) | length(variant) != 1) {
    stop("'variant' should be a character with the ID of a SNP.", call. = FALSE)
  }
  if (!is.character(diseaseUML) | length(diseaseUML) != 1) {
    stop(paste("'diseaseUML' must be a string of length 1",
               "containing the CUI of the disease. You can obtain the",
               "disease identifier with 'searchDisease()'"), call. = FALSE)
  }
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("UNIPROT", "CLINVAR", "GWASCAT", "GWASDB",
                       "CURATED", "BEFREE", "ALL")) {
    stop(paste("'database' must be one of: 'UNIPROT', 'CLINVAR', 'GWASCAT',
               'GWASDB', 'CURATED', 'BEFREE', 'ALL'."), call. = FALSE)
  }
  if (!is.numeric(score) |
      length(score) != 2) {
    stop(paste("'score' must be a numeric vector of length 2 with the initial",
               "and final values of score. Default: 'c(0, 1)'"), call. = FALSE)
  }

  ## set the DisGeNet api key
  disgenetApiKey <- "4cc61ac212f71bfefcdf375e31e879a135ffbb4d"

  ## explore the evidences in literature for a particular variant
  suppressMessages(
    disEv <- disgenet2r::disease2evidence(disease = diseaseUML,
                                          type = "VDA",
                                          database = database,
                                          variant = variant,
                                          api_key = disgenetApiKey,
                                          score = score,
                                          warnings = FALSE)
  )

  ## return evidence data.frame
  if (!is.character(disEv)) {
    return(disgenet2r::extract(disEv))
  } else {
    warning("No evidence was found for this disease-associated SNP")
  }
}

