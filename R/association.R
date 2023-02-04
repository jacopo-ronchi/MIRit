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
#' [disgenet2r] package.
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
#' This function allows to identify genomic variants, positioned within
#' differentially expressed miRNA genes, that are associated with a
#' particular disease of interest. To do so, this function uses the DisGeNet
#' database to retrieve SNPs-disease associations, and then retains only
#' SNPs that affect DE-miRNA genes.
#'
#' @details
#' SNPs occurring within miRNAs may have important effects on the biological
#' function of these transcripts. Indeed, a SNP present within a miRNA might
#' alter the spectrum of targets regulated by a particular microRNA.
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
#' After running this function, the user can use the [mirVariantPlot()]
#' function to produce a trackplot for visualizing the genomic presence of
#' SNPs within miRNA genes.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param diseaseUML The CUI identifier of a disease of interest. This can be
#' identified with the [searchDisease()] function
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
#' A `data.frame` containing details about disease-SNPs and the associated
#' differentially expressed miRNAs.
#'
#' The `variant` column contains SNP identifiers, the `disease` column
#' represents the disease name, the `mirnaGene` columns specifies the DE-miRNA
#' gene present, and the `allele` column shows the possible alleles of the
#' variant. Then, there are the `chromosome` and `position` columns, that give
#' the genomic coordinates of the SNP, and the `genesPresent` column that
#' indicates the genes present in the genomic trait in which the SNP is
#' located. Further, there are the `mirnaStart`, `mirnaEnd` and `mirnaStrand`
#' columns that provide genomic details of the differentially expressed miRNAs
#' affected by disease-SNPs. Finally, there are the `variantDSI`, `variantDPI`,
#' `ei` and `score` columns which contain different association metrics
#' obtained from DisGeNet.
#'
#' @examples
#' # search disease
#' cuis <- searchDisease("Alzheimer's disease")
#' disId <- cuis[1, 1]
#'
#' # retrieve associated SNPs
#' association <- findMirnaSNPs(obj, disId)
#' association
#'
#' @note
#' To retrieve disease-associated SNPs, this function makes use of the
#' [disgenet2r] package.
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
findMirnaSNPs <- function(mirnaObj,
                          diseaseUML,
                          database = "ALL",
                          score = c(0, 1)) {

  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
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

  ## retrieve all SNPs related to the disease
  disVar <- disgenet2r::disease2variant(disease = diseaseUML,
                                        database = database,
                                        score = score,
                                        api_key = disgenetApiKey)

  ## extract results
  varDf <- disgenet2r::extract(disVar)

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

  ## integrate variants information and miRNA details into one data.frame
  names(mirGenes)[which(names(mirGenes) == "hgnc_symbol")] <- "gene_symbol"
  mirMatches <- paste(mirGenes$gene_symbol, collapse = "|")
  resDf <- varDf[grepl(mirMatches, varDf$gene_symbol), ]
  resDf$genes <- resDf$gene_symbol
  resDf <- dplyr::mutate(resDf, gene_symbol = sapply(strsplit(gene_symbol, ";"), function(z)
    paste0("\\b(", paste(z, collapse = "|"), ")\\b")))
  resDf <- fuzzyjoin::regex_left_join(mirGenes, resDf, by = "gene_symbol")
  resDf <- resDf[, c("variantid", "disease_name", "gene_symbol.x", "genes",
                     "source", "variant_dsi", "variant_dpi", "ei", "score",
                     "chromosome_name", "strand",
                     "start_position", "end_position")]
  resDf <- na.omit(resDf)
  colnames(resDf) <- c("variant", "disease", "mirnaGene", "genesPresent",
                       "source", "variantDSI", "variantDPI", "ei", "score",
                       "chr", "mirnaStrand",
                       "mirnaStart", "mirnaEnd")
  resDf$mirnaStrand <- gsub("-1", "-", resDf$mirnaStrand)
  resDf$mirnaStrand <- gsub("1", "+", resDf$mirnaStrand)

  ## stop the code if no SNP is associated with DE-miRNAs
  if (nrow(resDf) == 0) {
    stop("No disease-associated SNPs occur at DE-miRNA loci!", call. = FALSE)
  }

  ## retrieve SNPs location and information
  snpData <- biomaRt::useEnsembl(biomart = "snps", dataset = "hsapiens_snp")
  snpInfo <- AnnotationDbi::select(snpData,
                                   keys = resDf$variant,
                                   keytype = "snp_filter", columns = c("refsnp_id",
                                                                       "chr_name",
                                                                       "chrom_start",
                                                                       "allele"))

  ## retain only SNPs on standard chromosomes
  snpInfo <- snpInfo[snpInfo$chr_name %in% c(as.character(seq(22)), "X", "Y"), ]
  colnames(snpInfo)[1] <- "variant"

  ## add SNPs details and return the data.frame
  resDf <- dplyr::left_join(resDf, snpInfo, by = "variant")
  resDf <- resDf[, c("variant", "disease", "mirnaGene", "allele", "chr",
                     "chrom_start", "genesPresent", "mirnaStart", "mirnaEnd",
                     "mirnaStrand", "variantDSI", "variantDPI", "ei", "score")]
  colnames(resDf)[c(5, 6)] <- c("chromosome", "position")

  ## display a message with miRNA-SNPs association results
  message(paste("After the analysis,", nrow(resDf), "variants associated with",
                "the disease were found within differentially",
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
#' dis <- searchDisease("Alzheimer's disease)
#' dis
#'
#' association <- findMirnaSNPs(obj, diseaseUML = "C0002395")
#' evidence <- getEvidence("rs2075650", diseaseUML = "C0002395")
#'
#' @note
#' To retrieve evidences of disease-SNP association, this function makes use of
#' the [disgenet2r] package.
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

