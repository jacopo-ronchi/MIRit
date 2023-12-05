## documentation for 'geneCounts' dataset


#' Count matrix for gene expression in thyroid cancer
#'
#' This dataset contains a gene expression matrix resulting from an RNA-Seq
#' analysis of thyroid cancer. Specifically, these data originate from
#' Riesco-Eizaguirre et al (2015), where they sequenced 8 papillary thyroid
#' carcinomas (PTC) together with paired samples of normal thyroid tissue.
#' The same thing was done for microRNAs in order to investigate the effects
#' on target genes. Data included in this package have been obtained through
#' the Gene Expression Omnibus (GEO accession: GSE63511).
#'
#' @format ## `geneCounts`
#' A `matrix` object containing samples as columns and genes as rows.
#'
#' @returns A `matrix` object with 23710 rows and 16 columns.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63511}
#'
#' @references
#' Garcilaso Riesco-Eizaguirre et al., “The MiR-146b-3p/PAX8/NIS Regulatory
#' Circuit Modulates the Differentiation Phenotype and Function of Thyroid Cells
#' during Carcinogenesis,” Cancer Research 75, no. 19 (September 30, 2015):
#' 4119–30, \url{https://doi.org/10.1158/0008-5472.CAN-14-3547}.
#'
#' @usage data(geneCounts)
#'
"geneCounts"


#' Count matrix for microRNA expression in thyroid cancer
#'
#' This dataset contains a gene expression matrix resulting from a miRNA-Seq
#' analysis of thyroid cancer. Specifically, these data originate from
#' Riesco-Eizaguirre et al (2015), where they sequenced 8 papillary thyroid
#' carcinomas (PTC) together with paired samples of normal thyroid tissue.
#' The same thing was done for mRNAs in order to investigate the effects
#' on target genes. Data included in this package have been obtained through
#' the Gene Expression Omnibus (GEO accession: GSE63511).
#'
#' @format ## `mirnaCounts`
#' A `matrix` object containing samples as columns and microRNAs as rows.
#'
#' @returns A `matrix` object with 2576 rows and 16 columns.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63511}
#'
#' @references
#' Garcilaso Riesco-Eizaguirre et al., “The MiR-146b-3p/PAX8/NIS Regulatory
#' Circuit Modulates the Differentiation Phenotype and Function of Thyroid Cells
#' during Carcinogenesis,” Cancer Research 75, no. 19 (September 30, 2015):
#' 4119–30, \url{https://doi.org/10.1158/0008-5472.CAN-14-3547}.
#'
#' @usage data(mirnaCounts)
#'
"mirnaCounts"
