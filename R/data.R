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
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63511}
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
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63511}
#' 
"mirnaCounts"
