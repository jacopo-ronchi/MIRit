#' Search for disease EFO identifiers
#'
#' This function allows to retrieve the Experimental Factor Ontology (EFO)
#' identifier of a particular disease. This ID is then needed to use the
#' function [findMirnaSNPs()].
#'
#' @param diseaseName The name of a particular disease
#' (ex. `Alzheimer disease`).
#'
#' @returns
#' A `character` object containing EFO identifiers.
#'
#' @examples
#' \donttest{
#' # search the EFO identifier of Alzheimer disease
#' searchDisease("Alzheimer disease")
#' }
#'
#' @note
#' To retrieve EFO IDs for specific diseases, this function makes use of the
#' `gwasrapidd` package.
#'
#' @references
#' Ramiro Magno, Ana-Teresa Maia, gwasrapidd: an R package to query, download
#' and wrangle GWAS catalog data, Bioinformatics, Volume 36, Issue 2, January
#' 2020, Pages 649–650, \url{https://doi.org/10.1093/bioinformatics/btz605}.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
searchDisease <- function(diseaseName) {
    ## check input
    if (!is.character(diseaseName) | length(diseaseName) != 1) {
        stop("'diseaseName' must be a string of length 1 ",
             "containing the name of the disease",
             call. = FALSE
        )
    }
    
    ## check that gwasrapidd is installed
    if (!requireNamespace("gwasrapidd", quietly = TRUE)) {
        stop("The gwasrapidd package is needed to use ",
             "this function. Install it through: ",
             paste("`install.packages('gwasrapidd')`.",
                   sep = ""
             ),
             call. = FALSE
        )
    }
    
    ## load cache
    bfc <- .get_cache()
    
    ## check if EFO traits are cached
    message("Checking for cached EFO traits...")
    cache <- BiocFileCache::bfcquery(bfc, "EFO_ids")
    if ("EFO_ids" %in% cache$rname) {
        ## load cached EFO ids
        message("Reading EFO traits from cache...")
        ids <- readRDS(BiocFileCache::bfcrpath(bfc, "EFO_ids")[1])
    } else {
        ## inform the user about retrieving IDs
        message("Downloading EFO traits, this may take some minutes...")
        
        ## retrieve disease EFO traits
        ids <- gwasrapidd::get_traits()
        ids <- ids@traits$trait
        
        ## save EFO traits to cache
        savepath <- BiocFileCache::bfcnew(bfc, rname = "EFO_ids", ext = ".RDS")
        saveRDS(ids, file = savepath)
    }
    
    ## return results
    message("Searching for disease: ", diseaseName)
    
    ## search for the disease specified by the user
    disList <- ids[agrep(diseaseName,
                         ids,
                         ignore.case = TRUE,
                         max.distance = 0.2
    )]
    
    ## return disease list
    return(disList)
}





#' Find disease-associated SNPs occurring at DE-miRNA loci
#'
#' This function allows to identify disease-associated genomic variants
#' affecting differentially expressed miRNA genes or their host genes.
#' To do so, this function uses `gwasrapidd` to retrieve SNPs-disease
#' associations, and then retains only SNPs that affect DE-miRNA genes or their
#' relative host genes (for intronic miRNAs).
#'
#' @details
#' SNPs occurring within miRNAs may have important effects on the biological
#' function of these transcripts. Indeed, a SNP present within a miRNA gene
#' might alter its expression or the spectrum of miRNA targets.
#'
#' To retrieve disease-SNPs, this function uses `gwasrapidd` package, which
#' directly queries the NHGRI-EBI Catalog of published genome-wide association
#' studies. After running this function, the user can use the [mirVariantPlot()]
#' function to produce a trackplot for visualizing the genomic location of
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
#'
#' * `variant` contains SNP identifiers;
#' * `gene` defines the gene affected by a disease-SNP (it may be the miRNA
#' gene itself or the host gene of an intronic miRNA);
#' * `miRNA.gene` specifies the DE-miRNA gene present;
#' * `miRNA.precursor` specifies the name of the miRNA precursor affected
#' by disease-SNPs;
#' * `chr` indicates the chromosome of SNPs;
#' * `position` shows the SNP position;
#' * `allele` displays possible alleles for this SNPs;
#' * `distance` specifies the distance between SNPs and miRNAs;
#' * `is_upstream` indicates whether SNP is upstream of miRNA gene;
#' * `is_downstream` indicates whether SNP is downstream of miRNA gene;
#' * `mirnaStrand` shows the strand;
#' * `mirnaStartPosition` displays the start position of DE-miRNA gene;
#' * `mirnaEndPosition` displays the end position of DE-miRNA gene.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' \donttest{
#' # search disease
#' searchDisease("Alzheimer disease")
#' disId <- "Alzheimer disease"
#'
#' # retrieve associated SNPs
#' association <- findMirnaSNPs(obj, disId)
#' }
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
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @importFrom rlang .data
#' @export
findMirnaSNPs <- function(mirnaObj,
                          diseaseEFO) {
    ## check inputs
    if (!is(mirnaObj, "MirnaExperiment")) {
        stop("'mirnaObj' should be of class MirnaExperiment! ",
             "See ?MirnaExperiment",
             call. = FALSE
        )
    }
    if (nrow(mirnaDE(mirnaObj, onlySignificant = FALSE)) == 0) {
        stop("MiRNA differential expression results are not present in ",
             "'mirnaObj'. Please, use 'performMirnaDE()' before using ",
             "this function. See ?performMirnaDE",
             call. = FALSE
        )
    }
    if (nrow(mirnaDE(mirnaObj)) == 0) {
        stop("No miRNAs in this object are defined as differentially ",
             "expressed. Consider using 'performMirnaDE()' with lower ",
             "thresholds. See ?performMirnaDE",
             call. = FALSE
        )
    }
    if (!is.character(diseaseEFO) | length(diseaseEFO) != 1) {
        stop("'diseaseEFO' must be a string of length 1 ",
             "containing the EFO ID of the disease. You can obtain the ",
             "disease identifier with 'searchDisease()'",
             call. = FALSE
        )
    }
    
    ## check that gwasrapidd is installed
    if (!requireNamespace("gwasrapidd", quietly = TRUE)) {
        stop(
            "The gwasrapidd package is needed to use ",
            "this function. Install it through: ",
            paste("`install.packages('gwasrapidd')`.",
                  sep = ""
            ),
            call. = FALSE
        )
    }
    
    ## check that biomaRt is installed
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
        stop("The biomaRt package is needed to use ",
             "this function. Install it through: ",
             paste("`BiocManager::install('biomaRt')`.",
                   sep = ""
             ),
             call. = FALSE
        )
    }
    
    ## inform the user about database querying
    message("Querying GWAS Catalog, this may take some time...")
    
    ## retrieve all SNPs related to the disease
    varDf <- gwasrapidd::get_variants(efo_trait = diseaseEFO)
    
    ## extract genomic context
    varDf <- varDf@genomic_contexts
    if (nrow(varDf) == 0) {
        stop("No variants exist for this biological trait!", call. = FALSE)
    }
    varDf <- as.data.frame(varDf)
    varDf <- varDf[varDf$chromosome_name %in% c(seq(22), "X", "Y"), seq(10)]
    varDf <- unique(varDf)
    
    ## only retain mapped genes
    varDf <- varDf[varDf$is_mapped_gene == "TRUE", ]
    
    ## find genomic details of differentially expressed miRNAs
    message("Finding genomic information of differentially expressed miRNAs...")
    ensembl <- biomaRt::useEnsembl(
        biomart = "ensembl",
        dataset = "hsapiens_gene_ensembl"
    )
    dem <- mirnaDE(mirnaObj)$ID
    deMirnas <- gsub("-5p|-3p", "", dem)
    mirGenes <- biomaRt::select(ensembl,
                                keys = deMirnas,
                                keytype = "mirbase_id",
                                columns = c(
                                    "hgnc_symbol",
                                    "mirbase_id",
                                    "gene_biotype",
                                    "chromosome_name",
                                    "strand",
                                    "start_position",
                                    "end_position"
                                )
    )
    mirGenes <- mirGenes[mirGenes$chromosome_name %in% c(seq(22), "X", "Y"), ]
    mirGenes$chromosome_name <- NULL
    
    ## intersect differentially expressed miRNAs with genes mapped with SNPs
    hostNames <- paste(mirGenes$hgnc_symbol, "HG", sep = "")
    mirMatches <- c(mirGenes$hgnc_symbol, hostNames)
    resDf <- varDf[which(varDf$gene_name %in% mirMatches), ]
    if (nrow(resDf) == 0) {
        message("No disease-related SNPs are present within DE-miRNA loci.")
        return(NULL)
    }
    idx2 <- lapply(mirGenes$hgnc_symbol, grep, resDf$gene_name)
    idx1 <- lapply(seq_along(idx2), function(i) rep(i, length(idx2[[i]])))
    resDf <- cbind(
        mirGenes[unlist(idx1), , drop = FALSE],
        resDf[unlist(idx2), , drop = FALSE]
    )
    
    ## add allele information
    snpData <- biomaRt::useEnsembl(biomart = "snps", dataset = "hsapiens_snp")
    snpInfo <- AnnotationDbi::select(snpData,
                                     keys = resDf$variant,
                                     keytype = "snp_filter",
                                     columns = c(
                                         "refsnp_id",
                                         "allele"
                                     )
    )
    colnames(snpInfo)[1] <- "variant_id"
    resDf <- merge(resDf, snpInfo, by = "variant_id")
    
    ## prepare resulting data.frame
    resDf <- resDf[, c(
        "variant_id", "gene_name", "hgnc_symbol", "mirbase_id",
        "chromosome_name", "chromosome_position", "allele",
        "distance", "is_upstream", "is_downstream", "strand",
        "start_position", "end_position"
    )]
    colnames(resDf) <- c(
        "variant", "gene", "miRNA.gene", "miRNA.precursor",
        "chr", "position", "allele", "distance", "is_upstream",
        "is_downstream", "mirnaStrand", "mirnaStartPosition",
        "mirnaEndPosition"
    )
    
    ## display a message with miRNA-SNPs association results
    message(
        "After the analysis, ", nrow(resDf), " variants associated with ",
        diseaseEFO, " were found within differentially ",
        "expressed miRNA genes"
    )
    
    ## return SNPs association data.frame
    return(resDf)
}





#' Get the scientific evidence for a particular disease-SNP association
#'
#' This function returns the biomedical evidence that supports the association
#' between a particular SNP and a phenotypic trait.
#'
#' @param variant The SNP ID of a particular variant of interest
#' (e.g. 'rs394581')
#' @param diseaseEFO The EFO identifier of a disease of interest. This can be
#' identified with the function [searchDisease()]
#'
#' @returns
#' A `tbl_df` dataframe containing information about literature evidences for a
#' disease-SNP association.
#'
#' @examples
#' \donttest{
#' # searchDisease("Alzheimer disease")
#' evidence <- getEvidence("rs2075650", diseaseEFO = "Alzheimer disease")
#' }
#'
#' @note
#' To retrieve evidences of disease-SNP association, this function makes use of
#' the `gwasrapidd` package.
#'
#' @references
#' Ramiro Magno, Ana-Teresa Maia, gwasrapidd: an R package to query, download
#' and wrangle GWAS catalog data, Bioinformatics, Volume 36, Issue 2, January
#' 2020, Pages 649–650, \url{https://doi.org/10.1093/bioinformatics/btz605}.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
getEvidence <- function(variant,
                        diseaseEFO) {
    ## check inputs
    if (!is.character(variant) | length(variant) != 1) {
        stop("'variant' should be a character with the ID of a SNP.",
             call. = FALSE)
    }
    if (!is.character(diseaseEFO) | length(diseaseEFO) != 1) {
        stop("'diseaseEFO' must be a string of length 1 ",
             "containing the EFO trait of the disease. You can obtain the ",
             "disease identifier with 'searchDisease()'",
             call. = FALSE
        )
    }
    
    ## check that gwasrapidd is installed
    if (!requireNamespace("gwasrapidd", quietly = TRUE)) {
        stop("The gwasrapidd package is needed to use ",
             "this function. Install it through: ",
             paste("`install.packages('gwasrapidd')`.",
                   sep = ""
             ),
             call. = FALSE
        )
    }
    
    ## retrieve biomedical evidence for a disease-SNP association
    message(
        "Retrieving biomedical evidence for the association between ",
        diseaseEFO, " and ", variant, " variant..."
    )
    ev <- gwasrapidd::get_studies(
        variant_id = variant,
        efo_trait = diseaseEFO
    )
    ev <- ev@publications
    
    ## returning results
    message(
        length(unique(ev$title)), " studies reporting this association ",
        "were found!"
    )
    return(ev)
}
