#' Get microRNA targets
#'
#' This function allows to obtain human miRNA-target interactions using two
#' databases, namely miRTarBase v9, which contains experimentally validated
#' interactions, and the microRNA Data Integration Portal (mirDIP) database,
#' which aggregates miRNA target predictions from 24 different resources by
#' using an integrated score inferred from different prediction metrics. In
#' this way, as demonstrated by Tokar et al. 2018, mirDIP reports more accurate
#' predictions compared to those of individual tools. However, for species
#' other than Homo sapiens only validated interactions are returned, since
#' mirDIP is only available for human miRNAs.
#'
#' @details
#' To define miRNA target genes, we can consider both experimentally validated
#' and computationally predicted interactions. Interactions of the former type
#' are generally preferred, since they are corroborated by biomolecular
#' experiments. However, they are often not sufficient, thus making it necessary
#' to consider the predicted interactions as well. The downside of miRNA target
#' prediction algorithms is the scarce extend of overlap existing between the
#' different tools. To address this issue, several ensemble methods have been
#' developed, trying to aggregate the predictions obtained by different
#' algorithms. Initially, several researchers determined as significant
#' miRNA-target pairs those predicted by more than one tool (intersection
#' method). However, this method is not able to capture an important number of
#' meaningful interactions. Alternatively, other strategies used to merge
#' predictions from several algorithms (union method). Despite identifying more
#' true relationships, the union method leads to a higher proportion of false
#' discoveries. Therefore, other ensemble methods including mirDIP started using
#' other statistics to rank miRNA-target predictions obtained by multiple
#' algorithms. For additional information on mirDIP database and its ranking
#' metric check Tokar et al. 2018 and Hauschild et al. 2023.
#'
#' This function defines miRNA targets by considering both validated
#' interactions present in miRTarBase (version 9), and predicted interactions
#' identified by mirDIP. Please note that for species other than Homo sapiens,
#' only miRTarBase interactions are available.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param organism The specie for which you are retrieving miRNA target genes.
#' Available species are: `Homo sapiens` (default), `Mus musculus`,
#' `Rattus norvegicus`, `Arabidopsis thaliana`, `Bos taurus`,
#' `Caenorhabditis elegans`, `Danio rerio`, `Drosophila melanogaster`,
#' `Gallus gallus`, `Sus scrofa`
#' @param score The minimum mirDIP confidence score. It must be one of
#' `Very High`, `High` (default), `Medium`, `Low`, which correspond to ranks
#' among top 1%, top 5% (excluding top 1%), top 1/3 (excluding top 5%) and
#' remaining predictions, respectively
#' @param includeValidated Logical, whether to include validated interactions
#' from miRTarBase or not. Default is TRUE in order to retrieve both predicted
#' and validated targets. Note that for species other than Homo sapines only
#' validated interactions are considered.
#'
#' @returns
#' A [`MirnaExperiment`][MirnaExperiment-class] object containing miRNA targets
#' stored in the `targets` slot. Results can be accessed with the
#' [mirnaTargets()] function.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' \donttest{
#' # retrieve targets
#' obj <- getTargets(mirnaObj = obj)
#' }
#'
#' # access targets
#' tg <- mirnaTargets(obj)
#'
#' @references
#' Tomas Tokar and others, mirDIP 4.1—integrative database of human microRNA
#' target predictions, Nucleic Acids Research, Volume 46, Issue D1, 4 January
#' 2018, Pages D360–D370, \url{https://doi.org/10.1093/nar/gkx1144}.
#'
#' Anne-Christin Hauschild and others, MirDIP 5.2: tissue context annotation
#' and novel microRNA curation, Nucleic Acids Research, Volume 51, Issue D1,
#' 6 January 2023, Pages D217–D225, \url{https://doi.org/10.1093/nar/gkac1070}.
#'
#' Hsi-Yuan Huang and others, miRTarBase update 2022: an informative resource
#' for experimentally validated miRNA–target interactions, Nucleic Acids
#' Research, Volume 50, Issue D1, 7 January 2022, Pages D222–D230,
#' \url{https://doi.org/10.1093/nar/gkab1079}.
#'
#' @note
#' To access mirDIP database at \url{https://ophid.utoronto.ca/mirDIP/}, this
#' function directly use mirDIP API through R.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
getTargets <- function(
        mirnaObj,
        organism = "Homo sapiens",
        score = "High",
        includeValidated = TRUE) {
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
    if (!is.character(organism) |
        length(organism) != 1 |
        !organism %in% c(
            "Homo sapiens", "Mus musculus", "Rattus norvegicus",
            "Arabidopsis thaliana", "Bos taurus",
            "Caenorhabditis elegans", "Danio rerio",
            "Drosophila melanogaster", "Gallus gallus",
            "Sus scrofa"
        )) {
        stop("'organism' must be  one of: 'Homo sapiens' (default), ",
            "'Mus musculus', 'Rattus norvegicus', 'Arabidopsis thaliana', ",
            "'Bos taurus', 'Caenorhabditis elegans', 'Danio rerio', ",
            "'Drosophila melanogaster', 'Gallus gallus', ",
            "and 'Sus scrofa'.",
            call. = FALSE
        )
    }
    if (!is.character(score) |
        length(score) != 1 |
        !score %in% c("Very High", "High", "Medium", "Low")) {
        stop("'score' must be one of 'Very High', 'High', 'Medium', 'Low'. ",
            "For additional details, see ?getTargets",
            call. = FALSE
        )
    }
    if (!is.logical(includeValidated) |
        length(includeValidated) != 1) {
        stop("'includeValidated' must be logical (TRUE/FALSE)!", call. = FALSE)
    }

    ## define miRNAs
    allMirnas <- mirnaDE(mirnaObj)$ID

    ## use only miRTarBase for organisms other than Homo sapiens
    if (organism != "Homo sapiens") {
        use.mirDIP <- FALSE
        includeValidated <- TRUE
        message(
            "For specie ", organism, " only miRTarBase database is ",
            "available..."
        )
    } else {
        use.mirDIP <- TRUE
    }

    ## use mirDIP for human target prediction
    if (use.mirDIP == TRUE) {
        ## collapse miRNA names
        microRNAs <- paste(allMirnas, collapse = ", ")

        ## set mirDIP database url
        url <- "http://ophid.utoronto.ca/mirDIP/Http_U"

        ## set mirDIP mapping score
        mapScore <- list("0", "1", "2", "3")
        names(mapScore) <- c("Very High", "High", "Medium", "Low")

        ## set API required parameters
        parameters <- list(
            genesymbol = "",
            microrna = microRNAs,
            scoreClass = mapScore[score]
        )

        ## send http POST throug 'getURL' and 'mirDIP.query' helper functions
        message("Retrieving targets from mirDIP (this may take a while)...")
        res <- getURL(url, mirDIP.query, body = parameters, encode = "form")

        ## extract results from query
        response <- httr::content(res, "text", encoding = "UTF-8")
        arr <- unlist(strsplit(response, "\001", fixed = TRUE))

        ## convert results to a list object
        listMap <- lapply(arr, function(str) {
            arrKeyValue <- unlist(strsplit(str, "\002", fixed = TRUE))
            if (length(arrKeyValue) > 1) {
                arrKeyValue[2]
            }
        })

        ## define the names of the retrieved values
        names(listMap) <- vapply(arr, function(str) {
            arrKeyValue <- unlist(strsplit(str, "\002", fixed = TRUE))
            if (length(arrKeyValue) > 1) {
                arrKeyValue[1]
            } else {
                item <- ""
                item
            }
        }, FUN.VALUE = character(1), USE.NAMES = FALSE)

        ## build a data.frame with miRNA-target pairs
        tg <- read.table(text = listMap$results, sep = "\t", header = TRUE)

        ## maintain only targets that are present in gene expression matrix
        tg <- tg[tg$Gene.Symbol %in% rownames(mirnaObj[["genes"]]), ]

        ## retain only interesting columns
        tg <- tg[, c(
            "Gene.Symbol", "MicroRNA", "Integrated.Score",
            "Number.of.Sources", "Score.Class"
        )]
        tg$Type <- "Predicted"
    }

    ## add validated interactions from miRTarBase
    if (includeValidated == TRUE) {
        ## define miRTarBase v9 link
        mtUrl <- paste("https://mirtarbase.cuhk.edu.cn/~miRTarBase/",
            "miRTarBase_2022/cache/download/9.0/miRTarBase_MTI.xlsx",
            sep = ""
        )

        ## load cache
        bfc <- .get_cache()

        ## check if miRTarBase is cached
        rid <- BiocFileCache::bfcquery(bfc, "miRTarBase", "rname")$rid
        if (!length(rid)) {
            ## download miRTarBase and add it to the cache directory
            message(
                "\nDownloading validated interactions ",
                "from miRTarBase v9.0..."
            )
            rid <- names(BiocFileCache::bfcadd(bfc, "miRTarBase", mtUrl))
        } else {
            message("\nLoading miRTarBase from cache...")
        }

        ## check if cached file needs to be updated
        if (!isFALSE(BiocFileCache::bfcneedsupdate(bfc, rid))) {
            BiocFileCache::bfcdownload(bfc, rid, ask = FALSE)
        }

        ## load miRTarBase
        mt <- quiet(readxl::read_xlsx(BiocFileCache::bfcrpath(bfc, rids = rid)))

        ## keep interactions involving measured miRNAs
        mt <- mt[mt$miRNA %in% allMirnas, c("miRNA", "Target Gene")]
        mt <- unique(mt)
        colnames(mt) <- c("MicroRNA", "Gene.Symbol")

        ## create resulting data.frame
        if (use.mirDIP == TRUE) {
            ## merge mirDIP and miRTarBase results
            message("Merging predicted and validated results...")
            tg <- merge(tg, mt, by = c("MicroRNA", "Gene.Symbol"), all = TRUE)
            tg$Type[tg$Type == "Predicted" & !is.na(tg$References)] <- "Both"
            tg$Type[is.na(tg$Type)] <- "Validated"
        } else {
            ## create a data.frame with only validated interactions
            tg <- mt
            tg$Type <- "Validated"
        }
    }

    ## add miRNA-target pairs to the MirnaExperiment object
    mirnaTargets(mirnaObj) <- tg

    ## print the results of target retrieval
    message(
        nrow(tg), " miRNA-target pairs have been identified for the ",
        length(allMirnas), " differentially expressed miRNAs."
    )

    ## return mirnaObj with targets
    return(mirnaObj)
}





## helper function to query mirDIP and check for status
mirDIP.query <- function(URL, ...) {
    response <- httr::POST(URL, config = httr::progress(), ...)
    httr::stop_for_status(response)
    response
}





## helper function to query web resources
getURL <- function(URL, FUN, ..., N.TRIES = 3L) {
    ## check that attempts are correctly defined
    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

    ## attempt to download the resource
    while (N.TRIES > 0L) {
        result <- tryCatch(FUN(URL, ...), error = identity)
        if (!inherits(result, "error")) {
            break
        }
        N.TRIES <- N.TRIES - 1L
        message("\nAttempting again to reach the resource...")
    }

    ## if attempts are finished, print error message
    if (N.TRIES == 0L) {
        stop(
            "\n'getURL()' failed:",
            "\n  URL: ", URL,
            "\n  code: ", conditionMessage(result)
        )
    }

    ## return results
    return(result)
}
