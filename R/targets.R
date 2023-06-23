#' Get microRNA targets
#'
#' This function allows to obtain human miRNA-target interactions using the
#' microRNA Data Integration Portal (mirDIP) database, which aggregates miRNA
#' target predictions from 24 different resources by using an integrated score
#' inferred from different prediction metrics. In this way, as demonstrated by
#' Tokar et al. 2018, mirDIP reports more accurate predictions compared to
#' those of individual tools.
#' 
#' @details
#' The downside of miRNA target prediction algorithms is the scarce extend of
#' overlap existing between the different tools. To address this issue,
#' several ensemble methods have been developed, trying to aggregate the
#' predictions obtained by different algorithms. Initially, several researchers
#' determined as significant miRNA-target pairs those predicted by more than
#' one tool (intersection method). However, this method is not able to capture
#' an important number of meaningful interactions. Alternatively, other
#' strategies used to merge predictions from several algorithms (union method).
#' Despite identifying more true relationships, the union method leads to a
#' higher proportion of false discoveries. Therefore, other ensemble methods
#' including mirDIP started using other statistics to rank miRNA-target
#' predictions obtained by multiple algorithms. For additional information on
#' mirDIP database and its ranking metric check Tokar et al. 2018 and
#' Hauschild et al. 2023.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param score The minimum mirDIP confidence score. It must be one of
#' `Very High`, `High`, `Medium` (default), `Low`, which correspond to ranks
#' among top 1%, top 5% (excluding top 1%), top 1/3 (excluding top 5%) and
#' remaining predictions, respectively.
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
#' # retrieve targets 
#' # obj <- getTargets(mirnaObj = obj)
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
#' @note
#' To access mirDIP database at \url{https://ophid.utoronto.ca/mirDIP/}, this
#' function directly use mirDIP API through R.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
getTargets <- function(mirnaObj,
                       score = "Medium") {
  
  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (nrow(geneDE(mirnaObj, onlySignificant = FALSE)) == 0) {
    stop(paste("Gene differential expression results are not present in",
               "'mirnaObj'. Please, use 'performGeneDE()' before using",
               "this function. See ?performGeneDE"), call. = FALSE)
  }
  if (!is.character(score) |
      length(score) != 1 |
      !score %in% c("Very High", "High", "Medium", "Low")) {
    stop(paste("'score' must be one of 'Very High', 'High', 'Medium', 'Low'.",
               "For additional details, see ?getTargets"),
         call. = FALSE)
  }
  
  ## define miRNAs
  allMirnas <- rownames(mirnaObj[["microRNA"]])
  
  ## collapse miRNA names
  microRNAs <- paste(allMirnas, collapse = ", ")
  
  ## set mirDIP database url
  url <- "http://ophid.utoronto.ca/mirDIP/Http_U"
  
  ## set mirDIP mapping score
  mapScore <- list("0", "1", "2", "3");
  names(mapScore) <- c("Very High", "High", "Medium", "Low")
  
  ## set API required parameters
  parameters <- list(
    genesymbol = "",
    microrna = microRNAs,
    scoreClass = mapScore[score]
  )
  
  ## send http POST throug 'getURL' and 'mirDIP.query' helper functions
  message("Retrieving targets from mirDIP database (this may take a while)...")
  res <- getURL(url, mirDIP.query, body = parameters, encode = "form")
  
  ## extract results from query
  response = httr::content(res, "text", encoding = "UTF-8")
  arr = unlist(strsplit(response, "\001", fixed = TRUE))
  
  ## convert results to a list object
  listMap <- lapply(arr, function(str) {
    arrKeyValue = unlist(strsplit(str, "\002", fixed = TRUE))
    if (length(arrKeyValue) > 1) {
      arrKeyValue[2]
    }
  })
  
  ## define the names of the retrieved values
  names(listMap) <- vapply(arr, function(str) {
    arrKeyValue = unlist(strsplit(str, "\002", fixed = TRUE))
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
  
  ## add miRNA-target pairs to the MirnaExperiment object
  mirnaTargets(mirnaObj) <- tg
  
  ## print the results of target retrieval
  message(paste(nrow(tg), "miRNA-target pairs have been identified for the",
                length(allMirnas), "miRNAs in study! Notably,",
                length(unique(mirnaTargets(mirnaObj)$Gene.Symbol)),
                "genes are targeted by differentially expressed miRNAs."))
  
  ## return mirnaObj with targets
  return(mirnaObj)
  
}





## helper function to query mirDIP and check for status
mirDIP.query <- function(URL, ...) {
  response <- httr::POST(URL, ...)
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
    if (!inherits(result, "error"))
      break
    N.TRIES <- N.TRIES - 1L
    message(paste("Attempting again to reach the resource..."))
  }
  
  ## if attempts are finished, print error message
  if (N.TRIES == 0L) {
    stop("'getURL()' failed:",
         "\n  URL: ", URL,
         "\n  error: ", conditionMessage(result))
  }
  
  ## return results
  return(result)
  
}

