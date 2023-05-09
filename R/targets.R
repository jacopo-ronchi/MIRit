#' Get microRNA targets
#'
#' This function allows to obtain the targets of differentially expressed
#' miRNAs using the `multiMiR` package, which outputs microRNA targets from
#' several databases, including 8 prediction databases and 3 databases with
#' experimentally validated targets. The microRNA-mRNA interactions reported
#' by this function can derive by both predicted and validated interactions.
#' The user can decide whether to retrieve only validated targets or a
#' combination of predicted and validated targets. Moreover, for predicted
#' interactions, the user can decide either to consider as targets all
#' predicted targets or to include only those targets predicted by a minimum
#' of n prediction databases.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param organism The name of the organism under consideration. It must be one
#' of: `Homo sapiens`, `Mus musculus`, `Rattus norvegicus`
#' @param onlyValidated Logical, if `TRUE`, the function will retrieve only
#' experimentally validated targets; if `FALSE`, a combination of validated
#' and predicted targets will be used (default is `FALSE`)
#' @param minPredicted Considers as true targets only those predicted by at
#' least n prediction algorithms. It should be within 1-8 (default is `3`)
#' @param minValidated Considers as true targets only those present in more
#' than n databases with validated interactions. It should be within 1-3
#' (default is `1`)
#' @param topCutoff This value corresponds to the percentage of top results
#' retrieved by the different prediction databases present in `multiMiR`.
#' Default is `20`
#'
#' @returns
#' A [`MirnaExperiment`][MirnaExperiment-class] object containing targets of
#' differentially expressed miRNAs. Targets are stored in the `targets` slot
#' and can be accessed with the [mirnaTargets()] function.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # retrieve validated targets and targets predicted by at lest 3 databases
#' # obj <- getTargets(mirnaObj = obj, onlyValidated = FALSE,
#' # minPredicted = 3, minValidated = 1, topCutoff = 20)
#'
#' @references
#' Yuanbin Ru, Katerina J. Kechris, Boris Tabakoff, Paula Hoffman, Richard A.
#' Radcliffe, Russell Bowler, Spencer Mahaffey, Simona Rossi, George A. Calin,
#' Lynne Bemis, Dan Theodorescu; The multiMiR R package and database:
#' integration of microRNA–target interactions along with their disease and
#' drug associations. Nucleic Acids Res 2014; 42 (17): e133. doi:
#' \url{10.1093/nar/gku631}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
getTargets <- function(mirnaObj,
                       organism = "Homo sapiens",
                       onlyValidated=FALSE,
                       minPredicted=3,
                       minValidated=1,
                       topCutoff=20) {

  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (!organism %in% convertOrganism("multiMiR", "all")) {
    stop(paste("'organism' must be one of:",
               paste(convertOrganism("multiMiR", "all"), collapse = ", ")),
         call. = FALSE)
  }
  if (!is.logical(onlyValidated) |
      length(onlyValidated) != 1) {
    stop("'onlyValidated' must be logical. See ?getTargets", call. = FALSE)
  }
  if (!is.numeric(minPredicted) |
      length(minPredicted) != 1 |
      minPredicted > 8 |
      minPredicted < 0 |
      minPredicted != round(minPredicted)) {
    stop("'minPredicted' must be a whole number within 1-8. See ?getTargets",
         call. = FALSE)
  }
  if (!is.numeric(minValidated) |
      length(minValidated) != 1 |
      minValidated > 3 |
      minValidated < 0 |
      minValidated != round(minValidated)) {
    stop("'minValidated' must be a whole number within 1-3. See ?getTargets",
         call. = FALSE)
  }
  if (!is.numeric(topCutoff) |
      length(topCutoff) != 1 |
      topCutoff > 100 |
      topCutoff < 0) {
    stop("'topCutoff' must be a number between 1 and 100. See ?getTargets",
         call. = FALSE)
  }

  ## set the databases for miRNA-target pairs
  if (onlyValidated == TRUE) {
    dbs <- c("mirecords", "mirtarbase", "tarbase")
  } else {
    dbs <- c("mirecords", "mirtarbase", "tarbase", "diana_microt",
             "elmmo", "microcosm", "miranda", "mirdb", "pictar",
             "pita", "targetscan")
  }
  
  ## select differentially expressed miRNAs from MirnaExperiment object
  mirnaTab <- mirnaDE(mirnaObj)
  
  ## obtain targets with multiMiR according to specified parameters
  message(paste("Retrieving targets of differentially expressed miRNAs",
                "(this may take some time) ..."))
  
  multi <- lapply(dbs, function(x) {
    mirTar <- suppressWarnings(
      suppressMessages(
        quiet(
          multiMiR::get_multimir(org=organism,
                                 mirna=mirnaTab$ID,
                                 table=x,
                                 predicted.cutoff=topCutoff,
                                 predicted.cutoff.type="p",
                                 summary=TRUE)
        )
      )
    )
    mirTar@summary[, c("mature_mirna_id", "target_symbol")]
  })
  
  ## rename list elements
  names(multi) <- dbs
  
  ## remove duplicated hits
  multi <- lapply(multi, unique)
  
  ## summarize validated interactions
  validated <- rbind(multi[[1]], multi[[2]], multi[[3]])
  validated <- dplyr::group_by_all(validated)
  validated <- dplyr::count(validated)
  validated <- validated[validated$n >= minValidated,
                         c("mature_mirna_id", "target_symbol")]
  
  ## summarize predicted interactions and merge them with validated ones
  if (onlyValidated == FALSE) {
    predicted <- rbind(multi[[4]], multi[[5]], multi[[6]], multi[[7]],
                       multi[[8]], multi[[9]], multi[[10]], multi[[11]])
    predicted <- dplyr::group_by_all(predicted)
    predicted <- dplyr::count(predicted)
    predicted <- predicted[predicted$n >= minPredicted,
                           c("mature_mirna_id", "target_symbol")]
    targDf <- merge(validated, predicted, all = TRUE)
  } else {
    targDf <- validated
  }
  
  ## remove empty interactions
  targDf <- targDf[targDf$target_symbol != "", ]

  ## append results to 'targets' slot and report number of targets found
  message(paste(length(unique(targDf$target_symbol)), "targets of the",
                length(mirnaTab$ID), "DE-miRNAs were found!"))
  mirnaTargets(mirnaObj) <- targDf[order(targDf$mature_mirna_id), ]
  return(mirnaObj)

}

