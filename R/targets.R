#' Get microRNA targets
#'
#' This function allows to obtain the targets of differentially expressed
#' miRNAs using the [multiMiR] package, which outputs microRNA targets from
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
#' @param onlyValidated Logical, if `TRUE`, the function will retrieve only
#' experimentally validated targets; if `FALSE`, a combination of validated
#' and predicted targets will be used (default is `FALSE`)
#' @param minPredicted Considers as true targets only those predicted by at
#' least n prediction algorithms. It should be within 1-8 (default is `3`)
#' @param minValidated Considers as true targets only those present in more
#' than n databases with validated interactions. It should be within 1-3
#' (default is `1`)
#' @param topCutoff This value corresponds to the percentage of top results
#' retrieved by the different prediction databases present in [multiMiR].
#' Default is `20`
#'
#' @returns
#' A [`MirnaExperiment`][MirnaExperiment-class] object containing targets of
#' differentially expressed miRNAs. Targets are stored in the `targets` slot
#' and can be accessed with the [mirnaTargets()] function.
#'
#' @examples
#' # retrieve validated targets and targets predicted by at lest 3 databases
#' mirnaObj <- getTargets(mirnaObj = mirnaObj, onlyValidated = FALSE,
#' minPredicted = 3, minValidated = 1, topCutoff = 20)
#'
#' @references
#' Yuanbin Ru, Katerina J. Kechris, Boris Tabakoff, Paula Hoffman, Richard A.
#' Radcliffe, Russell Bowler, Spencer Mahaffey, Simona Rossi, George A. Calin,
#' Lynne Bemis, Dan Theodorescu; The multiMiR R package and database:
#' integration of microRNA–target interactions along with their disease and
#' drug associations. Nucleic Acids Res 2014; 42 (17): e133. doi:
#' [10.1093/nar/gku631]
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
getTargets <- function(mirnaObj,
                       specie = "Homo sapiens",
                       onlyValidated=FALSE,
                       minPredicted=3,
                       minValidated=1,
                       topCutoff=20) {

  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (!specie %in% convertOrganism("multiMiR", "all")) {
    stop(paste("'specie' must be one of:",
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

  ## selecting differentially expressed miRNAs from MirnaExperiment object
  mirnaTab <- mirnaDE(mirnaObj)

  ## obtaining targets with multiMiR according to specified parameters
  message(paste("Retrieving targets of differentially expressed miRNAs",
                "(this may take some time) ..."))
  if (onlyValidated == TRUE) {

    mirTar <- suppressWarnings(
      suppressMessages(
        quiet(
          multiMiR::get_multimir(org=specie,
                                 mirna=mirnaTab$ID,
                                 table="validated",
                                 predicted.cutoff=topCutoff,
                                 predicted.cutoff.type="p",
                                 summary=TRUE)
        )
      )
    )

    multiRes <- mirTar@summary[mirTar@summary$validated.sum >= minValidated, ]

  } else {

    mirTar <- suppressWarnings(
      suppressMessages(
        quiet(
          multiMiR::get_multimir(org=specie,
                                 mirna=mirnaTab$ID,
                                 table="all",
                                 predicted.cutoff=topCutoff,
                                 predicted.cutoff.type="p",
                                 summary=TRUE)
        )
      )
    )

    multiRes <- mirTar@summary[mirTar@summary$predicted.sum >= minPredicted |
                                 mirTar@summary$validated.sum >= minValidated, ]

  }

  ## append results to 'targets' slot and report number of targets found
  multiRes <- na.omit(multiRes[, c(2, 3)])
  multiRes <- unique(multiRes)
  message(paste(length(unique(multiRes$target_symbol)), "targets of the",
                length(mirnaTab$ID), "DE-miRNAs were found!"))
  mirnaTargets(mirnaObj) <- multiRes[order(multiRes$mature_mirna_id), ]
  return(mirnaObj)

}

