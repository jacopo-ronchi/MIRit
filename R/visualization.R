#' Visualize the relationships between miRNAs and genes in a biological pathway
#'
#' This function is used to plot the network object created with the
#' [mirnaPathway()] function. The resulting network represents the interactions
#' occurring within a biological pathway alongside with regulating microRNAs
#' that affect pathway members. Moreover, this function visually displays the
#' network with gene and miRNA expression levels, so that the user can easily
#' understand the biological implications of miRNA dysregulations.
#'
#' @details
#' The `layout` parameter determines the layout of the integrated network. If
#' this parameter is set to `stress`, this function adopts an optimized version
#' of the stress-minimization algorithm (known as Kamada-Kawai algorithm).
#' Otherwise, if `circle` layout is used, biological pathways will be
#' displayed as a circular network. Both of these layout algorithms are
#' implemented in the `ggraph` package.
#'
#' The `stress` layout is suitable for small pathways, but it can get messed
#' for larger biological networks. In these cases, the `circle` layout
#' performs better.
#'
#' @param pathGraph A [`tidygraph::tbl_graph`] object created with the
#' [mirnaPathway()] function.
#' @param onlyIntegrated Logical. It must be set to `TRUE` for building the
#' network only for integrated targets, i.e. those targets whose expression is
#' statistically associated/correlated with miRNA expression. Instead, to build
#' the complete biological pathway, `onlyIntegrated` must be set to
#' `FALSE` (default)
#' @param layout It must be either `stress` (default), to build the network
#' through an optimized Kamada-Kawai algorithm, or `circle`, to represent the
#' pathway through a circular layout
#' @param node It must be either `box` (default), to represent nodes as
#' boxed labels, or `points`, to display nodes as points with the gene/miRNA
#' identifier
#' @param size The size of node elements in the graph. Default is `3`
#' @param title The title of the plot (e.g. "Thyroid hormone signaling
#' pathway"). Default is `NULL` not to include a plot title
#' @param edgesCol It must be an R color name that specifies the color of
#' interaction arrows in the network. Default is `gray`. All available colors
#' can be listed with [grDevices::colors()]
#' @param bindingEdgesCol It must be an R color name that specifies the color
#' of chemical binding arrows in the network. Default is `gray`. All available
#' colors can be listed with [grDevices::colors()]
#' @param stress.text.distance The distance between activation/inhibition
#' arrow caps and nodes for the `stress` design. Default is `0.2`
#' @param circle.text.distance The distance between activation/inhibition
#' arrow caps and nodes for the `circle` design. Default is `0.1`
#'
#' @returns
#' A `ggraph` object representing the network graph. For instructions on how
#' to handle and edit this object, please refer to the `ggraph` package.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # explore Thyroid hormone synthesis pathway present in KEGG
#' net <- mirnaPathway(obj,
#' pathway = "Thyroid hormone synthesis", organism = "Homo sapiens",
#' database = "KEGG")
#'
#' # visualize miRNA-gene network
#' visualizeNetwork(net)
#'
#' # visualize network with a circular layout and with points
#' visualizeNetwork(net, layout = "circle", node = "point")
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @importFrom rlang .data
#' @export
visualizeNetwork <- function(pathGraph,
                             onlyIntegrated = FALSE,
                             layout = "stress",
                             node = "box",
                             size = 3,
                             title = NULL,
                             edgesCol = "gray",
                             bindingEdgesCol = "gray",
                             stress.text.distance = 0.2,
                             circle.text.distance = 0.1) {

  ## check inputs
  if (!is(pathGraph, "tbl_graph")) {
    stop("'pathGraph' should be of class tbl_graph! See ?mirnaPathway",
         call. = FALSE)
  }
  if (!is.logical(onlyIntegrated) |
      length(onlyIntegrated) != 1) {
    stop("'onlyIntegrated' must be logical (TRUE/FALSE)! See ?visualizeNetwork",
         call. = FALSE)
  }
  if (!is.character(layout) |
      length(layout) != 1 |
      !layout %in% c("stress", "circle")) {
    stop(paste("'layout' must be either 'stress' or 'circle'.",
               "For additional details see ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!is.character(node) |
      length(node) != 1 |
      !node %in% c("box", "point")) {
    stop(paste("'node' must be either 'box' or 'point'.",
               "For additional details see ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!is.numeric(size) |
      length(size) != 1 |
      size < 0) {
    stop(paste("'size' must be a positive number!",
               "For additional details see ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot (e.g. 'Apoptosis').",
               "For additional details see ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!is.character(edgesCol) |
      length(edgesCol) != 1 |
      !edgesCol %in% grDevices::colors()) {
    stop(paste("'edgesCol' must be an R color name. All available colors",
               "can be listed with 'colors()'."),
         call. = FALSE)
  }
  if (!is.character(bindingEdgesCol) |
      length(bindingEdgesCol) != 1 |
      !bindingEdgesCol %in% grDevices::colors()) {
    stop(paste("'bindingEdgesCol' must be an R color name. All available",
               "colors can be listed with 'colors()'."),
         call. = FALSE)
  }
  if (!is.numeric(stress.text.distance) |
      length(stress.text.distance) != 1 |
      stress.text.distance < 0) {
    stop(paste("'stress.text.distance' must be a positive number!",
               "For additional details see ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!is.numeric(circle.text.distance) |
      length(circle.text.distance) != 1 |
      circle.text.distance < 0) {
    stop(paste("'circle.text.distance' must be a positive number!",
               "For additional details see ?visualizeNetwork"),
         call. = FALSE)
  }
  
  ## listing the variables not present in the global environment
  edges <- NULL

  ## set accepted pathway processes
  activation <- "ACTIVATION|Activation|activation|expression"
  inhibition <- "INHIBITION|Inhibition|inhibition"
  binding <- "Binding|binding"

  ## retain only integrated targets if wanted by the user
  if (onlyIntegrated == TRUE) {
    
    tmpGraph <- tidygraph::to_subgraph(pathGraph,
                                       .data$integrated == "integrated",
                                       subset_by = "nodes")$subgraph
    
    ## check if integrated axes are present, otherwise report the whole pathway
    if (length(tmpGraph) == 0) {
      warning(paste("There are no integrated axes within this pathway.",
                    "The complete pathway will be reported..."), call. = FALSE)
    } else {
      pathGraph <- tmpGraph
    }
  }

  ## extract network edges
  graphEdges <- as.data.frame(tidygraph::activate(pathGraph, edges))
  
  ## create a ggraph network to visualize biological pathways with DE-miRNAs
  mirNet <- ggraph::ggraph(pathGraph, layout = layout)
  
  ## use different edges for different biological interactions
  mirNet <- mirNet +
    ggraph::geom_edge_link(
      ggplot2::aes(filter = grepl(binding, .data$edgeType),
                   linetype = "binding"),
      color = bindingEdgesCol) +
    
    ggraph::geom_edge_link(
      ggplot2::aes(filter = grepl(inhibition, .data$edgeType) &
                     !grepl(activation, .data$edgeType),
                   linetype = "inhibition",
                   end_cap = {
                     ## set the arrow distance from nodes for boxed labels/points
                     if (node == "box") {
                       ggraph::label_rect(.data$node2.name,
                                          padding = ggplot2::margin(2, 2, 2, 2,
                                                                    unit = "mm"))
                     } else if (node == "point") {
                       ggraph::circle(3, unit = "mm")
                     }
                   }),
      color = edgesCol,
      arrow = grid::arrow(length = ggplot2::unit(2, "mm"),
                          angle = 90,
                          type = "open")) +
    
    ggraph::geom_edge_link(
      ggplot2::aes(filter = grepl(activation, .data$edgeType),
                   linetype = "activation",
                   end_cap = {
                     ## set the arrow distance from nodes for boxed labels/points
                     if (node == "box") {
                       ggraph::label_rect(.data$node2.name,
                                          padding = ggplot2::margin(2, 2, 2, 2,
                                                                    unit = "mm"))
                     } else if (node == "point") {
                       ggraph::circle(3, unit = "mm")
                     }
                   }),
      color = edgesCol,
      arrow = grid::arrow(length = ggplot2::unit(2, "mm"),
                          type = "closed")) +
    
    ggraph::geom_edge_link(
      ggplot2::aes(filter = !grepl(paste(activation,
                                         inhibition,
                                         binding,
                                         sep = "|"), .data$edgeType),
                   linetype = "standard"),
      color = edgesCol)

  ## customize colors and labels with different scales for genes and miRNAs
  if (node == "box") {
    mirNet <- mirNet +
      ggraph::geom_node_label(ggplot2::aes(label = .data$name,
                                           fill = .data$geneLogFC),
                              size = size) +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                    guide = ggplot2::guide_colorbar(order = 1)) +
      ggnewscale::new_scale_fill() +
      ggraph::geom_node_label(data = function(x) {tidygraph::filter(x, is.na(.data$geneLogFC))},
                              ggplot2::aes(label = .data$name,
                                           fill = .data$mirnaLogFC),
                              size = size) +
      ggplot2::scale_fill_gradient2(low = "green", mid = "grey", high = "yellow",
                                    guide = ggplot2::guide_colorbar(order = 2))
  } else if (node == "point") {
    mirNet <- mirNet +
      ggraph::geom_node_point(ggplot2::aes(fill = .data$geneLogFC),
                              pch = 21, size = size) +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                    guide = ggplot2::guide_colorbar(order = 1)) +
      ggnewscale::new_scale_fill() +
      ggraph::geom_node_point(data = function(x) {tidygraph::filter(x, is.na(.data$geneLogFC))},
                              ggplot2::aes(fill = .data$mirnaLogFC),
                              pch = 21, size = size) +
      ggplot2::scale_fill_gradient2(low = "green", mid = "grey", high = "yellow",
                                    guide = ggplot2::guide_colorbar(order = 2))
  }

  ## add text elements for point nodes
  if (node == "point") {
    if (layout == "stress") {
      mirNet <- mirNet +
        ggraph::geom_node_text(ggplot2::aes(label = .data$name),
                               size = size, nudge_y = - stress.text.distance)
    } else if (layout == "circle") {
      mirNet <- mirNet +
        ggraph::geom_node_text(ggplot2::aes(label = .data$name), size = size,
                               nudge_x = mirNet$data$x * circle.text.distance,
                               nudge_y = mirNet$data$y * circle.text.distance)
    }
  }

  ## add linetype scale
  mirNet <- mirNet +
    ggraph::scale_edge_linetype_manual(name = "interaction",
                                       values = c("activation" = "solid",
                                                  "binding" = "dashed",
                                                  "inhibition" = "solid",
                                                  "standard" = "solid"),
                                       guide = "none")

  ## increase plot limits on the basis of maximum and minimum coordinates
  mirNet <- mirNet +
    ggplot2::coord_cartesian(xlim = c(-0.5 + min(mirNet$data$x),
                                      0.5 + max(mirNet$data$x)),
                             ylim = c(-0.5 + min(mirNet$data$y),
                                      0.5 + max(mirNet$data$y)))

  ## add graph and title themes
  mirNet <- mirNet +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "plain",
                                                      hjust = 0.5))

  ## add the title of the plot
  if (!is.null(title)) {
    mirNet <- mirNet +
      ggplot2::ggtitle(title)
  }

  ## plot network
  return(mirNet)

}





#' @rdname mirnaDotplot
#' @export
setMethod("mirnaDotplot", "MirnaEnrichment",
          function(object,
                   showTerms,
                   splitDir,
                   ordBy,
                   sizeBy,
                   colBy,
                   title) {
            mirnaDotplot.MirnaEnrichment(object,
                                         showTerms = showTerms,
                                         splitDir = splitDir,
                                         ordBy = ordBy,
                                         sizeBy = sizeBy,
                                         colBy = colBy,
                                         title = title)
          })





#' @rdname mirnaDotplot
#' @export
setMethod("mirnaDotplot", "MirnaGsea",
          function(object,
                   showTerms,
                   splitDir,
                   ordBy,
                   sizeBy,
                   colBy,
                   title) {
            mirnaDotplot.MirnaGsea(object,
                                   showTerms = showTerms,
                                   splitDir = splitDir,
                                   ordBy = ordBy,
                                   sizeBy = sizeBy,
                                   colBy = colBy,
                                   title = title)
          })





## create a dotplot for MirnaEnrichment objects
#' @importFrom rlang .data
mirnaDotplot.MirnaEnrichment <- function(mirnaEnr,
                                         showTerms,
                                         splitDir,
                                         ordBy,
                                         sizeBy,
                                         colBy,
                                         title) {

  ## check inputs
  if (!is(mirnaEnr, "MirnaEnrichment")) {
    stop("'mirnaEnr' should be of class MirnaEnrichment! See ?enrichMirnas.",
         call. = FALSE)
  }
  if (nrow(enrichmentResults(mirnaEnr)) < 1) {
    stop("'mirnaEnr' object does not contain any significant results!",
         call. = FALSE)
  }
  if (!is.character(showTerms) &
      !(is.numeric(showTerms) & length(showTerms) == 1)) {
    stop(paste("'showTerms' must be the number of top enriched terms",
               "to plot or, alternatively, a character vector containing",
               "the terms to be shown."),
         call. = FALSE)
  }
  if (is.character(showTerms) &
      !all(showTerms %in% enrichmentResults(mirnaEnr)$Subcategory)) {
    stop("the terms provided are not present in 'mirnaEnr' Subcategory column!",
         call. = FALSE)
  }
  if (!is.logical(splitDir) |
      length(splitDir) != 1) {
    stop("'splitDir' must be logical (TRUE/FALSE)! See ?mirnaDotplot",
         call. = FALSE)
  }
  if (!is.character(ordBy) |
      length(ordBy) != 1 |
      !ordBy %in% c("fold", "P.adjusted", "P.value", "Observed")) {
    stop(paste("'ordBy' must be one of: 'fold' (default),",
               "'P.adjusted', 'P.value', 'Observed'"),
         call. = FALSE)
  }
  if (!is.character(sizeBy) |
      length(sizeBy) != 1 |
      !sizeBy %in% c("fold", "P.adjusted", "P.value", "Observed")) {
    stop(paste("'sizeBy' must be one of: 'fold', 'P.adjusted',",
               "'P.value', 'Observed' (default)"),
         call. = FALSE)
  }
  if (!is.character(colBy) |
      length(colBy) != 1 |
      !colBy %in% c("fold", "P.adjusted", "P.value", "Observed")) {
    stop(paste("'colBy' must be one of: 'fold', 'P.adjusted'",
               "(default), 'P.value', 'Observed'"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot (e.g. 'Enrichment').",
               "For additional details see ?mirnaDotplot"),
         call. = FALSE)
  }

  ## set the right parameter value for 'fold'
  if (ordBy == "fold") ordBy <- "foldEnrichment"
  if (sizeBy == "fold") sizeBy <- "foldEnrichment"
  if (colBy == "fold") colBy <- "foldEnrichment"

  ## extract results from enrichment object
  res <- enrichmentResults(mirnaEnr)

  ## warning if there are more than one category
  if (length(unique(res$Category)) > 1) {
    warning(paste("The object contains results from more than one category.",
                  "If you want to plot enriched terms just from one miEAA",
                  "category you can use: mirnaEnr['category']"), call. = FALSE)
  }

  ## compute ratio between observed and expected hits
  res$foldEnrichment <- res$Observed / res$Expected

  ## set an x-axis label for 'foldEnrichment'
  if (ordBy == "foldEnrichment") {
    ordLabel <- "Fold Enrichment"
  } else {
    ordLabel <- ordBy
  }

  ## reorder results based on specified criterion
  if (ordBy != "P.adjusted") {
    res <- res[order(res[, ordBy], decreasing = TRUE), ]
  } else {
    res <- res[order(res[, ordBy], decreasing = FALSE), ]
  }

  ## select terms to be shown in the dotplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$Subcategory %in% showTerms), ]
  }

  ## create a dotplot
  dotRes <- ggplot2::ggplot(res,
                            ggplot2::aes(x = !!ggplot2::sym(ordBy),
                                         y = stats::reorder(.data$Subcategory,
                                                            !!ggplot2::sym(ordBy)),
                                         size = !!ggplot2::sym(sizeBy),
                                         color = !!ggplot2::sym(colBy))) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "red", high = "blue",
                                  guide = ggplot2::guide_colorbar(reverse = TRUE)) +
    ggplot2::ylab(NULL) +
    ggplot2::scale_y_discrete() +
    ggplot2::scale_size(range = c(3, 8)) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 1),
                    color = ggplot2::guide_colorbar(order = 2)) +
    ggplot2::xlab(ordLabel) +
    theme_enr()

  ## divide by enrichment direction
  if (splitDir == TRUE) {
    dotRes <- dotRes +
      ggplot2::facet_grid(~ Enrichment) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 12))
  }
  
  ## add the title of the plot
  if (!is.null(title)) {
    dotRes <- dotRes +
      ggplot2::ggtitle(title)
  }

  ## return ggplot2 graph
  return(dotRes)

}





## create a dotplot for MirnaGsea objects
#' @importFrom rlang .data
mirnaDotplot.MirnaGsea <- function(mirnaGsea,
                                   showTerms,
                                   splitDir,
                                   ordBy,
                                   colBy,
                                   sizeBy,
                                   title) {

  ## check inputs
  if (!is(mirnaGsea, "MirnaGsea")) {
    stop("'mirnaGsea' should be of class MirnaGsea! See ?enrichMirnas.",
         call. = FALSE)
  }
  if (nrow(enrichmentResults(mirnaGsea)) < 1) {
    stop("'mirnaGsea' object does not contain any significant results!",
         call. = FALSE)
  }
  if (!is.character(showTerms) &
      !(is.numeric(showTerms) & length(showTerms) == 1)) {
    stop(paste("'showTerms' must be the number of top enriched terms",
               "to plot or, alternatively, a character vector containing",
               "the terms to be shown."),
         call. = FALSE)
  }
  if (is.character(showTerms) &
      !all(showTerms %in% enrichmentResults(mirnaGsea)$Subcategory)) {
    stop("the terms provided are not present in 'mirnaGsea' Subcategory column!",
         call. = FALSE)
  }
  if (!is.logical(splitDir) |
      length(splitDir) != 1) {
    stop("'splitDir' must be logical (TRUE/FALSE)! See ?mirnaDotplot",
         call. = FALSE)
  }
  if (!is.character(ordBy) |
      length(ordBy) != 1 |
      !ordBy %in% c("fold", "P.adjusted", "P.value", "Observed")) {
    stop(paste("'ordBy' must be one of: 'fold' (default),",
               "'P.adjusted', 'P.value', 'Observed'"),
         call. = FALSE)
  }
  if (!is.character(sizeBy) |
      length(sizeBy) != 1 |
      !sizeBy %in% c("fold", "P.adjusted", "P.value", "Observed")) {
    stop(paste("'sizeBy' must be one of: 'fold', 'P.adjusted', 'P.value',",
               "'Observed' (default)"),
         call. = FALSE)
  }
  if (!is.character(colBy) |
      length(colBy) != 1 |
      !colBy %in% c("fold", "P.adjusted", "P.value", "Observed")) {
    stop(paste("'colBy' must be one of: 'fold', 'P.adjusted' (default),",
               "'P.value', 'Observed'"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot (e.g. 'Enrichment').",
               "For additional details see ?mirnaDotplot"),
         call. = FALSE)
  }

  ## set the right parameter value for 'fold'
  if (ordBy == "fold") ordBy <- "meanFC"
  if (sizeBy == "fold") sizeBy <- "meanFC"
  if (colBy == "fold") colBy <- "meanFC"

  ## extract results from enrichment object
  res <- enrichmentResults(mirnaGsea)

  ## warning if there are more than one category
  if (length(unique(res$Category)) > 1) {
    warning(paste("The object contains results from more than one category.",
                  "If you want to plot enriched terms just from one miEAA",
                  "category you can use: mirnaGsea['category']"),
            call. = FALSE)
  }

  ## extract mirna fold changes
  mirnaID <- mirnaIdEnrichment(mirnaGsea)
  lfc <- lfcEnrichment(mirnaGsea)

  ## compute mean absolute fold change in gene sets
  res$meanFC <- vapply(rownames(res), function(x) {
    mirSet <- res[x, "miRNAs/precursors"]
    mirSet <- stringr::str_split(mirSet, "; ", simplify = TRUE)
    mean(lfc[mirnaID %in% mirSet])
  }, FUN.VALUE = numeric(1))
  res$meanFC <- abs(res$meanFC)

  ## set an x-axis label for 'meanFC'
  if (ordBy == "meanFC") {
    ordLabel <- "Mean Absolute logFC"
  } else {
    ordLabel <- ordBy
  }

  ## reorder results based on specified criterion
  if (ordBy != "P.adjusted") {
    res <- res[order(res[, ordBy], decreasing = TRUE), ]
  } else {
    res <- res[order(res[, ordBy], decreasing = FALSE), ]
  }

  ## select terms to be shown in the dotplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$Subcategory %in% showTerms), ]
  }

  ## create a dotplot
  dotRes <- ggplot2::ggplot(res,
                            ggplot2::aes(x = !!ggplot2::sym(ordBy),
                                         y = stats::reorder(.data$Subcategory,
                                                            !!ggplot2::sym(ordBy)),
                                         size = !!ggplot2::sym(sizeBy),
                                         color = !!ggplot2::sym(colBy))) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "red", high = "blue",
                                  guide = ggplot2::guide_colorbar(reverse = TRUE)) +
    ggplot2::ylab(NULL) +
    ggplot2::scale_y_discrete() +
    ggplot2::scale_size(range = c(3, 8)) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 1),
                    color = ggplot2::guide_colorbar(order = 2)) +
    ggplot2::xlab(ordLabel) +
    theme_enr()

  ## divide by enrichment direction
  if (splitDir == TRUE) {
    dotRes <- dotRes +
      ggplot2::facet_grid(~ Enrichment) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 12))
  }
  
  ## add the title of the plot
  if (!is.null(title)) {
    dotRes <- dotRes +
      ggplot2::ggtitle(title)
  }

  ## return ggplot2 graph
  return(dotRes)

}





#' Create a ridgeplot to display the results of miRNA GSEA analysis
#'
#' This function creates a ridgeplot that is useful for showing the results
#' of a miRNA GSEA analysis. The output of this function is a plot where
#' enrichmed terms found with [gseaMirnas()] function are visualized on the
#' basis of miRNA fold changes. In particular, the plot displays enriched
#' terms on the y-axis and miRNA fold changes on the x-axis. The resulting
#' area represents the density of fold changes belonging to miRNAs annotated
#' to that category.
#'
#' @param mirnaGsea A [`MirnaGsea`][MirnaGsea-class] object containing the
#' results of gene set enrichment analysis obtained with the [gseaMirnas()]
#' function
#' @param showTerms It is the number of top enriched terms to show or,
#' alternatively, a character vector indicating the enriched terms to plot.
#' Default is `10`
#' @param colBy The parameter used to set the color scale. It must be one of
#' `P.adjusted` (default), `P.value` and `Observed`
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' An object of class `ggplot` containing the ridgeplot of miRNA GSEA results.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform miRNA GSEA analysis
#' gse_res <- gseaMirnas(obj, organism = "Homo sapiens",
#' category = "GO_Annotations_mature")
#'
#' # plot results as a ridgeplot
#' mirnaRidgeplot(gse_res)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @importFrom rlang .data
#' @export
mirnaRidgeplot <- function(mirnaGsea,
                           showTerms = 10,
                           colBy = "P.adjusted",
                           title = NULL) {

  ## check inputs
  if (!is(mirnaGsea, "MirnaGsea")) {
    stop("'mirnaGsea' should be of class MirnaGsea! See ?gseaMirnas.",
         call. = FALSE)
  }
  if (nrow(enrichmentResults(mirnaGsea)) < 1) {
    stop("'mirnaGsea' object does not contain any significant results!",
         call. = FALSE)
  }
  if (!is.character(showTerms) &
      !(is.numeric(showTerms) & length(showTerms) == 1)) {
    stop(paste("'showTerms' must be the number of top enriched terms",
               "to plot or, alternatively, a character vector containing",
               "the terms to be shown."),
         call. = FALSE)
  }
  if (is.character(showTerms) &
      !all(showTerms %in% enrichmentResults(mirnaGsea)$Subcategory)) {
    stop("the terms provided are not present in 'mirnaGsea' Subcategory column!",
         call. = FALSE)
  }
  if (!is.character(colBy) |
      length(colBy) != 1 |
      !colBy %in% c("P.adjusted", "P.value", "Observed")) {
    stop(paste("'colBy' must be one of: P.adjusted' (default),",
               "'P.value', 'Observed'"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot (e.g. 'Enrichment').",
               "For additional details see ?mirnaRidgeplot"),
         call. = FALSE)
  }

  ## extract results from enrichment object
  res <- enrichmentResults(mirnaGsea)

  ## warning if there are more than one category
  if (length(unique(res$Category)) > 1) {
    warning(paste("The object contains results from more than one category.",
                  "If you want to plot enriched terms just from one miEAA",
                  "category you can use: mirnaEnr['category']"), call. = FALSE)
  }

  ## select terms to be shown in the dotplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$Subcategory %in% showTerms), ]
  }

  ## extract mirna fold changes
  mirnaID <- mirnaIdEnrichment(mirnaGsea)
  lfc <- lfcEnrichment(mirnaGsea)

  ## retrieve logFCs for miRNAs in each category
  lfcSet <- lapply(res[, "miRNAs/precursors"], function(x) {
    mirSet <- stringr::str_split(x, "; ", simplify = TRUE)
    lfc[mirnaID %in% mirSet]
  })

  ## create a dataframe with logFCs reported for each category
  setLen <- vapply(lfcSet, length, FUN.VALUE = numeric(1))
  meanLfc <- vapply(lfcSet, mean, FUN.VALUE = numeric(1))
  ridgeDf <- data.frame(term = rep(res$Subcategory, setLen),
                        P.adjusted = rep(res$P.adjusted, setLen),
                        P.value = rep(res$P.value, setLen),
                        Observed = rep(res$Observed, setLen),
                        meanLFC = rep(meanLfc, setLen),
                        val = unlist(lfcSet))

  ## create a ridgeplot
  ridgeRes <- ggplot2::ggplot(ridgeDf, ggplot2::aes(x = .data$val,
                                                    y = stats::reorder(.data$term,
                                                                       .data$meanLFC),
                                                    fill = !!ggplot2::sym(colBy))) +
    ggridges::geom_density_ridges() +
    ggplot2::scale_fill_continuous(low = "red",
                                   high = "blue",
                                   name = "P.adjusted",
                                   guide = ggplot2::guide_colorbar(reverse = TRUE)) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("miRNAs logFC") +
    theme_enr()
  
  ## add the title of the plot
  if (!is.null(title)) {
    ridgeRes <- ridgeRes +
      ggplot2::ggtitle(title)
  }

  ## return ggplot2 graph
  return(ridgeRes)

}





## ggplot2 theme for graphs
theme_enr <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(colour = "black",
                                          size = 12,
                                          vjust = 1),
      axis.text.y = ggplot2::element_text(colour = "black",
                                          size = 12,
                                          hjust = 1),
      axis.title = ggplot2::element_text(colour = "black",
                                         size = 12,
                                         margin = ggplot2::margin(10, 5, 0, 0)),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}





#' Create a trackplot tho show association between miRNAs and disease-SNPs
#'
#' This function plots a trackplot that shows the genomic position of
#' disease-associated SNPs that affect miRNA genes. This is useful to visualize
#' the genomic position and context of disease-associated variants that may
#' affect miRNA expression.
#'
#' @param variantId A valid name of a SNP variant! (e.g. `"rs394581"`)
#' @param snpAssociation A `data.frame` object containing the results of
#' [findMirnaSNPs()] function
#' @param showContext Logical, if `TRUE` a complete genomic context with genes
#' present in the region will be shown. Default is `FALSE` to just display the
#' variant and the miRNA gene
#' @param showSequence Logical, whether to display a color-coded sequence at
#' the bottom of the trackplot. Default is `TRUE`. This parameter will be set
#' to `FALSE` if `showContext` is `TRUE`
#' @param snpFill It must be an R color name that specifies the fill color of
#' the SNP locus. Default is `lightblue`. All available colors can be listed
#' with [grDevices::colors()]
#' @param mirFill It must be an R color name that specifies the fill color of
#' the miRNA locus. Default is `orange`. All available colors can be listed
#' with [grDevices::colors()]
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#' @param ... Other parameters that can be passed to [Gviz::plotTracks()]
#' function
#'
#' @note
#' This function retrieves genomic coordinates from the output of
#' [findMirnaSNPs()] function and then uses `Gviz` package to build
#' the trackplot.
#'
#' @references
#' Hahne, F., Ivanek, R. (2016). Visualizing Genomic Data Using Gviz and
#' Bioconductor. In: Mathé, E., Davis, S. (eds) Statistical Genomics. Methods
#' in Molecular Biology, vol 1418. Humana Press, New York, NY.
#' \url{https://doi.org/10.1007/978-1-4939-3578-9_16}
#'
#' @returns
#' A trackplot with information about chromosome, SNP and miRNA gene location.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # retrieve associated SNPs
#' # association <- findMirnaSNPs(obj, disId)
#'
#' # visualize association as a trackplot
#' # mirVariantPlot(variantId = varId, snpAssociation = association)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
mirVariantPlot <- function(variantId,
                           snpAssociation,
                           showContext = FALSE,
                           showSequence = TRUE,
                           snpFill = "lightblue",
                           mirFill = "orange",
                           title = NULL,
                           ...) {

  ## check inputs
  if (!is.character(variantId) |
      length(variantId) != 1 |
      !startsWith(variantId, "rs")) {
    stop("'variantId' must be a valid name of a SNP variant! (e.g. 'rs394581')",
         call. = FALSE)
  }
  if (!is.data.frame(snpAssociation) |
      !identical(colnames(snpAssociation), c("variant", "disease", "mirnaGene",
                                             "allele", "chromosome", "position",
                                             "genesPresent", "mirnaStart",
                                             "mirnaEnd", "mirnaStrand",
                                             "variantDSI", "variantDPI",
                                             "ei", "score"))) {
    stop(paste("'snpAssociation' must be a data.frame containing the",
               "list of SNPs occurring at DE-miRNA genes. To obtain this",
               "association you can use the 'findMirnaSNPs()' function.",
               "See ?findMirnaSNPs for details."),
         call. = FALSE)
  }
  if (nrow(snpAssociation) < 1) {
    stop("'snpAssociation' does not contain any significant results!",
         call. = FALSE)
  }
  if (!is.logical(showContext) |
      length(showContext) != 1) {
    stop("'showContext' must be logical (TRUE/FALSE)! See ?mirVariantPlot",
         call. = FALSE)
  }
  if (!is.logical(showSequence) |
      length(showSequence) != 1) {
    stop("'showSequence' must be logical (TRUE/FALSE)! See ?mirVariantPlot",
         call. = FALSE)
  }
  if (!is.character(snpFill) |
      length(snpFill) != 1 |
      !snpFill %in% grDevices::colors()) {
    stop(paste("'snpFill' must be an R color name. All available colors",
               "can be listed with 'colors()'. Default is: 'lightblue'"),
         call. = FALSE)
  }
  if (!is.character(mirFill) |
      length(mirFill) != 1 |
      !mirFill %in% grDevices::colors()) {
    stop(paste("'mirFill' must be an R color name. All available colors",
               "can be listed with 'colors()'. Default is: 'orange'"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot (e.g. 'SNPs overlap').",
               "For additional details see ?mirVariantPlot"),
         call. = FALSE)
  }

  ## select the specified variant
  snpAssociation <- snpAssociation[snpAssociation$variant == variantId, ]

  ## extract and format SNP genomic locations
  snpSeq <- snpAssociation[, c(5, 6, 6, 1)]
  colnames(snpSeq) <- c("chr", "start", "end", "variant")
  snpSeq$strand <- "*"

  ## create GRanges object containing SNP positions
  snpSeq <- GenomicRanges::makeGRangesFromDataFrame(snpSeq,
                                                    keep.extra.columns = TRUE)

  ## extract and format miRNA genomic coordinates
  mirSeq <- snpAssociation[, c(5, 8, 9, 10, 3)]
  colnames(mirSeq) <- c("chr", "start", "end", "strand", "miRNA_gene")

  ## create GRanges object containing miRNA positions
  mirSeq <- GenomicRanges::makeGRangesFromDataFrame(mirSeq,
                                                    keep.extra.columns = TRUE)

  ## set parameters
  chr <- paste("chr", snpAssociation$chromosome, sep = "")
  g <- "hg38"
  lf <- 0.1
  rf <- 0.1

  ## create ideomTrack
  iTrack <- Gviz::IdeogramTrack(genome = g, chromosome = chr, lwd = 4)

  ## create genome axis track
  gTrack <- Gviz::GenomeAxisTrack()

  ## create SNPs track
  snpTrack <- Gviz::GeneRegionTrack(snpSeq,
                                    name = "SNP",
                                    symbol = paste(variantId, "   "),
                                    fill = snpFill)

  ## create miRNA track
  mirTrack <- Gviz::GeneRegionTrack(mirSeq,
                                    genome = g,
                                    chromosome = chr,
                                    name = "miRNA",
                                    symbol = snpAssociation$mirnaGene)

  ## create sequence track
  if (showSequence == TRUE) {
    sTrack <- Gviz::SequenceTrack(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                                  chromosome = chr)
  }

  ## create genomic context track
  if (showContext == TRUE) {
    biomTrack <- Gviz::BiomartGeneRegionTrack(genome = g,
                                              symbol = snpAssociation$mirnaGene,
                                              name = "Gene context")
  }

  ## create trackplot element list
  if (showContext == TRUE) {
    pList <- list(iTrack, gTrack, biomTrack, snpTrack)
    lf <- 0
    rf <- 0
  } else {
    if (showSequence == TRUE) {
      pList <- list(iTrack, gTrack, snpTrack, mirTrack, sTrack)
    } else if (showSequence == FALSE) {
      pList <- list(iTrack, gTrack, snpTrack, mirTrack)
    }
  }

  ## create the trackplot object
  trackPlot <- Gviz::plotTracks(pList,
                                extend.left = lf,
                                extend.right = rf,
                                transcriptAnnotation = "symbol",
                                main = title,
                                ...)


}





#' Plot correlation between miRNAs and genes within biological groups
#' 
#' This function creates a scatter plot that shows the correlation between
#' miRNA and gene expression levels. This is useful after correlation
#' analysis performed through the [integrateMirnaTargets()] function, to
#' graphically visualize the quantitative effect of miRNA dysregulations on
#' target gene expression. Furthermore, this function performs linear/monotonic
#' regression to better represent the relationships between miRNA-target pairs.
#' 
#' When non-parametric correlation has been performed with the
#' [integrateMirnaTargets()] function, a regression line can be fitted through
#' monotonic regression on expression levels, or through linear regression
#' performed on rank-transformed data. Since, ranks do not correspond to real
#' expression values, the default option is to perform monotonic regression
#' to fit a monotonic curve. To do so, this function makes use of the `MonoPoly`
#' R package, which implements the algorithm proposed by Murray et al. in 2016.
#' 
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param mirna The name of the miRNA for which we want to observe the
#' correlation
#' @param gene The name of the gene for which we want to observe the
#' correlation
#' @param condition It must be the column name of a variable specified in the
#' metadata (colData) of a [`MirnaExperiment`][MirnaExperiment-class] object;
#' or, alternatively, it must be a character/factor object that specifies
#' group memberships (eg. c("healthy, "healthy", "disease", "disease"))
#' @param showCoeff Logical, whether to show the correlation coeffficient or
#' not. Note that the "R" is used for Pearson's correlation", "rho" for
#' Spearman's correlation, and "tau" for Kendall's correlation. Default is TRUE
#' @param regression Logical, whether to display a linear/monotonic regression
#' line that fits miRNA-gene correlation data. Default is TRUE
#' @param useRanks Logical, whether to represent non-parametric correlation
#' analyses (Spearman's and Kendall's correlations) through rank-transformed
#' data. Note that in this case, linear regression is performed on ranked
#' data instead of monotonic regression. Default is FALSE
#' @param lineCol It must be an R color name that specifies the color of
#' the regression line. Default is `red`. All available colors can be listed
#' with [grDevices::colors()]
#' @param lineType It specifies the line type used for th regression line. It
#' must be either 'blank', 'solid', 'dashed' (default), 'dotted', 'dotdash',
#' 'longdash' or 'twodash'
#' @param lineWidth The width of the fitted regression line (default is 0.8)
#' @param pointSize The size of points in the correlation plot (default is 3)
#' @param colorScale It must be a named character vector where values
#' correspond to R colors, while names coincide with the groups specified in
#' the `condition` parameter (eg. c("healthy" = "green", "disease" = "red")).
#' Default is NULL, in order to use the default color scale
#' 
#' @references
#' K. Murray, S. Müller & B. A. Turlach (2016) Fast and flexible methods for
#' monotone polynomial fitting, Journal of Statistical Computation and
#' Simulation, 86:15, 2946-2966, DOI: \url{10.1080/00949655.2016.1139582}.
#'
#' @returns
#' An object of class `ggplot` containing the correlation scatter plot.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform miRNA-target integration
#' obj <- integrateMirnaTargets(obj)
#'
#' # plot correlation between miR-34a and PAX8 with monotonic regression curve
#' plotCorrelation(obj, "hsa-miR-34a-5p", "PAX8", condition = "disease")
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
plotCorrelation <- function(mirnaObj,
                            mirna,
                            gene,
                            condition,
                            showCoeff = TRUE,
                            regression = TRUE,
                            useRanks = FALSE,
                            lineCol = "red",
                            lineType = "dashed",
                            lineWidth = 0.8,
                            pointSize = 3,
                            colorScale = NULL) {
  
  ## input checks
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (pairedSamples(mirnaObj) == FALSE) {
    stop(paste("Correlation analysis can only be performed for paired",
               "samples! See ?mirnaTargetsIntegration"), call. = FALSE)
  }
  if (!is.character(mirna) |
      length(mirna) != 1 |
      !mirna %in% rownames(mirnaObj[["microRNA"]])) {
    stop(paste("'mirna' must be a valid miRNA name that is present in",
               "miRNA expression matrix."),
         call. = FALSE)
  }
  if (!is.character(gene) |
      length(gene) != 1 |
      !gene %in% rownames(mirnaObj[["genes"]])) {
    stop(paste("'gene' must be a valid gene name that is present in",
               "gene expression matrix."),
         call. = FALSE)
  }
  if (length(condition) == 1) {
    if (!is.character(condition) |
        !condition %in%
        colnames(MultiAssayExperiment::colData(exampleObject))) {
      stop(paste("'condition' must be the column name of a variable specified",
                 "in the metadata (colData) of a MirnaExperiment object; or,",
                 "alternatively, it must be a character/factor object that",
                 "specifies group memberships."),
           call. = FALSE)
    }
  } else {
    if ((!is.character(condition) & !is.factor(condition)) |
        length(condition) != nrow(MultiAssayExperiment::colData(mirnaObj))) {
      stop(paste("'condition' must be the column name of a variable specified",
                 "in the metadata (colData) of a MirnaExperiment object; or,",
                 "alternatively, it must be a character/factor object that",
                 "specifies group memberships."),
           call. = FALSE)
    }
  }
  if (!is.logical(showCoeff) |
      length(showCoeff) != 1) {
    stop("'showCoeff' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.logical(regression) |
      length(regression) != 1) {
    stop("'regression' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.logical(useRanks) |
      length(useRanks) != 1) {
    stop("'useRanks' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.character(lineCol) |
      length(lineCol) != 1 |
      !lineCol %in% grDevices::colors()) {
    stop(paste("'lineCol' must be an R color name. All available colors",
               "can be listed with 'colors()'."),
         call. = FALSE)
  }
  if (!is.character(lineType) |
      length(lineType) != 1 |
      !lineType %in% c("blank", "solid", "dashed", "dotted", "dotdash",
                       "longdash", "twodash")) {
    stop(paste("'lineType' must be either 'blank', 'solid', 'dashed'",
               "(default), 'dotted', 'dotdash', 'longdash' or 'twodash'.",
               "For additional details see ?plotCorrelation"),
         call. = FALSE)
  }
  if (!is.numeric(lineWidth) |
      length(lineWidth) != 1 |
      lineWidth < 0) {
    stop("'lineWidth' must be a non-neagtive number! (default is 0.8)",
         call. = FALSE)
  }
  if (!is.numeric(pointSize) |
      length(pointSize) != 1 |
      pointSize < 0) {
    stop("'pointSize' must be a non-neagtive number! (default is 3)",
         call. = FALSE)
  }
  if (length(condition) == 1 & !is.null(colorScale)) {
    if ((!is.null(colorScale) & !is.character(colorScale)) |
        any(!colorScale %in% grDevices::colors()) |
        !identical(
          sort(names(colorScale)),
          sort(unique(MultiAssayExperiment::colData(mirnaObj)[, condition])))) {
      stop(paste("'colorScale' must be a named character vector where values",
                 "consist of R colors, whereas names coincide to the different",
                 "conditions. For additional details see ?plotCorrelation."),
           call. = FALSE)
    }
  } else if (length(condition) != 1 & !is.null(colorScale)) {
    if ((!is.null(colorScale) & !is.character(colorScale)) |
        any(!colorScale %in% grDevices::colors()) |
        !identical(sort(names(colorScale)),
                   as.character(sort(unique(condition))))) {
      stop(paste("'colorScale' must be a named character vector where values",
                 "consist of R colors, whereas names coincide to the different",
                 "conditions. For additional details see ?plotCorrelation."),
           call. = FALSE)
    }
  }
  
  ## get integration results
  intRes <- mirnaTargetsIntegration(mirnaObj)
  
  ## verify that correlation analysis has been performed
  if (colnames(intRes)[2] != "Target") {
    stop(paste("Correlation analysis must be performed before using this",
               "function! See ?mirnaTargetsIntegration"), call. = FALSE)
  }
  
  ## check if the specified miRNA-target pair is present
  if (nrow(intRes[intRes$microRNA == mirna &
                  intRes$Target == gene, ]) == 0) {
    stop("This miRNA-Target pair doesn't show any significant correlation!",
         call. = FALSE)
  }
  
  ## determine the correlation coefficient used
  corMethod <- colnames(mirnaTargetsIntegration(mirnaObj))[4]
  corMethod <- tolower(strsplit(corMethod, ".", fixed = TRUE)[[1]][1])
  
  ## inform the user about using Pearson's correlation with rank data
  if (useRanks == TRUE & corMethod == "pearson") {
    warning(paste("It's not possible to represent Pearson's correlation",
                  "using rank data! Ignoring 'useRanks'..."), call. = FALSE)
  }
  
  ## extract miRNA and gene expression values
  mirnaExpr <- mirnaObj[["microRNA"]]
  geneExpr <- mirnaObj[["genes"]]
  
  ## define condition vector
  if (is.character(condition) & length(condition) == 1) {
    cond <- MultiAssayExperiment::colData(mirnaObj)[, condition]
  } else if (is.factor(condition)) {
    cond <- as.character(condition)
  } else {
    cond <- condition
  }
  names(cond) <- MultiAssayExperiment::colData(mirnaObj)[, "primary"]
  
  ## check if samples are paired, otherwise exclude unpaired samples
  sMap <- MultiAssayExperiment::sampleMap(mirnaObj)
  mirnaSamples <- sMap$primary[sMap$assay == "microRNA"]
  geneSamples <- sMap$primary[sMap$assay == "genes"]
  
  if (!identical(mirnaSamples, geneSamples)) {
    
    ## determine common and uncommon samples
    common <- intersect(mirnaSamples, geneSamples)
    unpaired  <- setdiff(mirnaSamples, geneSamples)
    
    ## remove samples without measurments of both miRNAs and genes
    if (length(unpaired) > 0) {
      
      mirnaExpr <- mirnaExpr[, sMap$colname[sMap$assay == "microRNA" &
                                              sMap$primary %in% common]]
      geneExpr <- geneExpr[, sMap$colname[sMap$assay == "genes" &
                                            sMap$primary %in% common]]
      
    }
    
    ## order the columns of expression matrices in the same way
    mirnaMap = sMap[sMap$assay == "microRNA", ]
    geneMap = sMap[sMap$assay == "genes", ]
    mirnaOrder <- mirnaMap$primary[order(match(mirnaMap$colname,
                                               colnames(mirnaExpr)))]
    geneExpr <- geneExpr[, geneMap$colname[order(match(geneMap$primary,
                                                       mirnaOrder))]]
    
    ## re-define condition vector without removed samples
    cond <- cond[mirnaOrder]
    
    ## re-define colorScale without removed conditions
    colorScale <- colorScale[names(colorScale) %in% cond]
    
  }
  
  ## select the specified miRNA-target pair
  mirnaExpr <- mirnaExpr[mirna, ]
  geneExpr <- geneExpr[gene, ]
  
  ## convert to ranks if desired by the user
  if (useRanks == TRUE & corMethod != "pearson") {
    mirnaExpr <- rank(mirnaExpr)
    geneExpr <- rank(geneExpr)
  }
  
  ## create a dataframe with miRNA and gene expression
  corDf <- data.frame("miRNA" = mirnaExpr,
                      "gene" = geneExpr,
                      "Condition" = cond)
  
  ## define plot labels
  xlab <- paste(mirna, "expression")
  ylab <- paste(gene, "expression")
  
  ## define miRNA-mRNA direction
  selPair <- intRes[intRes$microRNA == mirna &
                      intRes$Target == gene, ]
  dir <- selPair$microRNA.Direction
  monoDir <- ifelse(dir == "upregulated", "decreasing", "increasing")
  coeffPos <- ifelse(dir == "upregulated", 1, 0)
  
  ## create correlation plot
  corPlot <- ggplot2::ggplot(corDf,
                             ggplot2::aes(miRNA, gene)) +
    ggplot2::geom_point(ggplot2::aes(color = Condition), size = pointSize)
  
  ## fit regression line/curve
  if (regression == TRUE) {
    if (useRanks == TRUE | corMethod == "pearson") {
      
      ## add regression line to the plot
      corPlot <- corPlot +
        ggplot2::geom_smooth(method = "lm",
                             formula = y ~ x,
                             se = FALSE, color = lineCol,
                             linewidth = lineWidth, linetype = lineType)
      
    } else {
      
      ## compute a monotonic regression curve
      x = mirnaExpr
      y = geneExpr
      monoFit <- MonoPoly::monpol(y ~ x,
                                  plot.it = FALSE,
                                  monotone = monoDir,
                                  algorithm = "Hawkins")
      
      ## create a sequence of values for x-axis
      xSeq <- seq(min(corDf$miRNA), max(corDf$miRNA), length.out = 100)
      
      ## predict the y values using the fitted model
      ySeq <- stats::predict(monoFit, newdata = data.frame(x = xSeq))
      
      ## create a data.frame for the fitted curve
      fitDf <- data.frame("xFit" = xSeq, "yFit" = ySeq[, "x"])
      
      ## add the monotonic curve to the plot
      corPlot <- corPlot +
        ggplot2::geom_line(data = fitDf,
                           ggplot2::aes(x = xFit, y = yFit),
                           color = lineCol,
                           linewidth = lineWidth,
                           linetype = lineType)
    }
  }
  
  ## add correlation coefficient to the plot
  if (showCoeff == TRUE) {
    
    ## determine coefficient symbol
    if (corMethod == "pearson") {
      corCoeff <- "R"
    } else if (corMethod == "spearman") {
      corCoeff <- "rho"
    } else if (corMethod == "kendall") {
      corCoeff <- "tau"
    }
    
    ## add correlation coefficient through ggpubr::stat_cor()
    corPlot <- corPlot +
      ggpubr::stat_cor(ggplot2::aes(label = ggplot2::after_stat(r.label)),
                       method = corMethod,
                       cor.coef.name = corCoeff,
                       label.x.npc = coeffPos, label.y.npc = 1,
                       hjust = coeffPos, vjust = 1)
    
  }
  
  ## add colorScale to ggplot2 graph
  if (!is.null(colorScale)) {
    corPlot <- corPlot +
      ggplot2::scale_color_manual(values = colorScale)
  }
  
  ## rename plot labels and set ggplot2 theme
  corPlot <- corPlot +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw()
  
  ## return the generated plot
  return(corPlot)
  
}

