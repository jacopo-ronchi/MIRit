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





#' Create a dotplot for functional enrichment analysis
#'
#' This function produces a dotplot to show the results of functional
#' enrichment analyses carried out through over-representation analysis (ORA),
#' gene set enrichment analysis (GSEA), and competitive gene set test accounting
#' for inter-gene correlation (CAMERA). In particular, this function can take
#' as input enrichment results generated by [enrichGenes()] and
#' [enrichMirnas()] functions.
#'
#' @param enrichment An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#' @param showTerms It is the number of top significant terms to show or,
#' alternatively, a character vector indicating the enriched terms to plot.
#' Default is `10`
#' @param splitDir Logical, if `TRUE` the resulting plot will be divided in
#' two columns on the basis of enrichment direction (Up and Down).
#' Default is `TRUE`. This only applies if enrichment method is GSEA or CAMERA
#' @param ordBy The parameter used to set the x-axis scale. It must be one of
#' `ratio` (default), `padj`, `pval` and `overlap`
#' @param sizeBy The parameter used to set the size scale. It must be one of
#' `ratio`, `padj`, `pval` and `overlap` (default)
#' @param colBy The parameter used to set the color scale. It must be one of
#' `ratio`, `padj` (default), `pval` and `overlap`
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' A `ggplot` graph with a dotplot of enrichment results.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform GSEA with KEGG database
#' gse <- enrichGenes(obj, method = "GSEA",
#' database = "KEGG", organism = "Homo sapiens")
#'
#' # extract results
#' res <- enrichmentResults(gse)
#'
#' # plot results
#' enrichmentDotplot(gse)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
enrichmentDotplot <- function(enrichment,
                              showTerms = 10,
                              splitDir = TRUE,
                              ordBy = "ratio",
                              sizeBy = "overlap",
                              colBy = "padj",
                              title = NULL) {
  
  ## check inputs
  if (!is(enrichment, "FunctionalEnrichment")) {
    stop("'enrichment' should be of class FunctionalEnrichment!",
         call. = FALSE)
  }
  if (nrow(enrichmentResults(enrichment)) < 1) {
    stop("'enrichment' object does not contain any significant results!",
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
      !all(showTerms %in% enrichmentResults(enrichment)$pathway)) {
    stop("The terms provided are not present in 'enrichment' pathway column!",
         call. = FALSE)
  }
  if (!is.logical(splitDir) |
      length(splitDir) != 1) {
    stop("'splitDir' must be logical (TRUE/FALSE)! See ?enrichmentDotplot",
         call. = FALSE)
  }
  if (!is.character(ordBy) |
      length(ordBy) != 1 |
      !ordBy %in% c("ratio", "padj", "pval", "overlap")) {
    stop(paste("'ordBy' must be one of: 'ratio' (default),",
               "'padj', 'pval', 'overlap'"),
         call. = FALSE)
  }
  if (!is.character(sizeBy) |
      length(sizeBy) != 1 |
      !sizeBy %in% c("ratio", "padj", "pval", "overlap")) {
    stop(paste("'sizeBy' must be one of: 'ratio', 'padj',",
               "'pval', 'overlap' (default)"),
         call. = FALSE)
  }
  if (!is.character(colBy) |
      length(colBy) != 1 |
      !colBy %in% c("ratio", "padj", "pval", "overlap")) {
    stop(paste("'colBy' must be one of: 'ratio', 'padj'",
               "(default), 'pval', 'overlap'"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot (e.g. 'Enrichment').",
               "For additional details see ?enrichmentDotplot"),
         call. = FALSE)
  }
  
  ## extract results from enrichment object
  res <- enrichmentResults(enrichment)
  
  ## compute gene ratio and direction for different analyses
  if (enrichmentMethod(enrichment) == "Gene-Set Enrichment Analysis (GSEA)") {
    ov <- as.numeric(lapply(res$leadingEdge, length))
    res$overlap <- ov
    res$ratio <- ov/res$size
    res$direction <- "Up"
    res$direction[which(res$NES < 0)] <- "Down"
  } else {
    res$ratio <- res$overlap/res$size
  }
  
  ## order results based on padj
  res <- res[order(res$padj), ]
  
  ## select terms to be shown in the dotplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$pathway %in% showTerms), ]
  }
  
  ## set an x-axis label for 'ratio'
  if (ordBy == "ratio") {
    ordLabel <- "Gene-Set Overlap"
  } else {
    ordLabel <- ordBy
  }
  
  ## reformat results based on specified criterion
  if (ordBy != "padj" & ordBy != "pval") {
    res <- res[order(res[, ordBy], decreasing = TRUE), ]
  } else {
    res[, ordBy] <- -log10(res[, ordBy])
    res <- res[order(res[, ordBy], decreasing = TRUE), ]
    ordLabel <- paste("-log10(", ordBy, ")", sep = "")
  }
  
  ## create a dotplot
  dotRes <- ggplot2::ggplot(res,
                            ggplot2::aes(x = !!ggplot2::sym(ordBy),
                                         y = stats::reorder(.data$pathway,
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
  if (splitDir == TRUE &
      enrichmentMethod(enrichment) != "Over-Representation Analysis (ORA)") {
    dotRes <- dotRes +
      ggplot2::facet_grid(~ direction) +
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





#' Create a barplot for functional enrichment analysis
#'
#' This function produces a barplot to show the results of functional
#' enrichment analyses carried out through over-representation analysis (ORA),
#' gene set enrichment analysis (GSEA), and competitive gene set test accounting
#' for inter-gene correlation (CAMERA). In particular, this function can take
#' as input enrichment results generated by [enrichGenes()] and
#' [enrichMirnas()] functions.
#'
#' @param enrichment An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#' @param showTerms It is the number of top significant terms to show or,
#' alternatively, a character vector indicating the enriched terms to plot.
#' Default is `10`
#' @param splitDir Logical, if `TRUE` the resulting plot will be divided in
#' two columns on the basis of enrichment direction (Up and Down).
#' Default is `TRUE`. This only applies if enrichment method is GSEA or CAMERA
#' @param ordBy The parameter used to set the x-axis scale. It must be one of
#' `ratio` (default), `padj`, `pval` and `overlap`
#' @param colBy The parameter used to set the color scale. It must be one of
#' `ratio`, `padj` (default), `pval` and `overlap`
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' A `ggplot` graph with a barplot of enrichment results.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform GSEA with KEGG database
#' gse <- enrichGenes(obj, method = "GSEA",
#' database = "KEGG", organism = "Homo sapiens")
#'
#' # extract results
#' res <- enrichmentResults(gse)
#'
#' # plot results
#' enrichmentBarplot(gse)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
enrichmentBarplot <- function(enrichment,
                              showTerms = 10,
                              splitDir = TRUE,
                              ordBy = "ratio",
                              colBy = "padj",
                              title = NULL) {
  
  ## check inputs
  if (!is(enrichment, "FunctionalEnrichment")) {
    stop("'enrichment' should be of class FunctionalEnrichment!",
         call. = FALSE)
  }
  if (nrow(enrichmentResults(enrichment)) < 1) {
    stop("'enrichment' object does not contain any significant results!",
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
      !all(showTerms %in% enrichmentResults(enrichment)$pathway)) {
    stop("The terms provided are not present in 'enrichment' pathway column!",
         call. = FALSE)
  }
  if (!is.logical(splitDir) |
      length(splitDir) != 1) {
    stop("'splitDir' must be logical (TRUE/FALSE)! See ?enrichmentBarplot",
         call. = FALSE)
  }
  if (!is.character(ordBy) |
      length(ordBy) != 1 |
      !ordBy %in% c("ratio", "padj", "pval", "overlap")) {
    stop(paste("'ordBy' must be one of: 'ratio' (default),",
               "'padj', 'pval', 'overlap'"),
         call. = FALSE)
  }
  if (!is.character(colBy) |
      length(colBy) != 1 |
      !colBy %in% c("ratio", "padj", "pval", "overlap")) {
    stop(paste("'colBy' must be one of: 'ratio', 'padj'",
               "(default), 'pval', 'overlap'"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot (e.g. 'Enrichment').",
               "For additional details see ?enrichmentBarplot"),
         call. = FALSE)
  }
  
  ## extract results from enrichment object
  res <- enrichmentResults(enrichment)
  
  ## compute gene ratio and direction for different analyses
  if (enrichmentMethod(enrichment) == "Gene-Set Enrichment Analysis (GSEA)") {
    ov <- as.numeric(lapply(res$leadingEdge, length))
    res$overlap <- ov
    res$ratio <- ov/res$size
    res$direction <- "Up"
    res$direction[which(res$NES < 0)] <- "Down"
  } else {
    res$ratio <- res$overlap/res$size
  }
  
  ## order results based on padj
  res <- res[order(res$padj), ]
  
  ## select terms to be shown in the barplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$pathway %in% showTerms), ]
  }
  
  ## set an x-axis label for 'ratio'
  if (ordBy == "ratio") {
    ordLabel <- "Gene-Set Overlap"
  } else {
    ordLabel <- ordBy
  }
  
  ## reformat results based on specified criterion
  if (ordBy != "padj" & ordBy != "pval") {
    res <- res[order(res[, ordBy], decreasing = TRUE), ]
  } else {
    res[, ordBy] <- -log10(res[, ordBy])
    res <- res[order(res[, ordBy], decreasing = TRUE), ]
    ordLabel <- paste("-log10(", ordBy, ")", sep = "")
  }
  
  ## create a dotplot
  barRes <- ggplot2::ggplot(res,
                            ggplot2::aes(x = !!ggplot2::sym(ordBy),
                                         y = stats::reorder(.data$pathway,
                                                            !!ggplot2::sym(ordBy)),
                                         fill = !!ggplot2::sym(colBy))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient(low = "red", high = "blue",
                                 guide = ggplot2::guide_colorbar(reverse = TRUE)) +
    ggplot2::ylab(NULL) +
    ggplot2::scale_y_discrete() +
    ggplot2::xlab(ordLabel) +
    theme_enr()
  
  ## divide by enrichment direction
  if (splitDir == TRUE &
      enrichmentMethod(enrichment) != "Over-Representation Analysis (ORA)") {
    barRes <- barRes +
      ggplot2::facet_grid(~ direction) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 12))
  }
  
  ## add the title of the plot
  if (!is.null(title)) {
    barRes <- barRes +
      ggplot2::ggtitle(title)
  }
  
  ## return ggplot2 graph
  return(barRes)
  
}





#' Create a ridgeplot to display the results of GSEA analysis
#'
#' This function creates a ridgeplot that is useful for showing the results
#' of GSEA analyses. The output of this function is a plot where enriched
#' terms/pathways found with [enrichGenes()] function are visualized on the
#' basis of the ranking metric used for the analysis. The resulting
#' areas represent the density of signed p-values, log2 fold changes, or
#' log.p-values belonging to genes annotated to that category.
#'
#' @param enrichment An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#' @param showTerms It is the number of top significant terms to show or,
#' alternatively, a character vector indicating the enriched terms to plot.
#' Default is `10`
#' @param colBy The parameter used to set the color scale. It must be one of
#' `ratio`, `padj` (default), `pval` and `overlap`
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' An object of class `ggplot` containing the ridgeplot of GSEA results.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform GSEA with KEGG database
#' gse <- enrichGenes(obj, method = "GSEA",
#' database = "KEGG", organism = "Homo sapiens")
#'
#' # extract results
#' res <- enrichmentResults(gse)
#'
#' # plot results as a ridgeplot
#' gseaRidgeplot(gse)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @importFrom rlang .data
#' @export
gseaRidgeplot <- function(enrichment,
                          showTerms = 10,
                          colBy = "padj",
                          title = NULL) {
  
  if (!is(enrichment, "FunctionalEnrichment")) {
    stop("'enrichment' should be of class FunctionalEnrichment!",
         call. = FALSE)
  }
  if (nrow(enrichmentResults(enrichment)) < 1) {
    stop("'enrichment' object does not contain any significant results!",
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
      !all(showTerms %in% enrichmentResults(enrichment)$pathway)) {
    stop("The terms provided are not present in 'enrichment' pathway column!",
         call. = FALSE)
  }
  if (!is.character(colBy) |
      length(colBy) != 1 |
      !colBy %in% c("ratio", "padj", "pval", "overlap")) {
    stop(paste("'colBy' must be one of: 'ratio', 'padj'",
               "(default), 'pval', 'overlap'"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot (e.g. 'Enrichment').",
               "For additional details see ?enrichmentDotplot"),
         call. = FALSE)
  }
  
  ## check that GSEA has been performed
  if (enrichmentMethod(enrichment) != "Gene-Set Enrichment Analysis (GSEA)") {
    stop(paste("To use this function, 'enrichment' must originate from a",
               "gene-set enrichment analysis (GSEA).",
               "For additional details see ?enrichmentBarplot"),
         call. = FALSE)
  }
  
  ## extract results from enrichment object
  res <- enrichmentResults(enrichment)
  
  ## extract ranked list from FunctionalEnrichment object
  rankedList <- enrichmentMetric(enrichment)
  names(rankedList) <- enrichedFeatures(enrichment)
  
  ## compute gene ratio
  ov <- as.numeric(lapply(res$leadingEdge, length))
  res$overlap <- ov
  res$ratio <- ov/res$size
  
  ## order results based on padj
  res <- res[order(res$padj), ]
  
  ## select terms to be shown in the ridgeplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$pathway %in% showTerms), ]
  }
  
  ## retrieve ranking metric for genes in each category
  metricSet <- lapply(res[, "leadingEdge"], function(x) {
    rankedList[intersect(x, names(rankedList))]
  })
  names(metricSet) <- res$pathway
  
  ## create a data.frame with metric reported for each category
  setLen <- vapply(metricSet, length, FUN.VALUE = numeric(1))
  meanMetric <- vapply(metricSet, mean, FUN.VALUE = numeric(1))
  ridgeDf <- data.frame(pathway = rep(res$pathway, setLen),
                        padj = rep(res$padj, setLen),
                        pval = rep(res$pval, setLen),
                        overlap = rep(res$overlap, setLen),
                        metric = rep(meanMetric, setLen),
                        val = unlist(metricSet))
  
  ## create a ridgeplot
  ridPlot <- ggplot2::ggplot(ridgeDf,
                             ggplot2::aes(x = .data$val,
                                          y = stats::reorder(.data$pathway,
                                                             .data$metric),
                                          fill = !!ggplot2::sym(colBy))) +
    ggridges::geom_density_ridges() +
    ggplot2::scale_fill_continuous(low = "red",
                                   high = "blue",
                                   name = colBy,
                                   guide = ggplot2::guide_colorbar(reverse = TRUE)) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("Ranking Metric") +
    theme_enr()
  
  ## add the title of the plot
  if (!is.null(title)) {
    ridPlot <- ridPlot +
      ggplot2::ggtitle(title)
  }
  
  ## return ggplot2 graph
  return(ridPlot)
  
}





#' Create a GSEA plot that displays the running enrichment score (ES) for a
#' given pathway
#'
#' This function creates a classic enrichment plot to show the results of
#' gene set enrichment analyses (GSEA). In particular, this function takes as
#' input GSEA results originating from the [enrichGenes()] function, and
#' returns a [ggplot2] object with GSEA plot. In this kind of plots, the
#' running enrichment score (ES) for a given pathway is shown on the y-axis,
#' whereas gene positions in the ranked list are reported on the x-axis.
#'
#' @param enrichment An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#' @param pathway It must be the name of a significantly enriched term/pathway
#' for which we want to produce a GSEA plot (e.g. 'Thyroid hormone synthesis') 
#' @param showTitle Logical, whether to add the name of the pathway/term as
#' plot title. Default is TRUE
#' @param rankingMetric Logical, whether to show the variations of the ranking
#' metric below the plot. Default is FALSE
#' @param lineColor It must be an R color name that specifies the color of
#' the running score line. Default is `green`. All available colors
#' can be listed with [grDevices::colors()]
#' @param lineSize The line width of the running score line. Default is `1`
#' @param vlineColor It must be an R color name that specifies the color of
#' the vertical line indicating the enrichment score (ES). Default is `red`.
#' All available colors can be listed with [grDevices::colors()]
#' @param vlineSize The line width of the vertical line indicating the
#' enrichment score (ES). Default is `0.6`
#'
#' @returns
#' An object of class `ggplot` containing the GSEA plot.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform GSEA with KEGG database
#' gse <- enrichGenes(obj, method = "GSEA",
#' database = "KEGG", organism = "Homo sapiens")
#'
#' # extract results
#' res <- enrichmentResults(gse)
#'
#' # plot results
#' gseaPlot(gse, pathway = "Thyroid hormone synthesis")
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
gseaPlot <- function(enrichment,
                     pathway,
                     showTitle = TRUE,
                     rankingMetric = FALSE,
                     lineColor = "green",
                     lineSize = 1,
                     vlineColor = "red",
                     vlineSize = 0.6) {
  
  ## check inputs
  if (!is(enrichment, "FunctionalEnrichment")) {
    stop("'enrichment' should be of class FunctionalEnrichment!",
         call. = FALSE)
  }
  if (nrow(enrichmentResults(enrichment)) < 1) {
    stop("'enrichment' object does not contain any significant results!",
         call. = FALSE)
  }
  if (is.character(pathway) &
      !pathway %in% enrichmentResults(enrichment)$pathway) {
    stop("'pathway' is not present in 'enrichment' pathway column!",
         call. = FALSE)
  }
  if (!is.logical(showTitle) |
      length(showTitle) != 1) {
    stop("'showTitle' must be logical (TRUE/FALSE)! See ?gseaPlot",
         call. = FALSE)
  }
  if (!is.logical(rankingMetric) |
      length(rankingMetric) != 1) {
    stop("'rankingMetric' must be logical (TRUE/FALSE)! See ?gseaPlot",
         call. = FALSE)
  }
  if (!is.character(lineColor) |
      length(lineColor) != 1 |
      !lineColor %in% grDevices::colors()) {
    stop(paste("'lineColor' must be an R color name. All available colors",
               "can be listed with 'colors()'."),
         call. = FALSE)
  }
  if (!is.numeric(lineSize) |
      length(lineSize) != 1 |
      lineSize < 0) {
    stop("'lineSize' must be a non-neagtive number! (default is 1)",
         call. = FALSE)
  }
  if (!is.character(vlineColor) |
      length(vlineColor) != 1 |
      !vlineColor %in% grDevices::colors()) {
    stop(paste("'vlineColor' must be an R color name. All available colors",
               "can be listed with 'colors()'."),
         call. = FALSE)
  }
  if (!is.numeric(vlineSize) |
      length(vlineSize) != 1 |
      vlineSize < 0) {
    stop("'vlineSize' must be a non-neagtive number! (default is 0.6)",
         call. = FALSE)
  }
  
  ## check that GSEA has been performed
  if (enrichmentMethod(enrichment) != "Gene-Set Enrichment Analysis (GSEA)") {
    stop(paste("To use this function, 'enrichment' must originate from a",
               "gene-set enrichment analysis (GSEA).",
               "For additional details see ?enrichmentBarplot"),
         call. = FALSE)
  }
  
  ## extract ranked list and gene sets from FunctionalEnrichment object
  geneSet <- enrichment@geneSet
  rankedList <- enrichmentMetric(enrichment)
  names(rankedList) <- enrichedFeatures(enrichment)
  
  ## extract the gene set of interest
  gs <- geneSet[[pathway]]
  
  ## keep only the genes present in ranked list
  gs <- intersect(gs, names(rankedList))
  
  ## define the lengths of gene set and ranked list
  lengthRank <- length(rankedList)
  lengthSet <- length(gs)
  
  ## create vectors for positive and negative hits
  pos <- numeric(lengthRank)
  misses <- numeric(lengthRank)
  
  ## identify the positions of the genes in ranked list belonging to the pathway
  hits <- names(rankedList) %in% gs
  
  ## calculate scores
  pos[hits] <- abs(rankedList[hits])
  misses[!hits] <-  1/(lengthRank-lengthSet)
  
  ## determine cumulative sums
  pSum <- sum(pos)
  pos <- cumsum(pos/pSum)
  misses <- cumsum(misses)
  
  ## compute GSEA running score
  runScore <- pos - misses
  
  ## set the maximum absolute deviance from zero (ES)
  if(abs(max(runScore)) < abs(min(runScore))) {
    es <- min(runScore)
  } else {
    es <- max(runScore)
  }
  
  ## create data.frame for ggplot2
  gseData <- data.frame(position = seq(length(runScore)),
                        runScore = runScore,
                        hits = hits,
                        genes = names(rankedList),
                        metric = rankedList)
  
  ## add ymin and ymax for ranges
  gseData$ymin <- 0
  gseData$ymax <- 0
  gseData$ymin[hits] <- - diff(range(runScore))/20
  gseData$ymax[hits] <- diff(range(runScore))/20
  
  ## create GSEA plot with ggplot2
  gsePlot <- ggplot2::ggplot(gseData, ggplot2::aes(x = position)) +
    ggplot2::geom_line(ggplot2::aes(y = runScore),
                       linewidth = lineSize,
                       color = lineColor) +
    ggplot2::geom_vline(xintercept = which(runScore == es),
                        color = vlineColor,
                        linewidth = vlineSize,
                        linetype = "dashed") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = ymin,
                                         ymax = ymax)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::xlab("Position in Ranked List of Genes") +
    ggplot2::ylab("Running Score") +
    theme_enr()
  
  ## add pathway title if desired
  if (showTitle == TRUE) {
    gsePlot <- gsePlot +
      ggplot2::ggtitle(pathway)
  }
  
  ## show ranking metric if wanted by the user
  if (rankingMetric == TRUE) {
    
    ## remove x axis element from GSEA plot
    gsePlot <- gsePlot +
      ggplot2::xlab(NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    
    ## create rank plot
    rankPlot <- ggplot2::ggplot(gseData, ggplot2::aes(x = position)) +
      ggplot2::geom_segment(ggplot2::aes(xend = position,
                                         y = metric,
                                         yend=0),
                            color="black") +
      ggplot2::ylab("Ranking Metric") +
      ggplot2::xlab("Position in Ranked List of Genes") +
      theme_enr()
    
    ## merge GSEA plot and rank plot
    gsePlot <- ggpubr::ggarrange(gsePlot,
                                 rankPlot,
                                 ncol = 1,
                                 heights = c(1, 0.7),
                                 align = "v")
    
  }
  
  ## return ggplot2 object
  return(gsePlot)
  
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





#' Create a trackplot to show the association between miRNAs and disease-SNPs
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
#' analysis performed through the [mirnaIntegration()] function, to
#' graphically visualize the quantitative effect of miRNA dysregulations on
#' target gene expression. Furthermore, this function performs linear/monotonic
#' regression to better represent the relationships between miRNA-target pairs.
#' 
#' When non-parametric correlation has been performed with the
#' [mirnaIntegration()] function, a regression line can be fitted through
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
#' @param lineType It specifies the line type used for the regression line. It
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
#' obj <- mirnaIntegration(obj)
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
  if (integratedTargets == TRUE &
      max(dim(integration(mirnaObj))) == 0) {
    stop(paste("Integration analysis is not detected in 'mirnaObj'!",
               "Before using this function, expression levels of miRNAs and",
               "genes must be integrated with the 'mirnaIntegration()'",
               "function. See '?mirnaIntegration' for details."),
         call. = FALSE)
  }
  if (pairedSamples(mirnaObj) == FALSE) {
    stop(paste("Correlation analysis can only be performed for paired",
               "samples! See ?mirnaIntegration"), call. = FALSE)
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
        !(condition %in% colnames(MultiAssayExperiment::colData(mirnaObj)) &
          !condition %in% c("primary", "mirnaCol", "geneCol"))) {
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
  intRes <- integration(mirnaObj)
  
  ## verify that correlation analysis has been performed
  if (colnames(intRes)[2] != "Target") {
    stop(paste("Correlation analysis must be performed before using this",
               "function! See ?integration"), call. = FALSE)
  }
  
  ## check if the specified miRNA-target pair is present
  if (nrow(intRes[intRes$microRNA == mirna &
                  intRes$Target == gene, ]) == 0) {
    stop("This miRNA-Target pair doesn't show any significant correlation!",
         call. = FALSE)
  }
  
  ## determine the correlation coefficient used
  corMethod <- integration(mirnaObj, param = TRUE)$method
  corMethod <- tolower(sub("([A-Za-z]+).*", "\\1", corMethod))
  
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





#' Represent differentially expressed miRNAs/genes as boxplots, barplots or
#' violinplots
#' 
#' This function is able to produce boxplots, barplots and violinplots that are
#' useful to visualize miRNA and gene differential expression. The user just
#' has to provide a vector of interesting miRNA/genes that he wants to plot
#' (e.g. "hsa-miR-34a-5p", "hsa-miR-146b-5p", "PAX8"). The chart type can be
#' specified through the `graph` parameter.
#' 
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param features A character vector containing the genes/miRNAs to plot
#' @param condition It must be the column name of a variable specified in the
#' metadata (colData) of a [`MirnaExperiment`][MirnaExperiment-class] object;
#' or, alternatively, it must be a character/factor object that specifies
#' group memberships (eg. c("healthy, "healthy", "disease", "disease"))
#' @param graph The type of plot to produce. It must be one of `boxplot`
#' (default), `barplot`, `violinplot`
#' @param showSignificance Logical, whether to display statistical significance
#' or not. Default is TRUE
#' @param starSig Logical, whether to represent statistical significance through
#' stars. Default is TRUE, and the significance scale is: \* for $p < 0.05$, \**
#' for $p < 0.01$, \*** for $p < 0.001$, and \**** for $p < 0.0001$. If
#' `starSig` is set to FALSE, p-values or adjusted p-values will be reported on
#' the plot as numbers
#' @param pCol The statistics used to evaluate comparison significance. It must
#' be one of `P.Value`, to use unadjusted p-values, and `adj.P.Val` (default),
#' to use p-values corrected for multiple testing
#' @param sigOffset The distance between the different brackets used to
#' show statistical significance. Default is 0.9, but the user can increment it
#' to enlarge the distance between significance brackets
#' @param sigLabelSize The size for the labels used to show statistical
#' significance. Default is 7, which is well suited for representing p-values
#' as significance stars. However, if `starSig` is set to FALSE, the user might
#' have to downsize this parameter
#' @param digits The number of digits to show when p-values are reported as
#' numbers (when `starSig` is FALSE). Default is 3
#' @param colorScale It must be a named character vector where values
#' correspond to R colors, while names coincide with the genes/miRNAs specified
#' in the `condition` parameter (eg. c("hsa-miR-34a-5p" = "blue",
#' "PAX8" = "red")). Default is NULL, in order to use the default color scale
#'
#' @returns
#' An object of class `ggplot` containing the plot.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # produce a boxplot for PAX8 and miR-34a-5p
#' plotDE(obj, features = c("hsa-miR-34a-5p", "PAX8"), condition = "disease")
#' 
#' # produce a barplot for PAX8 and miR-34a-5p without significance
#' plotDE(obj, features = c("hsa-miR-34a-5p", "PAX8"), condition = "disease",
#' graph = "barplot", showSignificance = FALSE)
#' 
#' # produce a violinplot for BCL2
#' plotDE(obj, features = "BCL2", condition = "disease", graph = "violinplot")
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @importFrom ggpubr mean_sd
#' @export
plotDE <- function(mirnaObj,
                   features,
                   condition,
                   graph = "boxplot",
                   showSignificance = TRUE,
                   starSig = TRUE,
                   pCol = "adj.P.Val",
                   sigOffset = 0.9,
                   sigLabelSize = 7,
                   digits = 3,
                   colorScale = NULL) {
  
  ## input checks
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (nrow(mirnaDE(mirnaObj, onlySignificant = FALSE)) == 0) {
    stop(paste("MiRNA differential expression results are not present in",
               "'mirnaObj'. Please, use 'performMirnaDE()' before using",
               "this function. See ?performMirnaDE"), call. = FALSE)
  }
  if (nrow(geneDE(mirnaObj, onlySignificant = FALSE)) == 0) {
    stop(paste("Gene differential expression results are not present in",
               "'mirnaObj'. Please, use 'performGeneDE()' before using",
               "this function. See ?performGeneDE"), call. = FALSE)
  }
  if (!is.character(features) |
      length(features) > 10 |
      any(!features %in% rownames(mirnaObj[["microRNA"]]) &
          !features %in% rownames(mirnaObj[["genes"]]))) {
    stop(paste("'features' must be a character vector containing miRNA and/or",
               "gene symbols (e.g. c('hsa-miR-34a-5p', 'PAX8').",
               "For additional details see ?plotDE"),
         call. = FALSE)
  }
  if (length(condition) == 1) {
    if (!is.character(condition) |
        !(condition %in% colnames(MultiAssayExperiment::colData(mirnaObj)) &
          !condition %in% c("primary", "mirnaCol", "geneCol"))) {
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
  if (!is.character(graph) |
      length(graph) != 1 |
      !graph %in% c("boxplot", "barplot", "violinplot")) {
    stop(paste("'graph' must be either 'boxplot' (default), `barplot`",
               "or 'violinplot'. For additional details see ?plotDE"),
         call. = FALSE)
  }
  if (!is.logical(showSignificance) |
      length(showSignificance) != 1) {
    stop("'showSignificance' must be logical (TRUE/FALSE)! See ?plotDE",
         call. = FALSE)
  }
  if (!is.logical(starSig) |
      length(starSig) != 1) {
    stop("'starSig' must be logical (TRUE/FALSE)! See ?plotDE",
         call. = FALSE)
  }
  if (!is.character(pCol) |
      length(pCol) != 1 |
      !pCol %in% c("P.Value", "adj.P.Val")) {
    stop(paste("'pCol' must be either 'P.Value' or 'adj.P.Val' (default).",
               "For additional details see ?plotDE"),
         call. = FALSE)
  }
  if (!is.numeric(sigOffset) |
      length(sigOffset) != 1 |
      sigOffset < 0) {
    stop("'sigOffset' must be a non-neagtive number! (default is 0.9)",
         call. = FALSE)
  }
  if (!is.numeric(sigLabelSize) |
      length(sigLabelSize) != 1 |
      sigLabelSize < 0) {
    stop("'sigLabelSize' must be a non-neagtive number! (default is 7)",
         call. = FALSE)
  }
  if (!is.numeric(digits) |
      length(digits) != 1 |
      digits < 0 |
      !digits%%1==0) {
    stop("'digits' must be a non-neagtive integer! (default is 3)",
         call. = FALSE)
  }
  if (!is.null(colorScale)) {
    if (!is.character(colorScale) |
        any(!colorScale %in% grDevices::colors()) |
        !identical(sort(names(colorScale)),
                   sort(features))) {
      stop(paste("'colorScale' must be a named character vector where values",
                 "consist of R colors, whereas names coincide to the different",
                 "genes specified in 'features'. For additional details",
                 "see ?plotDE"),
           call. = FALSE)
    }
  }
  
  ## extract miRNA and gene expression values
  mirnaExpr <- mirnaObj[["microRNA"]]
  geneExpr <- mirnaObj[["genes"]]
  
  ## load differential expression results
  statTest <- rbind(mirnaDE(mirnaObj, onlySignificant = FALSE),
                    geneDE(mirnaObj, onlySignificant = FALSE))
  
  ## define condition vector
  if (is.character(condition) & length(condition) == 1) {
    cond <- MultiAssayExperiment::colData(mirnaObj)[, condition]
  } else if (is.factor(condition)) {
    cond <- as.character(condition)
  } else {
    cond <- condition
  }
  names(cond) <- MultiAssayExperiment::colData(mirnaObj)[, "primary"]
  
  ## create a dataframe with miRNA and gene expression
  exprDf <- data.frame("Expression" = numeric(),
                       "Gene" = character(),
                       "Condition" = character())
  
  ## add entries to this data.frame for each gene/miRNA
  for (gene in features) {
    
    ## retrieve feature expression
    m <- mirnaExpr[rownames(mirnaExpr) == gene, ]
    g <- geneExpr[rownames(geneExpr) == gene, ]
    if (is.null(nrow(m))) {
      featExpr <- m
    } else {
      featExpr <- g
    }
    
    ## return feature expression, name and condition
    newDf <- data.frame("Expression" = featExpr,
                        "Gene" = rep(gene, length(cond)),
                        "Condition" = cond)
    exprDf <- rbind(exprDf, newDf)
    
  }
  
  ## restrict differential expression to the selected miRNAs/genes
  statTest <- statTest[statTest$ID %in% features, ]
  
  ## add conditions to differential expression results
  statTest$group1 <- unique(cond)[1]
  statTest$group2 <- unique(cond)[2]
  statTest$Gene <- features
  
  ## add y position for p-value labels
  yCo <- max(exprDf$Expression) + 0.1*max(exprDf$Expression)
  statTest$y.position <- seq(yCo,
                             yCo +
                               sigOffset*(length(unique(exprDf$Gene)) - 1),
                             sigOffset)
  
  ## round p-values if plotting numbers
  if (starSig == FALSE) {
    statTest[, pCol] <- round(statTest[, pCol], digits = digits)
  }
  
  ## use stars to show statistical significance
  if (starSig == TRUE) {
    statTest$star <- "ns"
    statTest$star[statTest[, pCol] < 0.05] <- "*"
    statTest$star[statTest[, pCol] < 0.01] <- "**"
    statTest$star[statTest[, pCol] < 0.001] <- "***"
    statTest$star[statTest[, pCol] < 0.0001] <- "****"
    pCol <- "star"
  }
  
  ## produce the desired plot
  if (graph == "boxplot") {
    
    ## create a grouped boxplot
    dePlot <- ggpubr::ggboxplot(data = exprDf, x = "Condition",
                                y = "Expression", fill = "Gene") +
      ggplot2::ylab(expression(paste(log[2], " expression")))
    
  } else if (graph == "barplot") {
    
    ## create a grouped barplot
    dePlot <- ggpubr::ggbarplot(data = exprDf,
                                x = "Condition",
                                y = "Expression",
                                fill = "Gene",
                                position = ggplot2::position_dodge(0.8),
                                add = "mean_sd",
                                error.plot = "upper_errorbar") +
      ggplot2::ylab(expression(paste(log[2], " expression")))
    
  } else if (graph == "violinplot") {
    
    ## create a grouped violinplot
    dePlot <- ggpubr::ggviolin(data = exprDf,
                               x = "Condition",
                               y = "Expression",
                               fill = "Gene",
                               add = "boxplot")
    
  }
  
  ## add significance levels
  if (showSignificance == TRUE) {
    
    dePlot <- dePlot +
      ggpubr::stat_pvalue_manual(data = statTest, label = pCol,
                                 step.increase = 0.1, size = sigLabelSize,
                                 step.group.by = "Gene", color = "Gene")
    
  }
  
  ## add colorScale to ggplot2 graph
  if (!is.null(colorScale)) {
    dePlot <- dePlot +
      ggplot2::scale_color_manual(values = colorScale) +
      ggplot2::scale_fill_manual(values = colorScale)
  }
  
  ## return the plot object
  return(dePlot)
  
}





#' Produce volcano plots to display miRNA/gene differential expression
#' 
#' This function allows the user to create publication-quality volcano plots to
#' represent the results of miRNA/gene differential expression. In this kind of
#' plots, the x-axis is relative to the log2 fold change between biological
#' conditions, while the y-axis contains the negative base-10 logarithm of the
#' p-value. Note, that even if volcano plots display unadjusted p-values on the
#' y-axis, the cutoff level shown in this plot derive from the adjusted p-value
#' cutoff used for differential expression analysis.
#' 
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param assay The results to display. It must be either 'microRNA', to plot
#' miRNA differential expression, or 'genes', to show the results for genes
#' @param labels The labels to show on the graph. Default is NULL not to
#' include labels. This parameter can be a character vector containing the IDs
#' of the features that you want to display. Alternatively, this parameter can
#' also be the number of most significant features for which we want to plot
#' labels
#' @param boxedLabel Logical, whether to show labels inside a rectangular shape
#' (default) or just as text elements
#' @param pointSize The size of points in the volcano plot (default is 3)
#' @param pointAlpha The transparency of points in the volcano plot (default
#' is 0.4)
#' @param interceptWidth The width of cutoff intercepts (default is 0.6)
#' @param interceptColor It must be an R color name that specifies the color of
#' cutoff intercepts. Default is `black`. All available colors can be listed
#' with [grDevices::colors()]
#' @param interceptType It specifies the line type used for cutoff intercepts.
#' It must be either 'blank', 'solid', 'dashed' (default), 'dotted', 'dotdash',
#' 'longdash' or 'twodash'
#' @param borderWidth The width of plot borders (default is 1)
#' @param colorScale It must be a character vector of length 3 containing valid
#' R color names for downregulated, non significant, and upregulated features,
#' respectively. Default value is `c('blue', 'grey', 'red')`. All available
#' colors can be listed with [grDevices::colors()]
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' An object of class `ggplot` containing the plot.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # produce a volcano plot for miRNAs with labels
#' plotVolcano(obj, "microRNA", labels = 5)
#' 
#' # produce a volcano plot for genes
#' plotVolcano(obj, "genes")
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
plotVolcano <- function(mirnaObj,
                        assay,
                        labels = NULL,
                        boxedLabel = TRUE,
                        pointSize = 3,
                        pointAlpha = 0.4,
                        interceptWidth = 0.6,
                        interceptColor = "black",
                        interceptType = "dashed",
                        borderWidth = 1,
                        colorScale = c("blue", "grey","red"),
                        title = NULL) {
  
  ## input checks
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (!assay %in% c("microRNA", "genes")) {
    stop("'assay' must be either 'microRNA' or 'genes'!", call. = FALSE)
  }
  if (nrow(mirnaDE(mirnaObj, onlySignificant = FALSE)) == 0 &
      assay == "microRNA") {
    stop(paste("MiRNA differential expression results are not present in",
               "'mirnaObj'. Please, use 'performMirnaDE()' before using",
               "this function. See ?performMirnaDE"), call. = FALSE)
  }
  if (nrow(geneDE(mirnaObj, onlySignificant = FALSE)) == 0 &
      assay == "genes") {
    stop(paste("Gene differential expression results are not present in",
               "'mirnaObj'. Please, use 'performGeneDE()' before using",
               "this function. See ?performGeneDE"), call. = FALSE)
  }
  if (!is.null(labels)) {
    if (!is.numeric(labels) &
        !is.character(labels)) {
      stop(paste("'labels' must be a character vector indicating the features",
                 "to plot; or, alternatively, it must be a number indicating",
                 "the top significant features to show. Default is NULL, not",
                 "to include labels. For details see ?plotVolcano"),
           call. = FALSE)
    }
  }
  if (!is.logical(boxedLabel) |
      length(boxedLabel) != 1) {
    stop("'boxedLabel' must be logical (TRUE/FALSE)! See ?plotVolcano",
         call. = FALSE)
  }
  if (!is.numeric(pointSize) |
      length(pointSize) != 1 |
      pointSize < 0) {
    stop("'pointSize' must be a non-neagtive number! (default is 3)",
         call. = FALSE)
  }
  if (!is.numeric(pointAlpha) |
      length(pointAlpha) != 1 |
      pointAlpha < 0 |
      pointAlpha > 1) {
    stop("'pointAlpha' must be a number between 0 and 1! (default is 0.4)",
         call. = FALSE)
  }
  if (!is.numeric(interceptWidth) |
      length(interceptWidth) != 1 |
      interceptWidth < 0) {
    stop("'interceptWidth' must be a non-neagtive number! (default is 0.6)",
         call. = FALSE)
  }
  if (!is.character(interceptColor) |
      length(interceptColor) != 1 |
      !interceptColor %in% grDevices::colors()) {
    stop(paste("'interceptColor' must be an R color name. All available colors",
               "can be listed with 'colors()'."),
         call. = FALSE)
  }
  if (!is.character(interceptType) |
      length(interceptType) != 1 |
      !interceptType %in% c("blank", "solid", "dashed", "dotted", "dotdash",
                            "longdash", "twodash")) {
    stop(paste("'interceptType' must be either 'blank', 'solid', 'dashed'",
               "(default), 'dotted', 'dotdash', 'longdash' or 'twodash'.",
               "For additional details see ?plotCorrelation"),
         call. = FALSE)
  }
  if (!is.numeric(borderWidth) |
      length(borderWidth) != 1 |
      borderWidth < 0) {
    stop("'borderWidth' must be a non-neagtive number! (default is 1)",
         call. = FALSE)
  }
  if (length(colorScale) != 3 |
      any(!colorScale %in% grDevices::colors())) {
    stop(paste("'colorScale' must be a vector with R color names for",
               "downregulated features, non significant features, and",
               "upregulated features. The default value is",
               "c('blue', 'grey', 'red'). All available colors",
               "can be listed with 'colors()'."), call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot.",
               "For additional details see ?plotVolcano"),
         call. = FALSE)
  }
  
  ## extract miRNA/gene differential expression and cutoffs
  if (assay == "microRNA") {
    featDE <- mirnaDE(mirnaObj, onlySignificant = FALSE)
    pCut <- mirnaObj@mirnaDE$pCutoff
    lCut <- mirnaObj@mirnaDE$logFC
  } else if (assay == "genes") {
    featDE <- geneDE(mirnaObj, onlySignificant = FALSE)
    pCut <- mirnaObj@geneDE$pCutoff
    lCut <- mirnaObj@geneDE$logFC
  }
  
  ## determine p-value cutoff based on adjusted p-value
  if (identical(featDE$P.Value, featDE$adj.P.Val)) {
    pCutoff <- pCut
  } else {
    pCutoff <- max(featDE$P.Value[featDE$adj.P.Val <= pCut])
  }
  
  ## determine Up and Downregulated features
  featDE$Change <- "Stable"
  featDE$Change[which(featDE$logFC > lCut &
                        featDE$adj.P.Val < pCut)] <- "Up"
  featDE$Change[which(featDE$logFC < -lCut &
                        featDE$adj.P.Val < pCut)] <- "Down"
  
  ## determine the labels to show
  if (is.character(labels)) {
    if (!all(labels %in% featDE$ID)) {
      stop("Not all specified labels are present in this assay!",
           call. = FALSE)
    } else {
      featDE$ID[which(!featDE$ID %in% labels)] <- ""
    }
  } else if (is.numeric(labels)) {
    fcFeat <- featDE[abs(featDE$logFC) > lCut, ]
    featDE$ID[which(!featDE$ID %in% fcFeat$ID[seq(labels)])] <- ""
  }
  
  ## produce a volcano plot
  pVol <- ggplot2::ggplot(data = featDE, 
                          ggplot2::aes(x = logFC, 
                                       y = -log10(P.Value), 
                                       colour = Change,
                                       label = ID)) +
    ggplot2::geom_point(alpha = pointAlpha, size = pointSize) +
    ggplot2::scale_color_manual(values = colorScale) +
    ggplot2::geom_vline(xintercept = c(-lCut, lCut), lty = interceptType,
                        col = interceptColor, lwd = interceptWidth) +
    ggplot2::geom_hline(yintercept = -log10(pCutoff), lty = interceptType,
                        col = interceptColor, lwd = interceptWidth) +
    ggplot2::labs(x = "log2(fold change)",
                  y = "-log10 (p-value)")  +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                   legend.position = "right", 
                   legend.title = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(
                     colour = "black", linewidth = borderWidth))
  
  ## add desired labels through ggrepel
  if (!is.null(labels)) {
    if (boxedLabel == TRUE) {
      pVol <- pVol +
        ggrepel::geom_label_repel(show.legend = FALSE)
    } else {
      pVol <- pVol +
        ggrepel::geom_text_repel(show.legend = FALSE)
    }
  }
  
  ## add title if desired by the user
  if (!is.null(title)) {
    pVol <- pVol +
      ggplot2::ggtitle(title)
  }
  
  ## return plot
  return(pVol)
  
}





#' Generate multidimensional scaling (MDS) plots to explore miRNA/gene
#' expression distances
#' 
#' This function performs multidimensional scaling in order to produce a simple
#' scatterplot that shows miRNA/gene expression variations among samples. In
#' particular, starting from a [`MirnaExperiment`][MirnaExperiment-class]
#' object, this functions allows to visualize both miRNA and gene expression in
#' the multidimensional space. Moreover, it is possible to color samples on the
#' basis of specific variables, and this is extremely useful to assess
#' miRNA/gene expression variations between distinct biological groups.
#' 
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param assay The results to display. It must be either 'microRNA', to plot
#' miRNA expression, or 'genes', to produce MDS plot for genes
#' @param condition It must be the column name of a variable specified in the
#' metadata (colData) of a [`MirnaExperiment`][MirnaExperiment-class] object;
#' or, alternatively, it must be a character/factor object that specifies
#' group memberships (eg. c("healthy, "healthy", "disease", "disease"))
#' @param dimensions It is a numeric vector of length 2 that indicates the two
#' dimensions to represent on the plot. Default is `c(1, 2)` to plot the two
#' dimensions that account for the highest portion of variability
#' @param labels Logical, whether to display labels or not. Default is FALSE
#' @param boxedLabel Logical, whether to show labels inside a rectangular shape
#' (default) or just as text elements
#' @param pointSize The size of points in the MDS plot (default is 3)
#' @param pointAlpha The transparency of points in the MDS plot (default
#' is 1)
#' @param borderWidth The width of plot borders (default is 1)
#' @param colorScale It must be a named character vector where values
#' correspond to R colors, while names coincide with the groups specified in
#' the `condition` parameter (eg. c("healthy" = "green", "disease" = "red")).
#' Default is NULL, in order to use the default color scale
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#' @param ... Other parameters that can be passed to [limma::plotMDS()] function
#'
#' @returns
#' An object of class `ggplot` containing the plot.
#' 
#' @note
#' To perform multidimensional scaling, this function internally uses
#' [limma::plotMDS()] function provided by `limma` package.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # produce MDS plot for miRNA expression with labels
#' plotDimensions(obj, "microRNA", condition = "disease", labels = TRUE)
#' 
#' # produce MDS plot for genes without condition color
#' plotDimensions(obj, "genes")
#' 
#' @references
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma
#' powers differential expression analyses for RNA-sequencing and microarray
#' studies.” Nucleic Acids Research, 43(7), e47. \url{doi:10.1093/nar/gkv007}.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
plotDimensions <- function(mirnaObj,
                           assay,
                           condition = NULL,
                           dimensions = c(1, 2),
                           labels = FALSE,
                           boxedLabel = TRUE,
                           pointSize = 3,
                           pointAlpha = 1,
                           borderWidth = 1,
                           colorScale = NULL,
                           title = NULL,
                           ...) {
  
  ## input checks
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (!assay %in% c("microRNA", "genes")) {
    stop("'assay' must be either 'microRNA' or 'genes'!", call. = FALSE)
  }
  if (length(condition) == 1) {
    if (!is.character(condition) |
        !(condition %in% colnames(MultiAssayExperiment::colData(mirnaObj)) &
          !condition %in% c("primary", "mirnaCol", "geneCol"))) {
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
  if (!is.numeric(dimensions) |
      length(dimensions) != 2 |
      any(dimensions%%1!=0) |
      any(dimensions <= 0)) {
    stop(paste("'dimensions' must be a numeric vector of length 2, that",
               "specifies the dimensions to be represented in the MDS plot!",
               "(e.g. 'c(1, 2)')"),
         call. = FALSE)
  }
  if (dimensions[1] > dimensions[2]) {
    dimensions <- c(dimensions[2], dimensions[1])
  }
  if (!is.logical(labels) |
      length(labels) != 1) {
    stop("'labels' must be logical (TRUE/FALSE)! See ?plotDimensions",
         call. = FALSE)
  }
  if (!is.logical(boxedLabel) |
      length(boxedLabel) != 1) {
    stop("'boxedLabel' must be logical (TRUE/FALSE)! See ?plotDimensions",
         call. = FALSE)
  }
  if (!is.numeric(pointSize) |
      length(pointSize) != 1 |
      pointSize < 0) {
    stop("'pointSize' must be a non-neagtive number! (default is 3)",
         call. = FALSE)
  }
  if (!is.numeric(pointAlpha) |
      length(pointAlpha) != 1 |
      pointAlpha < 0 |
      pointAlpha > 1) {
    stop("'pointAlpha' must be a number between 0 and 1! (default is 1)",
         call. = FALSE)
  }
  if (!is.numeric(borderWidth) |
      length(borderWidth) != 1 |
      borderWidth < 0) {
    stop("'borderWidth' must be a non-neagtive number! (default is 1)",
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
                 "conditions. For additional details see ?plotDimensions"),
           call. = FALSE)
    }
  } else if (length(condition) != 1 & !is.null(colorScale)) {
    if ((!is.null(colorScale) & !is.character(colorScale)) |
        any(!colorScale %in% grDevices::colors()) |
        !identical(sort(names(colorScale)),
                   as.character(sort(unique(condition))))) {
      stop(paste("'colorScale' must be a named character vector where values",
                 "consist of R colors, whereas names coincide to the different",
                 "conditions. For additional details see ?plotDimensions"),
           call. = FALSE)
    }
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot.",
               "For additional details see ?plotDimensions"),
         call. = FALSE)
  }
  
  ## define assay column in sample map
  if (assay == "microRNA") {
    featCol <- "mirnaCol"
  } else if (assay == "genes") {
    featCol <- "geneCol"
  }
  
  ## extract expression values
  featExpr <- mirnaObj[[assay]]
  
  ## identify sample metadata
  samplesMetadata <- MultiAssayExperiment::colData(mirnaObj)
  meta <- samplesMetadata[!is.na(samplesMetadata[, featCol]), ]
  
  ## define condition vector
  if (is.character(condition) & length(condition) == 1) {
    cond <- meta[, condition]
  } else if (is.factor(condition)) {
    cond <- as.character(condition)
  } else if (is.character(condition) & length(condition) > 1) {
    cond <- condition
  } else {
    cond = rep(NA, nrow(meta))
  }
  
  ## define cpm values if expression values are count-based
  oldCounts <- MultiAssayExperiment::metadata(mirnaObj)[["oldCounts"]]
  if (is.null(oldCounts[[assay]])) {
    featExpr <- edgeR::cpm(featExpr, log = TRUE)
  }
  
  ## perform multidimensional scaling through limma
  mds <- limma::plotMDS(featExpr, dim.plot = dimensions, plot = FALSE, ...)
  
  ## extract variance explained by the selected dimensions
  var.exp <- mds$var.explained[dimensions]
  
  ## create a data.frame object with MDS coordinates and covariates
  mds <- data.frame(x = mds$x,
                    y = mds$y,
                    primary = meta$primary,
                    condition = cond)
  
  ## create MDS plot based on biological condition or not
  if (!is.null(condition)) {
    mdsPlot <- ggplot2::ggplot(mds, ggplot2::aes(x = x,
                                                 y = y,
                                                 color = condition,
                                                 label = primary))
  } else {
    mdsPlot <- ggplot2::ggplot(mds, ggplot2::aes(x = x,
                                                 y = y,
                                                 label = primary))
  }
  
  ## add points to the MDS scatterplot
  mdsPlot <- mdsPlot +
    ggplot2::geom_point(alpha = pointAlpha, size = pointSize)
  
  ## add the desired color scale
  if (!is.null(colorScale)) {
    mdsPlot <- mdsPlot +
      ggplot2::scale_color_manual(values = colorScale)
  }
  
  ## add plot theme
  mdsPlot <- mdsPlot +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                   legend.position = "right",
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(
                     colour = "black", linewidth = borderWidth))
  
  ## set axis labels
  mdsPlot <- mdsPlot +
    ggplot2::xlab(paste("Leading logFC dim ", dimensions[1], " (",
                        round(var.exp[1]*100), "%)", sep = "")) +
    ggplot2::ylab(paste("Leading logFC dim ", dimensions[2], " (",
                        round(var.exp[2]*100), "%)", sep = ""))
  
  ## add labels
  if (labels == TRUE) {
    if (boxedLabel == TRUE) {
      mdsPlot <- mdsPlot +
        ggrepel::geom_label_repel(show.legend = FALSE)
    } else {
      mdsPlot <- mdsPlot +
        ggrepel::geom_text_repel(show.legend = FALSE)
    }
  }
  
  ## add title
  if (!is.null(title)) {
    mdsPlot <- mdsPlot +
      ggplot2::ggtitle(title)
  }
  
  ## return plot
  return(mdsPlot)
  
}





#' Display integrated miRNA-mRNA augmented pathways in a dotplot
#'
#' This function produces a dotplot that depicts the results of a
#' topologically-aware integrative pathway analysis (TAIPA) carried out through
#' the [topologicalAnalysis()] function.
#' 
#' When producing the dotplot with default values, the significant pathways
#' are ordered on the x-axis on the basis of their pathway score computed by
#' [topologicalAnalysis()]. The higher is this score, and the more affected a
#' pathway result between biological condditions. Moroever, the size of each
#' dot is equal to the ratio between the number of nodes for which a
#' measurement is available, and the total number of nodes. Finally, the color
#' scale of dots is relative to the adjusted p-values of each pathway. However,
#' the user can change the behavior of this function, by changing `ordBy`,
#' `sizeBy`, and `colBy` parameters.
#'
#' @param object An object of class
#' [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class]
#' @param showTerms It is the number of top significant pathways to show or,
#' alternatively, a character vector indicating the pathways to plot.
#' Default is `10`
#' @param ordBy The parameter used to set the x-axis scale. It must be one of
#' `ratio`, `padj`, `pval` and `score` (default). See the *details* section
#' for further information
#' @param sizeBy The parameter used to set the size scale. It must be one of
#' `ratio` (default), `padj`, `pval` and `score`. See the *details* section
#' for further information
#' @param colBy The parameter used to set the color scale. It must be one of
#' `ratio`, `padj` (default), `pval` and `score`. See the *details* section
#' for further information
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' A `ggplot` graph with a dotplot of integrated pathways.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform integration analysis with default settings
#' obj <- mirnaIntegration(obj)
#'
#' # perform the integrative pathway analysis with default settings
#' #ipa <- topologicalAnalysis(obj)
#' 
#' # access the results of pathway analysis
#' #integratedPathways(ipa)
#' 
#' # create a dotplot of integrated pathways
#' #integrationDotplot(ipa)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
integrationDotplot <- function(object,
                               showTerms = 10,
                               ordBy = "score",
                               sizeBy = "ratio",
                               colBy = "padj",
                               title = NULL) {
  
  ## check inputs
  if (!is(object, "IntegrativePathwayAnalysis")) {
    stop("'object' should be of class IntegrativePathwayAnalysis!",
         call. = FALSE)
  }
  if (nrow(integratedPathways(object)) < 1) {
    stop("'object' object does not contain any significant results!",
         call. = FALSE)
  }
  if (!is.character(showTerms) &
      !(is.numeric(showTerms) & length(showTerms) == 1)) {
    stop(paste("'showTerms' must be the number of top significant pathways",
               "to plot or, alternatively, a character vector containing",
               "the pathways to be shown."),
         call. = FALSE)
  }
  if (is.character(showTerms) &
      !all(showTerms %in% integratedPathways(object)$pathway)) {
    stop("The pathways provided are not present in 'object' pathway column!",
         call. = FALSE)
  }
  if (!is.character(ordBy) |
      length(ordBy) != 1 |
      !ordBy %in% c("ratio", "padj", "pval", "score")) {
    stop(paste("'ordBy' must be one of: 'ratio', 'padj',",
               "'pval', 'score' (default)"),
         call. = FALSE)
  }
  if (!is.character(sizeBy) |
      length(sizeBy) != 1 |
      !sizeBy %in% c("ratio", "padj", "pval", "score")) {
    stop(paste("'sizeBy' must be one of: 'ratio' (default), 'padj',",
               "'pval', 'score'"),
         call. = FALSE)
  }
  if (!is.character(colBy) |
      length(colBy) != 1 |
      !colBy %in% c("ratio", "padj", "pval", "score")) {
    stop(paste("'colBy' must be one of: 'ratio', 'padj'",
               "(default), 'pval', 'score'"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot.",
               "For additional details see ?integrationDotplot"),
         call. = FALSE)
  }
  
  ## extract results from object object
  res <- integratedPathways(object)
  
  ## CHANGE FROM HERE
  
  ## compute gene ratio and direction for different analyses
  if (objectMethod(object) == "Gene-Set object Analysis (GSEA)") {
    ov <- as.numeric(lapply(res$leadingEdge, length))
    res$overlap <- ov
    res$ratio <- ov/res$size
    res$direction <- "Up"
    res$direction[which(res$NES < 0)] <- "Down"
  } else {
    res$ratio <- res$overlap/res$size
  }
  
  ## order results based on padj
  res <- res[order(res$padj), ]
  
  ## select terms to be shown in the dotplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$pathway %in% showTerms), ]
  }
  
  ## set an x-axis label for 'ratio'
  if (ordBy == "ratio") {
    ordLabel <- "Gene-Set Overlap"
  } else {
    ordLabel <- ordBy
  }
  
  ## reformat results based on specified criterion
  if (ordBy != "padj" & ordBy != "pval") {
    res <- res[order(res[, ordBy], decreasing = TRUE), ]
  } else {
    res[, ordBy] <- -log10(res[, ordBy])
    res <- res[order(res[, ordBy], decreasing = TRUE), ]
    ordLabel <- paste("-log10(", ordBy, ")", sep = "")
  }
  
  ## create a dotplot
  dotRes <- ggplot2::ggplot(res,
                            ggplot2::aes(x = !!ggplot2::sym(ordBy),
                                         y = stats::reorder(.data$pathway,
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
  
  ## divide by object direction
  if (splitDir == TRUE &
      objectMethod(object) != "Over-Representation Analysis (ORA)") {
    dotRes <- dotRes +
      ggplot2::facet_grid(~ direction) +
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

