#' Visualize the relationships between miRNAs and genes in a biological pathway
#' 
#' This function can be used to plot augmented pathways created by the
#' [topologicalAnalysis()] function. In particular, given a valid object of
#' class [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class],
#' this function allows to produce a network graph for a specified biological
#' pathway, alongside with expression fold changes. In this way, augmented
#' pathways made of both miRNAs and genes can be visually explored to better
#' investigate the consequences of miRNA/gene dysregulations.
#'
#' @details
#' The network created by this function is highly flexible, allowing to tweak
#' different parameters that can influence the resulting graph, including
#' node highlighting, layout algorithms, colors, and legends.
#' 
#' ## Nodes included in the plot
#' 
#' For huge messy networks, the user can specify the nodes to include in the
#' plot through the `subgraph` parameter, in order to represent only the
#' features that he wants to display. Alternatively, this parameter can be set
#' to NULL (default), to plot all nodes of that biological pathway.
#' 
#' ## Highlight nodes and edges
#' 
#' One interesting feature offered by this function consists in highlighting
#' specific nodes and edges within a network. This results particularly useful
#' when we want to put in evidence affected routes in a biological pathway.
#' To highlight nodes, you must provide the `highlightNodes` parameter with a
#' `character` vector that lists all the desired nodes. As a result, the
#' borders of highlighted nodes will be colored according to `highlightCol`
#' parameter (default is 'gold'), and will have a width specified by
#' `highlightWidth` (default is 2). Notably, this function automatically
#' highlights in the same way the edges connecting the selected nodes.
#' 
#' ## Layout algorithms
#' 
#' Furthermore, this function allows to use different methods to lay out nodes
#' in the network by setting the `algorithm` parameter. In this regard, several
#' algorithms in `Rgraphviz` package can be used, namely:
#' 
#' * `dot` (default), which is an algorithm attributed to Sugiyama et al. and
#' described by Gansner et al., that creates a ranked layout that is
#' particularly suited to display hierarchies and complex pathways;
#' * `circo`, which uses a recursive radial algorithm resulting in a circular
#' layout;
#' * `fdp`, which adopts a force-directed approach similar to that of
#' Fruchterman and Reingold;
#' * `neato`, which relies on a spring model where an iterative solver finds
#' low energy configurations;
#' * `osage`, which is a layout engine that recursively draws cluster subgraphs;
#' * `twopi`, which places a node in the center of the network, and then
#' arranges the remaining nodes in a series of concentric circles around the
#' center.
#' 
#' For additional information on these algorithms, refer to
#' [Rgraphviz::GraphvizLayouts].
#' 
#' ## Customization
#' 
#' To customize the look of the resulting plot, this function allows to change
#' different graphical parameters, including:
#' 
#' * the color scale for log2 fold changes, that can be set with `lfcScale`;
#' * the font size of nodes, which can be changed through `fontsize`;
#' * the border color for nodes, which can be edited with `nodeBorderCol`;
#' * the text color of nodes, which can be changed through `nodeTextCol`;
#' * the color used for edges, set by `edgeCol`;
#' * the width of edges, customizable with `edgeWidth`.
#' 
#' Additionally, this function allows to include handy legends that are useful
#' for interpreting the biological consequences of network alterations. In
#' particular:
#' 
#' * a color bar legend displaying the log2 fold changes corresponding to each
#' fill color can be included with `legendColorbar = TRUE` (default); and
#' * a legend that links the appearance of edges and arrow heads to the type of
#' biological interaction can be shown through `legendInteraction = TRUE`
#' (default).
#' 
#' Lastly, `title`, `titleCex` and `titleFace` parameters can be tweaked to
#' include a network title with the desired look.
#' 
#' @param object An object of class
#' [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class] containing
#' the results of a miRNA-mRNA pathway analysis
#' @param pathway The name of the biological pathway to show. The available
#' pathways for a given database can be seen through the [listPathways()]
#' function
#' @param algorithm The layout algorithm used to arrange nodes in the network.
#' It must be one of `dot` (default), `circo`, `fdp`, `neato`, `osage` or
#' `twopi`. For more information regarding these algorithms, please check the
#' *details* section
#' @param fontsize The font size of each node in the graph. Default is 14
#' @param lfcScale It must be a `character` vector of length 3 containing valid
#' R color names for creating a gradient of log2 fold changes. The first value
#' refers to downregulation, the middle one to stable expression, and the last
#' one to upregulation. Default value is `c('royalblue', 'white', 'red')`.
#' Available color formats include color names, such as 'blue' and 'red', and
#' hexadecimal colors specified as #RRGGBB
#' @param nodeBorderCol It must be an R color name that specifies the color
#' of node borders. Default is `black`. Available color formats include
#' color names, such as 'blue' and 'red', and hexadecimal colors specified
#' as #RRGGBB
#' @param nodeTextCol It must be an R color name that specifies the color of
#' miRNA/gene names. Default is `black`. Available color formats include
#' color names, such as 'blue' and 'red', and hexadecimal colors specified
#' as #RRGGBB
#' @param edgeCol It must be an R color name that specifies the color of
#' edges between nodes. Default is `darkgrey`. Available color formats include
#' color names, such as 'blue' and 'red', and hexadecimal colors specified
#' as #RRGGBB
#' @param edgeWidth The width of edges. Default is 1
#' @param subgraph An optional `character` vector containing the nodes that you
#' want to maintain in the final plot. All the other nodes will not be shown.
#' This is useful to display specific features of extremely messy graphs.
#' Default is NULL
#' @param highlightNodes A `character` vector containing the names of nodes
#' that you want to highlight. Default is NULL not to highlight any nodes.
#' See the *details* section for additional information
#' @param highlightCol It must be an R color name that specifies the color of
#' edges and borders for highlighted nodes. Default is `gold`. Available color
#' formats include color names, such as 'blue' and 'red', and hexadecimal colors
#' specified as #RRGGBB
#' @param highlightWidth The width of edges between highlighted nodes. Default
#' is 2
#' @param legendColorbar Logical, whether to add a legend with a color bar for
#' log2 fold changes. Default is TRUE
#' @param legendInteraction Logical, whether to add a legend that links edge
#' types to biological interactions. Default is TRUE
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#' @param titleCex The cex of the plot main title. Default is 2
#' @param titleFace An integer which specifies which font to use for title. 1
#' corresponds to plain text, 2 to bold face, 3 to italic, 4 to bold italic,
#' and 5 to symbol font. Default is 1
#'
#' @returns
#' A base R plot with the augmented pathway.
#'
#' @examples
#' # load example IntegrativePathwayAnalysis object
#' obj <- loadExamples("IntegrativePathwayAnalysis")
#' 
#' \donttest{
#' # explore a specific biological network
#' visualizeNetwork(obj, "Thyroid hormone synthesis")
#' }
#' 
#' @note
#' This function uses the `Rgraphviz` package to render the network object.
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
visualizeNetwork <- function(object,
                             pathway,
                             algorithm = "dot",
                             fontsize = 14,
                             lfcScale = c("royalblue", "white", "red"),
                             nodeBorderCol = "black",
                             nodeTextCol = "black",
                             edgeCol = "darkgrey",
                             edgeWidth = 1,
                             subgraph = NULL,
                             highlightNodes = NULL,
                             highlightCol = "gold",
                             highlightWidth = 2,
                             legendColorbar = TRUE,
                             legendInteraction = TRUE,
                             title = NULL,
                             titleCex = 2,
                             titleFace = 1) {
  
  ## input checks
  if (!is(object, "IntegrativePathwayAnalysis")) {
    stop(paste("'object' should be of class IntegrativePathwayAnalysis!",
               "See ?IntegrativePathwayAnalysis-class"),
         call. = FALSE)
  }
  if (!is.character(pathway) |
      length(pathway) != 1 |
      !pathway %in% names(augmentedPathways(object))) {
    stop(paste("'pathway' must be the name of a valid pathway included in",
               "'names(augmentedPathways(object))'. For additional details,",
               "see ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!is.character(algorithm) |
      length(algorithm) != 1 |
      !algorithm %in% c("dot", "circo", "fdp", "neato", "osage", "twopi")) {
    stop(paste("'algorithm' must be either 'dot' (default), 'circo', 'fdp'",
               "'neato', 'osage' or 'twopi'.",
               "For additional details see ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!is.numeric(fontsize) |
      length(fontsize) != 1 |
      fontsize < 1) {
    stop("'fontsize' must be a number higher than 1! (default is 14)",
         call. = FALSE)
  }
  if (length(lfcScale) != 3 |
      any(areColors(lfcScale) == FALSE)) {
    stop(paste("'lfcScale' must be a vector with R color names for",
               "downregulated features, non significant features, and",
               "upregulated features. The default value is",
               "c('royalblue', 'white', 'red')."), call. = FALSE)
  }
  if (!is.character(nodeBorderCol) |
      length(nodeBorderCol) != 1 |
      areColors(nodeBorderCol) == FALSE) {
    stop(paste("'nodeBorderCol' must be an R color name."),
         call. = FALSE)
  }
  if (!is.character(nodeTextCol) |
      length(nodeTextCol) != 1 |
      areColors(nodeTextCol) == FALSE) {
    stop(paste("'nodeTextCol' must be an R color name."),
         call. = FALSE)
  }
  if (!is.character(edgeCol) |
      length(edgeCol) != 1 |
      areColors(edgeCol) == FALSE) {
    stop(paste("'edgeCol' must be an R color name."),
         call. = FALSE)
  }
  if (!is.numeric(edgeWidth) |
      length(edgeWidth) != 1 |
      edgeWidth < 0) {
    stop("'edgeWidth' must be a non-negative number! (default is 1)",
         call. = FALSE)
  }
  if (!is.null(subgraph) &
      !is.character(subgraph)) {
    stop(paste("'subgraph' must contain the names of nodes to",
               "maintain. For additional details see ?visualizeNetwork."),
         call. = FALSE)
  }
  if (!is.null(highlightNodes) &
      !is.character(highlightNodes)) {
    stop(paste("'highlightNodes' must contain the names of nodes to",
               "highlight. For additional details see ?visualizeNetwork."),
         call. = FALSE)
  }
  if (!is.character(highlightCol) |
      length(highlightCol) != 1 |
      areColors(highlightCol) == FALSE) {
    stop(paste("'highlightCol' must be an R color name."),
         call. = FALSE)
  }
  if (!is.numeric(highlightWidth) |
      length(highlightWidth) != 1 |
      highlightWidth < 0) {
    stop("'highlightWidth' must be a non-negative number! (default is 2)",
         call. = FALSE)
  }
  if (!is.logical(legendColorbar) |
      length(legendColorbar) != 1) {
    stop(paste("'legendColorbar' must be logical (TRUE/FALSE)!",
               "See ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!is.logical(legendInteraction) |
      length(legendInteraction) != 1) {
    stop(paste("'legendInteraction' must be logical (TRUE/FALSE)!",
               "See ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the network.",
               "For additional details see ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!is.numeric(titleCex) |
      length(titleCex) != 1 |
      titleCex <= 0) {
    stop("'titleCex' must be a positive number higher than 0! (default is 2)",
         call. = FALSE)
  }
  if (!is.numeric(titleFace) |
      length(titleFace) != 1 |
      titleFace < 1 |
      titleFace > 5 |
      titleFace%%1 != 0) {
    stop("'titleFace' must be an integer between 1 and 5! (default is 1)",
         call. = FALSE)
  }
  
  ## extract augmented pathways
  pList <- augmentedPathways(object)
  
  ## select the pathways of interest
  p <- pList[[pathway]]
  
  ## return NULL for invalid pathways
  if (is.null(p)) {
    stop("The selected pathway is NULL!", call. = FALSE)
  }
  
  ## retain just the nodes specified by the user
  if (!is.null(subgraph)) {
    
    ## check that supplied nodes are present
    if (any(!subgraph %in% graph::nodes(p))) {
      warning(paste(paste(subgraph[!subgraph %in% graph::nodes(p)],
                          collapse = ", "), "not belonging to this network!",
                    "Ignoring these nodes..."),
              call. = FALSE)
      subgraph <- subgraph[subgraph %in% graph::nodes(p)]
    }
    
    ## maintein only desired nodes
    p <- graph::subGraph(subgraph, p)
    
  }
  
  ## extract expression changes
  featDE <- object@expression
  
  ## extract weights and edges
  e <- vapply(Rgraphviz::edgeData(p), function(x) {
    x$weight
  }, FUN.VALUE = numeric(1))
  
  ## extract nodes
  nodes <- graph::nodes(p)
  
  ## retain expression changes just for nodes
  exprNodes <- featDE[nodes, ]
  
  ## change edge names according to Rgraphviz
  names(e) <- gsub("|", "~", names(e), fixed = TRUE)
  
  ## use different shapes for genes and miRNAs
  nS <- rep("ellipse", length(nodes))
  nS[exprNodes$type == "miRNA"] <- "box"
  names(nS) <- nodes
  
  ## use larger node shapes for miRNAs
  fixedsize <- rep(TRUE, length(nodes))
  fixedsize[exprNodes$type == "miRNA"] <- FALSE
  names(fixedsize) <- nodes
  
  ## set color scale for genes and miRNA fold changes
  filCol <- setColorScale(exprNodes$logFC,
                          cols = lfcScale,
                          numColors = 300)
  names(filCol) <- nodes
  
  ## change arrowheads on the basis of interaction type
  arrHead <- rep("open", length(e))
  arrHead[e == 0] <- "none"
  arrHead[e == -1] <- "tee"
  names(arrHead) <- names(e)
  
  ## set linetype dashed for binding
  lty <- rep("solid", length(e))
  lty[e == 0] <- "dashed"
  names(lty) <- names(e)
  
  ## arrange nodes in the network with respect to node shape
  g <- Rgraphviz::layoutGraph(p, layoutType = algorithm,
                              nodeAttr = list(shape = nS,
                                              fixedsize = fixedsize,
                                              fontsize = fontsize))
  
  ## define node parameters
  nAttr <- list(fill = filCol,
                fontsize = fontsize,
                col = nodeBorderCol,
                textCol = nodeTextCol)
  
  ## define edge parameters
  eAttr <- list(col = edgeCol,
                arrowhead = arrHead,
                lty = lty,
                lwd = edgeWidth)
  
  ## setting parameters for nodes and adges
  graph::nodeRenderInfo(g) <- nAttr
  graph::edgeRenderInfo(g) <- eAttr
  
  ## highlight desired nodes
  if (!is.null(highlightNodes)) {
    
    ## check that supplied nodes are present
    if (any(!highlightNodes %in% nodes)) {
      warning(paste(paste(highlightNodes[!highlightNodes %in% nodes],
                          collapse = ", "), "not belonging to this network!",
                    "Highlighting is ignored..."),
              call. = FALSE)
      highlightNodes <- highlightNodes[highlightNodes %in% nodes]
    }
    
    ## set the properties of selected nodes and edges
    highNodeCol <- rep(highlightCol, length(highlightNodes))
    highNodeWidth <- rep(highlightWidth, length(highlightNodes))
    names(highNodeCol) <- highlightNodes
    names(highNodeWidth) <- highlightNodes
    highEdges <- as.character(outer(highlightNodes, highlightNodes,
                                    paste, sep = "~"))
    highEdges <- highEdges[highEdges %in% graph::edgeNames(g)]
    highEdgeCol <- rep(highlightCol, length(highEdges))
    highEdgeWidth <- rep(highlightWidth, length(highEdges))
    names(highEdgeCol) <- highEdges
    names(highEdgeWidth) <- highEdges
    if (length(highNodeCol) > 0) {
      graph::nodeRenderInfo(g) <- list(col = highNodeCol, lwd = highNodeWidth)
    }
    if (length(highEdgeCol) > 0) {
      graph::edgeRenderInfo(g) <- list(col = highEdgeCol, lwd = highEdgeWidth)
    }
  }
  
  ## store current layout settings
  currentPar <- par(no.readonly = TRUE)
  
  ## define plotting layout
  if (legendColorbar == TRUE &
      legendInteraction == TRUE) {
    layout(matrix(c(1, 1, 1, 0, 2, 3), ncol = 2),
           widths = c(4, 1),
           heights = c(0.5, 2, 2))
  } else if (legendColorbar == TRUE &
             legendInteraction == FALSE) {
    layout(matrix(c(1, 1, 1, 0, 2, 0), ncol = 2),
           widths = c(4, 1),
           heights = c(1.2, 2, 1.2))
  } else if (legendColorbar == FALSE &
             legendInteraction == TRUE) {
    layout(matrix(c(1, 1, 1, 0, 2, 0), ncol = 2),
           widths = c(4, 1),
           heights = c(1.4, 2, 1))
  } else {
    layout(1)
  }
  
  ## add extra space for title
  if (!is.null(title)) {
    par(oma = c(0, 0, 2, 0))
  }
  
  ## make sure to restore graphical parameters on exit
  on.exit({
    par(currentPar)
  })
  
  ## produce the network graph
  Rgraphviz::renderGraph(g)
  
  ## produce logFC color bar legend
  if (legendColorbar == TRUE) {
    par(mar = c(2, 4.2, 4, 4.2))
    legMat <- matrix(seq(min(exprNodes$logFC), max(exprNodes$logFC),
                         length.out = 300))
    image(x = 1,
          y = legMat,
          z = matrix(1:300, nrow=1),
          col = setColorScale(seq(min(exprNodes$logFC),
                                  max(exprNodes$logFC),
                                  len = 300),
                              cols = lfcScale,
                              numColors = 300),
          axes = FALSE,
          xlab = "",
          ylab = "",
          main = "LogFC",
          font.main = 1)
    axis(2, las = 1)
  }
  
  ## create a legend for interaction types
  if (legendInteraction == TRUE) {
    par(mar = c(4, 2.5, 2, 2.5))
    plot.new()
    arrows(0, 0.8, 1, 0.8, code = 0, lty = "dashed", col = edgeCol)
    arrows(0, 0.5, 1, 0.5, angle = 30, col = edgeCol, length = 0.07)
    arrows(0, 0.2, 1, 0.2, angle = 90, col = edgeCol, length = 0.07)
    text(x = 0, y = 0.9, labels = "Binding", adj = 0)
    text(x = 0, y = 0.6, labels = "Activation", adj = 0)
    text(x = 0, y = 0.3, labels = "Inhibition", adj = 0)
    title(main = "Interaction", font.main = 1)
  }
  
  ## add a title to the network plot
  if (!is.null(title)) {
    title(title, line = 0, outer = TRUE,
          cex.main = titleCex, font.main = titleFace)
  }
}





## helper function for creating and assigning color scales
setColorScale <- function(values, cols, numColors) {
  
  ## set colors for positive and negative logFCs
  if (all(values > 0)) {
    cols <- cols[c(2, 3)]
  } else if (all(values < 0)) {
    cols <- cols[c(1, 2)]
  }
  
  ## create a function that generates the color scale
  colorScale <- colorRampPalette(cols)
  
  ## generate the color scale with the desired number of colors
  colorVector <- colorScale(numColors)
  
  ## define breaks for values
  posRatio <- sum(values > 0) / length(values)
  breaks <- c(-Inf, rev(seq(0, min(values), len = numColors / 2)),
              seq(0, max(values), len = numColors / 2)[-1], Inf)
  
  ## assign color values based on another vector
  colorValues <- colorVector[cut(values, breaks = breaks)]
  
  ## set the names of color values equal to those of vector
  names(colorValues) <- names(values)
  
  ## return color values
  return(colorValues)
  
}





#' Create a dotplot for functional enrichment analysis
#'
#' This function produces a dotplot to show the results of functional
#' enrichment analyses carried out through over-representation analysis (ORA),
#' gene set enrichment analysis (GSEA), and competitive gene set test accounting
#' for inter-gene correlation (CAMERA). In particular, this function can take
#' as input enrichment results generated by the [enrichGenes()] function.
#' 
#' When producing a dotplot with this function, significant pathways are
#' ordered on the x-axis on the basis of the ratio between the number of
#' overlapping genes in that set, and the total number of genes in the set.
#' Moreover, the size of each dot is proportional to the number of overlapping
#' features. Finally, the color scale of dots is relative to the adjusted
#' p-values of each category.
#'
#' @param enrichment An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#' @param showTerms It is the number of terms to be shown, based on the
#' order determined by the parameter `showTermsParam`; or, alternatively, a
#' character vector indicating the terms to plot. Default is `10`
#' @param showTermsParam The order in which the top terms are selected as
#' specified by the `showTerms` parameter. It must be one of `ratio` (default),
#' `padj`, `pval` and `overlap`
#' @param splitDir Logical, if `TRUE` the resulting plot will be divided in
#' two columns on the basis of enrichment direction (Up and Down).
#' Default is `TRUE`. This only applies if enrichment method is GSEA or CAMERA
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' A `ggplot` graph with a dotplot of enrichment results.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # extract results
#' res <- enrichmentResults(obj)
#'
#' # plot results
#' enrichmentDotplot(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
enrichmentDotplot <- function(enrichment,
                              showTerms = 10,
                              showTermsParam = "ratio",
                              splitDir = TRUE,
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
  if (!is.character(showTermsParam) |
      length(showTermsParam) != 1 |
      !showTermsParam %in% c("ratio", "padj", "pval", "overlap")) {
    stop(paste("'showTermsParam' must be one of: 'ratio' (default), 'padj',",
               "'pval', 'overlap'"),
         call. = FALSE)
  }
  if (!is.logical(splitDir) |
      length(splitDir) != 1) {
    stop("'splitDir' must be logical (TRUE/FALSE)! See ?enrichmentDotplot",
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
  
  ## order results on the basis of showTermsParam
  if (showTermsParam != "padj" & showTermsParam != "pval") {
    res <- res[order(res[, showTermsParam], decreasing = TRUE), ]
  } else {
    res <- res[order(res[, showTermsParam], decreasing = FALSE), ]
  }
  
  ## select pathways to be shown in the dotplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$pathway %in% showTerms), ]
  }
  
  ## set the parameters for plotting
  ordBy <- "ratio"
  sizeBy <- "overlap"
  colBy <- "padj"
  
  ## set an x-axis label for 'ratio'
  ordLabel <- "Gene-Set Overlap"
  
  ## create a dotplot
  dotRes <- ggplot2::ggplot(res,
                            ggplot2::aes(x = !!ggplot2::sym(ordBy),
                                         y = reorder(.data$pathway,
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
#' as input enrichment results generated by the [enrichGenes()] function.
#' 
#' When producing a barplot with this function, significant pathways are
#' ordered on the x-axis on the basis of the ratio between the number of
#' overlapping genes in that set, and the total number of genes in the set.
#' Moreover, the color scale of dots is relative to the adjusted p-values of
#' each category.
#'
#' @param enrichment An object of class
#' [`FunctionalEnrichment`][FunctionalEnrichment-class] containing
#' enrichment results
#' @param showTerms It is the number of terms to be shown, based on the
#' order determined by the parameter `showTermsParam`; or, alternatively, a
#' character vector indicating the terms to plot. Default is `10`
#' @param showTermsParam The order in which the top terms are selected as
#' specified by the `showTerms` parameter. It must be one of `ratio` (default),
#' `padj`, `pval` and `overlap`
#' @param splitDir Logical, if `TRUE` the resulting plot will be divided in
#' two columns on the basis of enrichment direction (Up and Down).
#' Default is `TRUE`. This only applies if enrichment method is GSEA or CAMERA
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' A `ggplot` graph with a barplot of enrichment results.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # extract results
#' res <- enrichmentResults(obj)
#'
#' # plot results
#' enrichmentBarplot(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
enrichmentBarplot <- function(enrichment,
                              showTerms = 10,
                              showTermsParam = "ratio",
                              splitDir = TRUE,
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
  if (!is.character(showTermsParam) |
      length(showTermsParam) != 1 |
      !showTermsParam %in% c("ratio", "padj", "pval", "overlap")) {
    stop(paste("'showTermsParam' must be one of: 'ratio' (default), 'padj',",
               "'pval', 'overlap'"),
         call. = FALSE)
  }
  if (!is.logical(splitDir) |
      length(splitDir) != 1) {
    stop("'splitDir' must be logical (TRUE/FALSE)! See ?enrichmentBarplot",
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
  
  ## order results on the basis of showTermsParam
  if (showTermsParam != "padj" & showTermsParam != "pval") {
    res <- res[order(res[, showTermsParam], decreasing = TRUE), ]
  } else {
    res <- res[order(res[, showTermsParam], decreasing = FALSE), ]
  }
  
  ## select pathways to be shown in the dotplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$pathway %in% showTerms), ]
  }
  
  ## set the parameters for plotting
  ordBy <- "ratio"
  colBy <- "padj"
  
  ## set an x-axis label for 'ratio'
  ordLabel <- "Gene-Set Overlap"
  
  ## create a dotplot
  barRes <- ggplot2::ggplot(res,
                            ggplot2::aes(x = !!ggplot2::sym(ordBy),
                                         y = reorder(.data$pathway,
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
#' @param showTerms It is the number of terms to be shown, based on the
#' order determined by the parameter `showTermsParam`; or, alternatively, a
#' character vector indicating the terms to plot. Default is `10`
#' @param showTermsParam The order in which the top terms are selected as
#' specified by the `showTerms` parameter. It must be one of `ratio`, `padj`
#' (default), `pval` and `overlap`
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' An object of class `ggplot` containing the ridgeplot of GSEA results.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # extract results
#' res <- enrichmentResults(obj)
#'
#' # plot results
#' gseaRidgeplot(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#' 
#' @export
gseaRidgeplot <- function(enrichment,
                          showTerms = 10,
                          showTermsParam = "padj",
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
  if (!is.character(showTermsParam) |
      length(showTermsParam) != 1 |
      !showTermsParam %in% c("ratio", "padj", "pval", "overlap")) {
    stop(paste("'showTermsParam' must be one of: 'ratio', 'padj' (default),",
               "'pval', 'overlap'"),
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
  
  ## order results on the basis of showTermsParam
  if (showTermsParam != "padj" & showTermsParam != "pval") {
    res <- res[order(res[, showTermsParam], decreasing = TRUE), ]
  } else {
    res <- res[order(res[, showTermsParam], decreasing = FALSE), ]
  }
  
  ## select pathways to be shown in the dotplot
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
  
  ## set the parameters for plotting
  colBy <- "padj"
  
  ## create a ridgeplot
  ridPlot <- ggplot2::ggplot(ridgeDf,
                             ggplot2::aes(x = .data$val,
                                          y = reorder(.data$pathway,
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
#' the running score line. Default is `green`. Available color formats include
#' color names, such as 'blue' and 'red', and hexadecimal colors specified
#' as #RRGGBB
#' @param lineSize The line width of the running score line. Default is `1`
#' @param vlineColor It must be an R color name that specifies the color of
#' the vertical line indicating the enrichment score (ES). Default is `red`.
#' Available color formats include color names, such as 'blue' and 'red', and
#' hexadecimal colors specified as #RRGGBB
#' @param vlineSize The line width of the vertical line indicating the
#' enrichment score (ES). Default is `0.6`
#'
#' @returns
#' An object of class `ggplot` containing the GSEA plot.
#'
#' @examples
#' # load example FunctionalEnrichment object
#' obj <- loadExamples("FunctionalEnrichment")
#'
#' # extract results
#' res <- enrichmentResults(obj)
#'
#' # plot results
#' gseaPlot(obj, pathway = "Thyroid hormone synthesis")
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
      areColors(lineColor) == FALSE) {
    stop(paste("'lineColor' must be an R color name."),
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
      areColors(vlineColor) == FALSE) {
    stop(paste("'vlineColor' must be an R color name."),
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
  gsePlot <- ggplot2::ggplot(gseData,
                             ggplot2::aes(x = .data$position)) +
    ggplot2::geom_line(ggplot2::aes(y = runScore),
                       linewidth = lineSize,
                       color = lineColor) +
    ggplot2::geom_vline(xintercept = which(runScore == es),
                        color = vlineColor,
                        linewidth = vlineSize,
                        linetype = "dashed") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$ymin,
                                         ymax = .data$ymax)) +
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
    rankPlot <- ggplot2::ggplot(gseData,
                                ggplot2::aes(x = .data$position)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$position,
                                         y = .data$metric,
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





## helper function for the definition of a ggplot2 theme for MIRit
#' @import ggplot2
theme.MIRit <- function(base_size = 12,
                        base_family = "",
                        legend = "top",
                        borderWidth = 1,
                        allBorders = TRUE,
                        grid = FALSE) {
  
  ## set plot borders
  if(allBorders == TRUE){
    br <- element_rect(size = borderWidth,
                       fill = NA,
                       color = "black")
    al <- element_blank()
  } else {
    br <- element_blank()
    al <- element_line(linewidth = borderWidth/2, color = "black")
  }
  
  ## define ggplot2 theme
  th <- theme_bw(base_size = base_size, base_family = base_family) +
    theme(panel.border = br,
          axis.line = al,
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.key = element_blank(),
          legend.position = legend,
          plot.title = ggplot2::element_text(hjust = 0.5))
  
  ## set grid lines
  if (grid == FALSE) {
    th <- th + theme(panel.grid = element_blank())
  }
  
  ## return the theme
  return(th)
  
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
#' the SNP locus. Default is `lightblue`. Available color formats include
#' color names, such as 'blue' and 'red', and hexadecimal colors specified
#' as #RRGGBB
#' @param mirFill It must be an R color name that specifies the fill color of
#' the miRNA locus. Default is `orange`. Available color formats include
#' color names, such as 'blue' and 'red', and hexadecimal colors specified
#' as #RRGGBB
#' @param from The start position of the plotted genomic range. Default is
#' NULL to automatically determine an appropriate position
#' @param to The end position of the plotted genomic range. Default is
#' NULL to automatically determine an appropriate position
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
#' \dontrun{
#' # retrieve associated SNPs
#' association <- findMirnaSNPs(obj, disId)
#'
#' # visualize association as a trackplot
#' mirVariantPlot(variantId = varId, snpAssociation = association)
#' }
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
                           from = NULL,
                           to = NULL,
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
      !identical(colnames(snpAssociation),
                 c("variant", "gene", "miRNA.gene", "miRNA.precursor",
                   "chr", "position", "allele", "distance", "is_upstream",
                   "is_downstream", "mirnaStrand", "mirnaStartPosition",
                   "mirnaEndPosition"))) {
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
      areColors(snpFill) == FALSE) {
    stop(paste("'snpFill' must be an R color name. Default is: 'lightblue'"),
         call. = FALSE)
  }
  if (!is.character(mirFill) |
      length(mirFill) != 1 |
      areColors(mirFill) == FALSE) {
    stop(paste("'mirFill' must be an R color name. Default is: 'orange'"),
         call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot (e.g. 'SNPs overlap').",
               "For additional details see ?mirVariantPlot"),
         call. = FALSE)
  }
  if (!is.null(from)) {
    if (!is.numeric(from) |
        length(from) != 1 |
        from < 0) {
      stop("'from' must be a non-neagtive number!",
           call. = FALSE)
    }
  }
  if (!is.null(to)) {
    if (!is.numeric(to) |
        length(to) != 1 |
        to < 0) {
      stop("'to' must be a non-neagtive number!",
           call. = FALSE)
    }
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
  mirSeq <- snpAssociation[, c(5, 12, 13, 11, 2)]
  colnames(mirSeq) <- c("chr", "start", "end", "strand", "miRNA_gene")
  mirSeq$strand <- ifelse(mirSeq$strand == 1, "+", "-")
  
  ## create GRanges object containing miRNA positions
  mirSeq <- GenomicRanges::makeGRangesFromDataFrame(mirSeq,
                                                    keep.extra.columns = TRUE)
  
  ## set parameters
  chr <- paste("chr", snpAssociation$chr, sep = "")
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
                                    symbol = snpAssociation$miRNA.gene,
                                    fill = mirFill)
  
  ## create sequence track
  if (showSequence == TRUE) {
    sTrack <- Gviz::SequenceTrack(
      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
      chromosome = chr)
  }
  
  ## create genomic context track
  if (showContext == TRUE) {
    message("Retrieving genomic context from Ensembl...")
    mart <- biomaRt::useMart("ensembl")
    mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart)
    biomTrack <- Gviz::BiomartGeneRegionTrack(
      biomart = mart,
      genome = g,
      symbol = snpAssociation$miRNA.gene,
      name = "Genomic context")
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
  
  ## ignore title if missing
  if (is.null(title)) {
    title <- ""
  }
  
  ## create the trackplot object
  trackPlot <- Gviz::plotTracks(pList,
                                extend.left = lf,
                                extend.right = rf,
                                transcriptAnnotation = "symbol",
                                main = title,
                                collapseTranscripts = "meta",
                                from = from,
                                to = to,
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
#' @param condition It must be NULL (default) to plot expression based on the
#' group variable used for differential expression analysis. Alternatively, it
#' must be a character/factor object that specifies group memberships
#' (eg. c("healthy, "healthy", "disease", "disease"))
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
#' the regression line. Default is `red`. Available color formats include
#' color names, such as 'blue' and 'red', and hexadecimal colors specified
#' as #RRGGBB
#' @param lineType It specifies the line type used for the regression line. It
#' must be either 'blank', 'solid', 'dashed' (default), 'dotted', 'dotdash',
#' 'longdash' or 'twodash'
#' @param lineWidth The width of the fitted regression line (default is 0.8)
#' @param pointSize The size of points in the correlation plot (default is 3)
#' @param colorScale It must be a named character vector where values
#' correspond to R colors, while names coincide with the groups specified in
#' the `condition` parameter (eg. c("healthy" = "green", "disease" = "red")).
#' Default is NULL, in order to use the default color scale. Available color
#' formats include color names, such as 'blue' and 'red', and hexadecimal
#' colors specified as #RRGGBB
#' @param fontSize The base size for text elements within the plot.
#' Default is 12
#' @param fontFamily The base family for text elements within the plot
#' @param legend The position of the legend. Allowed values are `top`,
#' `bottom`, `right`, `left` and `none`. The default setting is `top` to show
#' a legend above the plot. If `none` is specified, the legend will not be
#' included in the graph.
#' @param borderWidth The width of plot borders (default is 1)
#' @param allBorders Logical, whetether to show all panel borders, or just the
#' bottom and left borders. Default is TRUE
#' @param grid Logical, whether to show grid lines or not. Default is TRUE
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
#' # plot correlation between miR-146b and PAX8 with monotonic regression curve
#' plotCorrelation(obj, "hsa-miR-146b-5p", "PAX8", condition = "disease")
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
plotCorrelation <- function(mirnaObj,
                            mirna,
                            gene,
                            condition = NULL,
                            showCoeff = TRUE,
                            regression = TRUE,
                            useRanks = FALSE,
                            lineCol = "red",
                            lineType = "dashed",
                            lineWidth = 0.8,
                            pointSize = 3,
                            colorScale = NULL,
                            fontSize = 12,
                            fontFamily = "",
                            legend = "top",
                            borderWidth = 1,
                            allBorders = TRUE,
                            grid = TRUE) {
  
  ## input checks
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (max(dim(integration(mirnaObj))) == 0) {
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
  if (is.null(condition)) {
    if (mirnaObj@mirnaDE$group != mirnaObj@geneDE$group) {
      stop(paste("For unpaired data, the 'group' variable used for",
                 "differential expression analysis must be the same for both",
                 "miRNAs and genes in order to use",
                 "this function with 'condition = NULL'. Instead, try to",
                 "supply 'condition' as a factor/character vector!"),
           call. = FALSE)
    }
    if (length(mirnaObj@mirnaDE$group) == 0 |
        length(mirnaObj@geneDE$group) == 0) {
      stop(paste("For objects where differential expression has been manually",
                 "added, 'condition' must be specified as a factor/character",
                 "vector!"),
           call. = FALSE)
    }
  }
  if (!is.null(condition)) {
    if (length(condition) == 1) {
      if (!is.character(condition) |
          !(condition %in% colnames(colData(mirnaObj)) &
            !condition %in% c("primary", "mirnaCol", "geneCol"))) {
        stop(paste("'condition' must be the column name of a variable present",
                   "in the metadata (colData) of a MirnaExperiment object; or,",
                   "alternatively, it must be a character/factor object that",
                   "specifies group memberships."),
             call. = FALSE)
      }
    } else {
      if ((!is.character(condition) & !is.factor(condition)) |
          length(condition) != nrow(colData(mirnaObj))) {
        stop(paste("'condition' must be the column name of a variable present",
                   "in the metadata (colData) of a MirnaExperiment object; or,",
                   "alternatively, it must be a character/factor object that",
                   "specifies group memberships."),
             call. = FALSE)
      }
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
      areColors(lineCol) == FALSE) {
    stop(paste("'lineCol' must be an R color name."),
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
  if (!is.numeric(fontSize) |
      length(fontSize) != 1 |
      fontSize < 0) {
    stop("'fontSize' must be a non-neagtive number! (default is 12)",
         call. = FALSE)
  }
  if (!is.character(fontFamily) |
      length(fontFamily) != 1) {
    stop("'fontFamily' must be a character of length 1",
         call. = FALSE)
  }
  if (!is.character(legend) |
      length(legend) != 1 |
      !legend %in% c("top", "bottom", "right", "left", "none")) {
    stop("'legend' must be one of 'top', 'bottom' 'right', 'left', and 'none'",
         call. = FALSE)
  }
  if (!is.numeric(borderWidth) |
      length(borderWidth) != 1 |
      borderWidth < 0) {
    stop("'borderWidth' must be a non-neagtive number! (default is 1)",
         call. = FALSE)
  }
  if (!is.logical(allBorders) |
      length(allBorders) != 1) {
    stop("'allBorders' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.logical(grid) |
      length(grid) != 1) {
    stop("'grid' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  
  ## get integration results
  intRes <- integration(mirnaObj)
  
  ## verify that correlation analysis has been performed
  if (grepl("correlation", mirnaObj@integration$method) == FALSE) {
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
  if (is.null(condition)) {
    depM <- mirnaDE(mirnaObj, param = TRUE)
    depG <- geneDE(mirnaObj, param = TRUE)
    cond <- as.character(colData(mirnaObj)[, depM$group])
    cond[is.na(cond)] <- as.character(
      colData(mirnaObj)[, depG$group])[is.na(cond)]
  } else if (is.character(condition) & length(condition) == 1) {
    cond <- colData(mirnaObj)[, condition]
  } else if (is.factor(condition)) {
    cond <- as.character(condition)
  } else {
    cond <- condition
  }
  names(cond) <- colData(mirnaObj)[, "primary"]
  
  ## check the validity of color scale
  if (!is.null(colorScale)) {
    if (!is.character(colorScale) |
        any(areColors(colorScale) == FALSE) |
        !identical(sort(names(colorScale)),
                   as.character(sort(unique(cond))))) {
      stop(paste("'colorScale' must be a named character vector where values",
                 "consist of R colors, whereas names coincide to the different",
                 "conditions. For additional details see ?plotCorrelation."),
           call. = FALSE)
    }
  }
  
  ## check if samples are paired, otherwise exclude unpaired samples
  sMap <- sampleMap(mirnaObj)
  mirnaSamples <- sMap$primary[sMap$assay == "microRNA"]
  geneSamples <- sMap$primary[sMap$assay == "genes"]
  
  if (!identical(mirnaSamples, geneSamples)) {
    
    ## determine common and uncommon samples
    common <- intersect(mirnaSamples, geneSamples)
    unpaired  <- setdiff(c(mirnaSamples, geneSamples), common)
    
    ## remove samples without measurments of both miRNAs and genes
    if (length(unpaired) > 0) {
      
      mirnaExpr <- mirnaExpr[, sMap$colname[sMap$assay == "microRNA" &
                                              sMap$primary %in% common]]
      geneExpr <- geneExpr[, sMap$colname[sMap$assay == "genes" &
                                            sMap$primary %in% common]]
      
    }
    
    ## order the columns of expression matrices in the same way
    mirnaMap = sMap[sMap$assay == "microRNA" &
                      sMap$primary %in% common, ]
    geneMap = sMap[sMap$assay == "genes" &
                     sMap$primary %in% common, ]
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
                             ggplot2::aes(.data$miRNA, .data$gene)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Condition),
                        size = pointSize)
  
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
      ySeq <- predict(monoFit, newdata = data.frame(x = xSeq))
      
      ## create a data.frame for the fitted curve
      fitDf <- data.frame("xFit" = xSeq, "yFit" = ySeq[, "x"])
      
      ## add the monotonic curve to the plot
      corPlot <- corPlot +
        ggplot2::geom_line(data = fitDf,
                           ggplot2::aes(x = .data$xFit,
                                        y = .data$yFit),
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
      ggpubr::stat_cor(ggplot2::aes(label = ggplot2::after_stat(.data$r.label)),
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
  
  ## rename plot labels
  corPlot <- corPlot +
    ggplot2::labs(x = xlab, y = ylab)
  
  ## apply MIRit ggplot2 theme
  corPlot <- corPlot +
    theme.MIRit(base_size = fontSize,
                base_family = fontFamily,
                legend = legend,
                borderWidth = borderWidth,
                allBorders = allBorders,
                grid = grid)
  
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
#' @param condition It must be NULL (default) to plot expression based on the
#' group variable used for differential expression analysis. Alternatively, it
#' must be a character/factor object that specifies group memberships
#' (eg. c("healthy, "healthy", "disease", "disease"))
#' @param graph The type of plot to produce. It must be one of `boxplot`
#' (default), `barplot`, `violinplot`
#' @param linear Logical, whether to plot expression levels in linear scale
#' or in log2 space. Default is TRUE in order to use the linear space
#' @param showSignificance Logical, whether to display statistical significance
#' or not. Default is TRUE
#' @param starSig Logical, whether to represent statistical significance through
#' stars. Default is TRUE, and the significance scale is: * for \eqn{p < 0.05},
#' ** for \eqn{p < 0.01}, *** for \eqn{p < 0.001}, and **** for
#' \eqn{p < 0.0001}. If `starSig` is set to FALSE, p-values or adjusted
#' p-values will be reported on the plot as numbers
#' @param pCol The statistics used to evaluate comparison significance. It must
#' be one of `P.Value`, to use unadjusted p-values, and `adj.P.Val` (default),
#' to use p-values corrected for multiple testing
#' @param sigLabelSize The size for the labels used to show statistical
#' significance. Default is 7, which is well suited for representing p-values
#' as significance stars. However, if `starSig` is set to FALSE, the user might
#' have to downsize this parameter
#' @param digits The number of digits to show when p-values are reported as
#' numbers (when `starSig` is FALSE). Default is 3
#' @param nameAsTitle Logical, if set to TRUE, the miRNA/gene name will be
#' added as plot title, and the x-axis and legend will be removed. Note that
#' this option is only considered if `features` contains just one miRNA/gene.
#' Default is FALSE
#' @param colorScale It must be a named character vector where values
#' correspond to R colors, while names coincide with the groups specified in
#' the `condition` parameter (eg. c("healthy" = "green", "disease" = "red")).
#' Default is NULL, in order to use the default color scale. Available color
#' formats include color names, such as 'blue' and 'red', and hexadecimal
#' colors specified as #RRGGBB
#' @param fontSize The base size for text elements within the plot.
#' Default is 12
#' @param fontFamily The base family for text elements within the plot
#' @param legend The position of the legend. Allowed values are `top`,
#' `bottom`, `right`, `left` and `none`. The default setting is `top` to show
#' a legend above the plot. If `none` is specified, the legend will not be
#' included in the graph.
#' @param borderWidth The width of plot borders (default is 1)
#' @param allBorders Logical, whetether to show all panel borders, or just the
#' bottom and left borders. Default is FALSE
#' @param grid Logical, whether to show grid lines or not. Default is FALSE
#'
#' @returns
#' An object of class `ggplot` containing the plot.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#'
#' # produce a boxplot for PAX8 and miR-34a-5p
#' plotDE(obj, features = c("hsa-miR-34a-5p", "PAX8"))
#' 
#' # produce a barplot for PAX8 and miR-34a-5p without significance
#' plotDE(obj, features = c("hsa-miR-34a-5p", "PAX8"),
#' graph = "barplot", showSignificance = FALSE)
#' 
#' # produce a violinplot for BCL2
#' plotDE(obj, features = "BCL2", graph = "violinplot")
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @importFrom ggpubr mean_sd
#' @export
plotDE <- function(mirnaObj,
                   features,
                   condition = NULL,
                   graph = "boxplot",
                   linear = TRUE,
                   showSignificance = TRUE,
                   starSig = TRUE,
                   pCol = "adj.P.Val",
                   sigLabelSize = 7,
                   digits = 3,
                   nameAsTitle = FALSE,
                   colorScale = NULL,
                   fontSize = 12,
                   fontFamily = "",
                   legend = "top",
                   borderWidth = 1,
                   allBorders = FALSE,
                   grid = FALSE) {
  
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
  if (is.null(condition)) {
    if (mirnaObj@mirnaDE$group != mirnaObj@geneDE$group) {
      stop(paste("For unpaired data, the 'group' variable used for",
                 "differential expression analysis must be the same for both",
                 "miRNAs and genes in order to use",
                 "this function with 'condition = NULL'. Instead, try to",
                 "supply 'condition' as a factor/character vector!"),
           call. = FALSE)
    }
    if (length(mirnaObj@mirnaDE$group) == 0 |
        length(mirnaObj@geneDE$group) == 0) {
      stop(paste("For objects where differential expression has been manually",
                 "added, 'condition' must be specified as a factor/character",
                 "vector!"),
           call. = FALSE)
    }
  }
  if (!is.null(condition)) {
    if (length(condition) == 1) {
      if (!is.character(condition) |
          !(condition %in% colnames(colData(mirnaObj)) &
            !condition %in% c("primary", "mirnaCol", "geneCol"))) {
        stop(paste("'condition' must be the column name of a variable present",
                   "in the metadata (colData) of a MirnaExperiment object; or,",
                   "alternatively, it must be a character/factor object that",
                   "specifies group memberships."),
             call. = FALSE)
      }
    } else {
      if ((!is.character(condition) & !is.factor(condition)) |
          length(condition) != nrow(colData(mirnaObj))) {
        stop(paste("'condition' must be the column name of a variable present",
                   "in the metadata (colData) of a MirnaExperiment object; or,",
                   "alternatively, it must be a character/factor object that",
                   "specifies group memberships."),
             call. = FALSE)
      }
    }
  }
  if (!is.character(graph) |
      length(graph) != 1 |
      !graph %in% c("boxplot", "barplot", "violinplot")) {
    stop(paste("'graph' must be either 'boxplot' (default), `barplot`",
               "or 'violinplot'. For additional details see ?plotDE"),
         call. = FALSE)
  }
  if (!is.logical(linear) |
      length(linear) != 1) {
    stop("'linear' must be logical (TRUE/FALSE)! See ?plotDE",
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
  if (!is.logical(nameAsTitle) |
      length(nameAsTitle) != 1) {
    stop("'nameAsTitle' must be logical (TRUE/FALSE)! See ?plotDE",
         call. = FALSE)
  }
  if (nameAsTitle == TRUE & length(features) > 1) {
    warning(paste("'nameAsTitle' can only be used for 1 miRNA/gene!",
                  "This option has benn set to FALSE. See ?plotDE"),
            call. = FALSE)
    nameAsTitle <- FALSE
  }
  if (!is.numeric(fontSize) |
      length(fontSize) != 1 |
      fontSize < 0) {
    stop("'fontSize' must be a non-neagtive number! (default is 12)",
         call. = FALSE)
  }
  if (!is.character(fontFamily) |
      length(fontFamily) != 1) {
    stop("'fontFamily' must be a character of length 1",
         call. = FALSE)
  }
  if (!is.character(legend) |
      length(legend) != 1 |
      !legend %in% c("top", "bottom", "right", "left", "none")) {
    stop("'legend' must be one of 'top', 'bottom' 'right', 'left', and 'none'",
         call. = FALSE)
  }
  if (!is.numeric(borderWidth) |
      length(borderWidth) != 1 |
      borderWidth < 0) {
    stop("'borderWidth' must be a non-neagtive number! (default is 1)",
         call. = FALSE)
  }
  if (!is.logical(allBorders) |
      length(allBorders) != 1) {
    stop("'allBorders' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.logical(grid) |
      length(grid) != 1) {
    stop("'grid' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  
  ## extract miRNA and gene expression values
  mirnaExpr <- mirnaObj[["microRNA"]]
  geneExpr <- mirnaObj[["genes"]]
  
  ## define condition vector
  if (is.null(condition)) {
    depM <- mirnaDE(mirnaObj, param = TRUE)
    depG <- geneDE(mirnaObj, param = TRUE)
    cond <- as.character(colData(mirnaObj)[, depM$group])
    cond[is.na(cond)] <- as.character(
      colData(mirnaObj)[, depG$group])[is.na(cond)]
    contrast <- strsplit(depG$contrast, "-")[[1]]
    lv1 <- contrast[1]
    lv2 <- contrast[2]
  } else if (is.factor(condition)) {
    cond <- as.character(condition)
    lv1 <- unique(cond)[1]
    lv2 <- unique(cond)[2]
  } else {
    cond <- condition
    lv1 <- unique(cond)[1]
    lv2 <- unique(cond)[2]
  }
  names(cond) <- colData(mirnaObj)[, "primary"]
  
  ## check the validity of color scale
  if (!is.null(colorScale)) {
    if (!is.character(colorScale) |
        any(areColors(colorScale) == FALSE) |
        !identical(sort(names(colorScale)),
                   as.character(sort(unique(cond))))) {
      stop(paste("'colorScale' must be a named character vector where values",
                 "consist of R colors, whereas names coincide to the different",
                 "conditions. For additional details see ?plotDE"),
           call. = FALSE)
    }
  }
  
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
    
    ## subset condition vector based on available samples
    subCond <- cond[names(featExpr)]
    
    ## return feature expression, name and condition
    newDf <- data.frame("Expression" = featExpr,
                        "Gene" = rep(gene, length(subCond)),
                        "Condition" = subCond)
    exprDf <- rbind(exprDf, newDf)
    
  }
  
  ## keep only the condition levels used for DE analysis
  exprDf <- exprDf[exprDf$Condition %in% c(lv1, lv2), ]
  
  ## transform data in linear space
  if (linear == TRUE) {
    exprDf$Expression <- 2^exprDf$Expression
  }
  
  ## produce the desired plot
  if (graph == "boxplot") {
    
    ## create a grouped boxplot
    dePlot <- ggpubr::ggboxplot(data = exprDf, x = "Gene",
                                y = "Expression", fill = "Condition")
    
  } else if (graph == "barplot") {
    
    ## create a grouped barplot
    dePlot <- ggpubr::ggbarplot(data = exprDf,
                                x = "Gene",
                                y = "Expression",
                                fill = "Condition",
                                position = ggplot2::position_dodge(0.8),
                                add = "mean_sd",
                                error.plot = "upper_errorbar")
    
  } else if (graph == "violinplot") {
    
    ## create a grouped violinplot
    dePlot <- ggpubr::ggviolin(data = exprDf,
                               x = "Gene",
                               y = "Expression",
                               fill = "Condition",
                               add = "boxplot")
    
  }
  
  ## change y-axis label for log2 expression data
  if (linear == FALSE) {
    dePlot <- dePlot + ggplot2::ylab(expression(paste(log[2], " expression")))
  }
  
  ## add significance levels
  if (showSignificance == TRUE) {
    
    ## load differential expression results
    statTest <- rbind(mirnaDE(mirnaObj, onlySignificant = FALSE),
                      geneDE(mirnaObj, onlySignificant = FALSE))
    
    ## restrict differential expression to the selected miRNAs/genes
    statTest <- statTest[features, ]
    
    ## add conditions to differential expression results
    statTest$group1 <- lv1
    statTest$group2 <- lv2
    
    ## add y position for p-value labels
    maxExpr <- vapply(features, function(g) {
      max(exprDf$Expression[exprDf$Gene == g])
    }, FUN.VALUE = numeric(1))
    maxExpr <- maxExpr + 0.1 * mean(maxExpr)
    statTest$y.position <- maxExpr
    
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
    
    ## add significance to the ggplot2 object
    dePlot <- dePlot +
      ggpubr::stat_pvalue_manual(data = statTest,
                                 label = pCol,
                                 x = "ID",
                                 size = sigLabelSize)
    
    ## expand y-limit
    dePlot <- dePlot +
      ggplot2::ylim(c(NA, max(statTest$y.position) +
                        0.05 * max(statTest$y.position)))
    
  }
  
  ## add colorScale to ggplot2 graph
  if (!is.null(colorScale)) {
    dePlot <- dePlot +
      ggplot2::scale_color_manual(values = colorScale) +
      ggplot2::scale_fill_manual(values = colorScale)
  }
  
  ## apply MIRit ggplot2 theme
  dePlot <- dePlot +
    theme.MIRit(base_size = fontSize,
                base_family = fontFamily,
                legend = legend,
                borderWidth = borderWidth,
                allBorders = allBorders,
                grid = grid)
  
  ## remove x-axis and set gene name as title
  if (nameAsTitle == TRUE) {
    dePlot <- dePlot + ggplot2::ggtitle(features) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
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
#' cutoff intercepts. Default is `black`. Available color formats include
#' color names, such as 'blue' and 'red', and hexadecimal colors specified
#' as #RRGGBB
#' @param interceptType It specifies the line type used for cutoff intercepts.
#' It must be either 'blank', 'solid', 'dashed' (default), 'dotted', 'dotdash',
#' 'longdash' or 'twodash'
#' @param colorScale It must be a character vector of length 3 containing valid
#' R color names for downregulated, non significant, and upregulated features,
#' respectively. Default value is `c('blue', 'grey', 'red')`. Available color
#' formats include color names, such as 'blue' and 'red', and hexadecimal
#' colors specified as #RRGGBB
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#' @param fontSize The base size for text elements within the plot.
#' Default is 12
#' @param fontFamily The base family for text elements within the plot
#' @param legend The position of the legend. Allowed values are `top`,
#' `bottom`, `right`, `left` and `none`. The default setting is `none` so that
#' the legend will not be included in the graph.
#' @param borderWidth The width of plot borders (default is 1)
#' @param allBorders Logical, whetether to show all panel borders, or just the
#' bottom and left borders. Default is TRUE
#' @param grid Logical, whether to show grid lines or not. Default is FALSE
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
                        colorScale = c("blue", "grey","red"),
                        title = NULL,
                        fontSize = 12,
                        fontFamily = "",
                        legend = "none",
                        borderWidth = 1,
                        allBorders = TRUE,
                        grid = FALSE) {
  
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
      areColors(interceptColor) == FALSE) {
    stop(paste("'interceptColor' must be an R color name."),
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
  if (length(colorScale) != 3 |
      any(areColors(colorScale) == FALSE)) {
    stop(paste("'colorScale' must be a vector with R color names for",
               "downregulated features, non significant features, and",
               "upregulated features. The default value is",
               "c('blue', 'grey', 'red')."), call. = FALSE)
  }
  if (!(is.character(title) | is.null(title)) |
      !length(title) %in% c(0, 1)) {
    stop(paste("'title' must be the title of the plot.",
               "For additional details see ?plotVolcano"),
         call. = FALSE)
  }
  if (!is.numeric(fontSize) |
      length(fontSize) != 1 |
      fontSize < 0) {
    stop("'fontSize' must be a non-neagtive number! (default is 12)",
         call. = FALSE)
  }
  if (!is.character(fontFamily) |
      length(fontFamily) != 1) {
    stop("'fontFamily' must be a character of length 1",
         call. = FALSE)
  }
  if (!is.character(legend) |
      length(legend) != 1 |
      !legend %in% c("top", "bottom", "right", "left", "none")) {
    stop("'legend' must be one of 'top', 'bottom' 'right', 'left', and 'none'",
         call. = FALSE)
  }
  if (!is.numeric(borderWidth) |
      length(borderWidth) != 1 |
      borderWidth < 0) {
    stop("'borderWidth' must be a non-neagtive number! (default is 1)",
         call. = FALSE)
  }
  if (!is.logical(allBorders) |
      length(allBorders) != 1) {
    stop("'allBorders' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.logical(grid) |
      length(grid) != 1) {
    stop("'grid' must be logical (TRUE/FALSE)!", call. = FALSE)
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
                          ggplot2::aes(x = .data$logFC, 
                                       y = -log10(.data$P.Value), 
                                       colour = .data$Change,
                                       label = .data$ID)) +
    ggplot2::geom_point(alpha = pointAlpha, size = pointSize) +
    ggplot2::scale_color_manual(values = colorScale) +
    ggplot2::geom_vline(xintercept = c(-lCut, lCut), lty = interceptType,
                        col = interceptColor, lwd = interceptWidth) +
    ggplot2::geom_hline(yintercept = -log10(pCutoff), lty = interceptType,
                        col = interceptColor, lwd = interceptWidth) +
    ggplot2::labs(x = "log2(fold change)",
                  y = "-log10 (p-value)")
  
  ## apply MIRit ggplot2 theme
  pVol <- pVol +
    theme.MIRit(base_size = fontSize,
                base_family = fontFamily,
                legend = legend,
                borderWidth = borderWidth,
                allBorders = allBorders,
                grid = grid)
  
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
#' @param colorScale It must be a named character vector where values
#' correspond to R colors, while names coincide with the groups specified in
#' the `condition` parameter (eg. c("healthy" = "green", "disease" = "red")).
#' Default is NULL, in order to use the default color scale. Available color
#' formats include color names, such as 'blue' and 'red', and hexadecimal
#' colors specified as #RRGGBB
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#' @param fontSize The base size for text elements within the plot.
#' Default is 12
#' @param fontFamily The base family for text elements within the plot
#' @param legend The position of the legend. Allowed values are `top`,
#' `bottom`, `right`, `left` and `none`. The default setting is `top` to show
#' a legend above the plot. If `none` is specified, the legend will not be
#' included in the graph.
#' @param borderWidth The width of plot borders (default is 1)
#' @param allBorders Logical, whetether to show all panel borders, or just the
#' bottom and left borders. Default is TRUE
#' @param grid Logical, whether to show grid lines or not. Default is FALSE
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
                           colorScale = NULL,
                           title = NULL,
                           fontSize = 12,
                           fontFamily = "",
                           legend = "top",
                           borderWidth = 1,
                           allBorders = TRUE,
                           grid = FALSE,
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
        !(condition %in% colnames(colData(mirnaObj)) &
          !condition %in% c("primary", "mirnaCol", "geneCol"))) {
      stop(paste("'condition' must be the column name of a variable specified",
                 "in the metadata (colData) of a MirnaExperiment object; or,",
                 "alternatively, it must be a character/factor object that",
                 "specifies group memberships."),
           call. = FALSE)
    }
  } else {
    if (!is.null(condition)) {
      if ((!is.character(condition) & !is.factor(condition)) |
          length(condition) != nrow(colData(mirnaObj))) {
        stop(paste("'condition' must be the column name of a variable present",
                   "in the metadata (colData) of a MirnaExperiment object; or,",
                   "alternatively, it must be a character/factor object that",
                   "specifies group memberships."),
             call. = FALSE)
      }
    }
  }
  if (!is.numeric(dimensions) |
      length(dimensions) != 2 |
      any(dimensions%%1 != 0) |
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
  if (length(condition) == 1 & !is.null(colorScale)) {
    if ((!is.null(colorScale) & !is.character(colorScale)) |
        any(areColors(colorScale) == FALSE) |
        !identical(
          sort(names(colorScale)),
          sort(unique(colData(mirnaObj)[, condition])))) {
      stop(paste("'colorScale' must be a named character vector where values",
                 "consist of R colors, whereas names coincide to the different",
                 "conditions. For additional details see ?plotDimensions"),
           call. = FALSE)
    }
  } else if (length(condition) != 1 & !is.null(colorScale)) {
    if ((!is.null(colorScale) & !is.character(colorScale)) |
        any(areColors(colorScale) == FALSE) |
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
  if (!is.numeric(fontSize) |
      length(fontSize) != 1 |
      fontSize < 0) {
    stop("'fontSize' must be a non-neagtive number! (default is 12)",
         call. = FALSE)
  }
  if (!is.character(fontFamily) |
      length(fontFamily) != 1) {
    stop("'fontFamily' must be a character of length 1",
         call. = FALSE)
  }
  if (!is.character(legend) |
      length(legend) != 1 |
      !legend %in% c("top", "bottom", "right", "left", "none")) {
    stop("'legend' must be one of 'top', 'bottom' 'right', 'left', and 'none'",
         call. = FALSE)
  }
  if (!is.numeric(borderWidth) |
      length(borderWidth) != 1 |
      borderWidth < 0) {
    stop("'borderWidth' must be a non-neagtive number! (default is 1)",
         call. = FALSE)
  }
  if (!is.logical(allBorders) |
      length(allBorders) != 1) {
    stop("'allBorders' must be logical (TRUE/FALSE)!", call. = FALSE)
  }
  if (!is.logical(grid) |
      length(grid) != 1) {
    stop("'grid' must be logical (TRUE/FALSE)!", call. = FALSE)
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
  samplesMetadata <- colData(mirnaObj)
  meta <- samplesMetadata[!is.na(samplesMetadata[, featCol]), ]
  
  ## reorder metadata based on expression matrix
  meta <- meta[order(match(meta[, featCol], colnames(featExpr))), ]
  
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
  oldCounts <- metadata(mirnaObj)[["oldCounts"]]
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
    Condition <- cond
    mdsPlot <- ggplot2::ggplot(mds, ggplot2::aes(x = .data$x,
                                                 y = .data$y,
                                                 color = Condition,
                                                 label = .data$primary))
  } else {
    mdsPlot <- ggplot2::ggplot(mds, ggplot2::aes(x = .data$x,
                                                 y = .data$y,
                                                 label = .data$primary))
  }
  
  ## add points to the MDS scatterplot
  mdsPlot <- mdsPlot +
    ggplot2::geom_point(alpha = pointAlpha, size = pointSize)
  
  ## add the desired color scale
  if (!is.null(colorScale)) {
    mdsPlot <- mdsPlot +
      ggplot2::scale_color_manual(values = colorScale)
  }
  
  ## apply MIRit ggplot2 theme
  mdsPlot <- mdsPlot +
    theme.MIRit(base_size = fontSize,
                base_family = fontFamily,
                legend = legend,
                borderWidth = borderWidth,
                allBorders = allBorders,
                grid = grid)
  
  ## set axis labels
  mdsPlot <- mdsPlot +
    ggplot2::xlab(paste("Dim ", dimensions[1], " (",
                        round(var.exp[1]*100), "%)", sep = "")) +
    ggplot2::ylab(paste("Dim ", dimensions[2], " (",
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
#' When producing the dotplot with this function, significant pathways
#' are ordered on the x-axis on the basis of their normalized pathway score
#' computed by [topologicalAnalysis()]. The higher is this score, and the more
#' affected a pathway is between biological conditions. Moreover, the size
#' of each dot is equal to the ratio between the number of nodes for which a
#' measurement is available, and the total number of nodes (pathway coverage).
#' Finally, the color scale of dots is relative to the adjusted p-values of
#' each pathway.
#'
#' @param object An object of class
#' [`IntegrativePathwayAnalysis`][IntegrativePathwayAnalysis-class]
#' @param showTerms It is the number of pathways to be shown, based on the
#' order determined by the parameter `showTermsParam`; or, alternatively, a
#' character vector indicating the pathways to plot. Default is `10`
#' @param showTermsParam The order in which the top pathways are selected as
#' specified by the `showTerms` parameter. It must be one of `coverage`,
#' `padj`, `pval`, `score` and `normalized.score` (default)
#' @param title The title of the plot. Default is `NULL` not to include a plot
#' title
#'
#' @returns
#' A `ggplot` graph with a dotplot of integrated pathways.
#'
#' @examples
#' # load example IntegrativePathwayAnalysis object
#' obj <- loadExamples("IntegrativePathwayAnalysis")
#' 
#' # access the results of pathway analysis
#' integratedPathways(obj)
#' 
#' # create a dotplot of integrated pathways
#' integrationDotplot(obj)
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @export
integrationDotplot <- function(object,
                               showTerms = 10,
                               showTermsParam = "normalized.score",
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
  if (!is.character(showTermsParam) |
      length(showTermsParam) != 1 |
      !showTermsParam %in% c("coverage", "padj", "pval", "score",
                             "normalized.score")) {
    stop(paste("'showTermsParam' must be one of: 'coverage', 'padj',",
               "'pval', 'score', 'normalized.score' (default)"),
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
  
  ## rename column p-value column names
  colnames(res)[colnames(res) %in% c("P.Val", "adj.P.Val")] <- c("pval", "padj")
  
  ## order results on the basis of showTermsParam
  if (showTermsParam != "padj" & showTermsParam != "pval") {
    res <- res[order(res[, showTermsParam], decreasing = TRUE), ]
  } else {
    res <- res[order(res[, showTermsParam], decreasing = FALSE), ]
  }
  
  ## select pathways to be shown in the dotplot
  if (is.numeric(showTerms)) {
    res <- res[seq(ifelse(showTerms <= nrow(res), showTerms, nrow(res))), ]
  } else if (is.character(showTerms)) {
    res <- res[which(res$pathway %in% showTerms), ]
  }
  
  ## set the parameters for plotting
  ordBy <- "normalized.score"
  sizeBy <- "coverage"
  colBy <- "padj"
  
  ## set x-axis label
  ordLabel <- "Normalized pathway Score"
  
  ## create a dotplot
  dotRes <- ggplot2::ggplot(res,
                            ggplot2::aes(x = !!ggplot2::sym(ordBy),
                                         y = reorder(.data$pathway,
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
  
  ## add the title of the plot
  if (!is.null(title)) {
    dotRes <- dotRes +
      ggplot2::ggtitle(title)
  }
  
  ## return ggplot2 graph
  return(dotRes)
  
}

