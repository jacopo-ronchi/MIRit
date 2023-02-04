#' List all the available biological pathways in KEGG, Reactome and
#' WikiPathways
#'
#' This function can be used to retrieve a list of valid biological pathways
#' present in KEGG, Reactome and WikiPathways. A valid pathway name is needed
#' to create and visualize miRNA-gene regulation networks through the
#' [mirnaPathway()] function.
#'
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function.
#' @param database The name of the database to use. It must be one of: `KEGG`,
#' `Reactome`, and `WikiPathways`.
#'
#' @returns
#' A `character` vector containing the pathway names present in the
#' specified database.
#'
#' @examples
#' # list the mouse pathways present in WikiPathways
#' listPathways("Mus musculus", "WikiPathways")
#'
#' @note
#' This function uses the [graphite] package to retrieve biological pathways
#' from KEGG, Reactome and WikiPathways.
#'
#' @references
#' Sales, G., Calura, E., Cavalieri, D. et al. graphite - a Bioconductor
#' package to convert pathway topology to gene network.
#' BMC Bioinformatics 13, 20 (2012). [https://doi.org/10.1186/1471-2105-13-20]
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
listPathways <- function(organism, database) {

  ## check inputs
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("KEGG", "Reactome", "WikiPathways")) {
    stop("Databases supported are: 'KEGG', 'Reactome' and 'WikiPathways'",
         call. = FALSE)
  }
  if (database == "KEGG" & !organism %in% convertOrganism("graph_kegg", "all")) {
    stop(paste("For KEGG database 'organism' must be one of:",
               paste(convertOrganism("graph_kegg", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "Reactome" &
             !organism %in% convertOrganism("graph_reactome", "all")) {
    stop(paste("For Reactome database 'organism' must be one of:",
               paste(convertOrganism("graph_reactome", "all"),
                     collapse = ", ")),
         call. = FALSE)
  } else if (database == "WikiPathways" &
             !organism %in% convertOrganism("graph_wikipathways", "all")) {
    stop(paste("For WikiPathways database 'organism' must be one of:",
               paste(convertOrganism("graph_wikipathways", "all"),
                     collapse = ", ")),
         call. = FALSE)
  }

  ## set organism name and database
  database <- tolower(database)
  organism <- convertOrganism(paste("graph", database, sep = "_"), organism)

  ## download pathways from specified database
  pathDb <- graphite::pathways(species = organism, database = database)

  ## return the names of the pathways present in the specified database
  return(names(pathDb))

}





#' Explore the relationships between miRNAs and genes within biological pathways
#'
#' This function is useful to integrate biological pathways with microRNAs. In
#' particular, this function retrieves a user-defined pathway from one of the
#' available databases (KEGG, Reactome and WikiPathways), and integrate every
#' gene of the pathway with the miRNAs that regulate those members. The
#' resulting network can then be used with the [visualizeNetwork] function
#' to visualize the relationships between miRNA and gene expression.
#'
#' @param mirnaObj A [`MirnaExperiment`][MirnaExperiment-class] object
#' containing miRNA and gene data
#' @param pathway The name of a biological pathway to integrate with miRNA
#' expression. It must be a valid pathway name contained in KEGG, Reactome or
#' WikiPathways. For example: `"Phosphorylation of CD3 and TCR zeta chains"`
#' contained in Reactome database
#' @param organism The name of the organism under consideration. The different
#' databases have different supported organisms. To see the list of supported
#' organisms for a given database, use the [supportedOrganisms()] function.
#' Default is `Homo sapiens`
#' @param database The name of the database used for retrieving the biological.
#' pathway. It must be one of: `KEGG`, `Reactome`, and `WikiPathways`.
#' Default is `KEGG`
#'
#' @returns
#' A [`tidygraph::tbl_graph`] object containing the biological pathway
#' integrated with gene and miRNA expression. To plot this network object, you
#' can use the [visualizeNetwork()] function.
#'
#' @seealso
#' For additional details on how to handle [`tidygraph::tbl_graph`] objects,
#' refer to the [tidygraph] package.
#'
#' @examples
#' # perform the integration analysis between miRNAs and their targets
#' obj <- integrateMirnaTargets(obj)
#'
#' # explore a TCR-associated pathway present in Reactome
#' net <- mirnaPathway(obj,
#' pathway = "Phosphorylation of CD3 and TCR zeta chains", organism = "Homo
#' sapiens", database = "Reactome")
#'
#' # visualize miRNA-gene network
#' visualizeNetwork(net)
#'
#' @note
#' This function uses the [graphite] package to retrieve biological pathways as
#' `graph` objects from KEGG, Reactome and WikiPathways.
#'
#' @references
#' Sales, G., Calura, E., Cavalieri, D. et al. graphite - a Bioconductor
#' package to convert pathway topology to gene network.
#' BMC Bioinformatics 13, 20 (2012). [https://doi.org/10.1186/1471-2105-13-20]
#'
#' Minoru Kanehisa, Susumu Goto, KEGG: Kyoto Encyclopedia of Genes and Genomes,
#' Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Pages 27–30,
#' [https://doi.org/10.1093/nar/28.1.27]
#'
#' G. Joshi-Tope, M. Gillespie, I. Vastrik, P. D'Eustachio, E. Schmidt, B.
#' de Bono, B. Jassal, G.R. Gopinath, G.R. Wu, L. Matthews, S. Lewis, E.
#' Birney, L. Stein, Reactome: a knowledgebase of biological pathways, Nucleic
#' Acids Research, Volume 33, Issue suppl_1, 1 January 2005, Pages D428–D432,
#' [https://doi.org/10.1093/nar/gki072]
#'
#' Marvin Martens, Ammar Ammar, Anders Riutta, Andra Waagmeester, Denise N
#' Slenter, Kristina Hanspers, Ryan A. Miller, Daniela Digles, Elisson N Lopes,
#' Friederike Ehrhart, Lauren J Dupuis, Laurent A Winckers, Susan L Coort, Egon
#' L Willighagen, Chris T Evelo, Alexander R Pico, Martina Kutmon,
#' WikiPathways: connecting communities, Nucleic Acids Research, Volume 49,
#' Issue D1, 8 January 2021, Pages D613–D621,
#' [https://doi.org/10.1093/nar/gkaa1024]
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
mirnaPathway <- function(mirnaObj,
                         pathway,
                         database = "KEGG",
                         organism = "Homo sapiens") {

  ## check inputs
  if (!is(mirnaObj, "MirnaExperiment")) {
    stop("'mirnaObj' should be of class MirnaExperiment! See ?MirnaExperiment",
         call. = FALSE)
  }
  if (max(dim(mirnaTargetsIntegration(mirnaObj))) == 0) {
    stop(paste("Integration analysis is not detected in 'mirnaObj'!",
               "Before using this function, expression levels of miRNAs and",
               "genes must be integrated with the 'integrateMirnaTargets()'",
               "function. See '?integrateMirnaTargets' for the details."),
         call. = FALSE)
  }
  if (!is.character(pathway) |
      length(pathway) != 1) {
    stop(paste("'pathway' must be a character of lenght 1 specifying",
               "the pathway name. See ?mirnaPathway. You can also use",
               "the 'listPathways()' function to retrieve all the pathways",
               "present in a given database."),
         call. = FALSE)
  }
  if (!is.character(database) |
      length(database) != 1 |
      !database %in% c("KEGG", "Reactome", "WikiPathways")) {
    stop("Databases supported are: 'KEGG', 'Reactome' and 'WikiPathways'",
         call. = FALSE)
  }
  if (database == "KEGG" & !organism %in% convertOrganism("graph_kegg", "all")) {
    stop(paste("For KEGG database 'organism' must be one of:",
               paste(convertOrganism("graph_kegg", "all"), collapse = ", ")),
         call. = FALSE)
  } else if (database == "Reactome" &
             !organism %in% convertOrganism("graph_reactome", "all")) {
    stop(paste("For Reactome database 'organism' must be one of:",
               paste(convertOrganism("graph_reactome", "all"),
                     collapse = ", ")),
         call. = FALSE)
  } else if (database == "WikiPathways" &
             !organism %in% convertOrganism("graph_wikipathways", "all")) {
    stop(paste("For WikiPathways database 'organism' must be one of:",
               paste(convertOrganism("graph_wikipathways", "all"),
                     collapse = ", ")),
         call. = FALSE)
  }

  ## define integration method used
  integration <- mirnaTargetsIntegration(mirnaObj)
  if (colnames(integration)[5] == "Fisher.P.Val") {
    intMeth <- "fisher"
  } else {
    intMeth <- "correlation"
  }

  ## set organism name and database
  database <- tolower(database)
  organism <- convertOrganism(paste("graph", database, sep = "_"), organism)

  ## download pathways from specified database
  pathDb <- graphite::pathways(species = organism, database = database)

  ## check if the pathway name provided is valid
  if (!pathway %in% names(pathDb)) {
    stop(paste("The selected 'pathway' is not present in ", database, "! ",
               "To retrieve the list of biological pathways available in ",
               "this database use the 'listPathways()' function.",
               sep = ""),
         call. = FALSE)
  }

  ## extract the selected pathway
  pathObj <- pathDb[[pathway]]
  pathObj <- graphite::convertIdentifiers(pathObj, "SYMBOL")
  pathGraph <- graphite::pathwayGraph(pathway = pathObj)
  graph::nodes(pathGraph) <- gsub("SYMBOL:", "", graph::nodes(pathGraph))

  ## extract integrated-targets
  if (intMeth == "fisher") {
    intTarg <- vapply(integration$DE_targets, stringr::str_split, pattern = "/",
                      list(1))
    targLengths <- vapply(intTarg, length, 1)
    targets <- data.frame("mirna" = rep(integration$microRNA, targLengths),
                          "targets" = unlist(intTarg), row.names = NULL)
    targets$correlation <- "negative"
  } else if (intMeth == "correlation") {
    targets <- integration[, c(1, 2, 4)]
    colnames(targets) <- c("mirna", "targets", "correlation")
  }

  ## keep only targets involved in the specified pathway
  pathTargs <- targets[targets$targets %in% graph::nodes(pathGraph), ]

  ## add miRNA-target pairs to the biological network
  pathGraph <- graph::addNode(unique(pathTargs$mirna), pathGraph)
  pathGraph <- graph::addEdge(pathTargs$mirna, pathTargs$targets, pathGraph)

  ## convert graphNEL to a tbl_graph object
  tblGraph <- tidygraph::as_tbl_graph(pathGraph)

  ## add expression fold changes to nodes in the pathway
  dem <- mirnaDE(mirnaObj, onlySignificant = FALSE)
  deg <- geneDE(mirnaObj, onlySignificant = FALSE)
  graphNodes <- as.data.frame(tidygraph::activate(tblGraph, nodes))
  graphNodes$type <- ifelse(graphNodes$name %in% pathTargs$mirna,
                            "miRNA",
                            "gene")
  graphNodes$mirnaLogFC <- dem$logFC[match(graphNodes$name, dem$ID)]
  graphNodes$geneLogFC <- deg$logFC[match(graphNodes$name, deg$ID)]
  graphNodes$geneLogFC[is.na(graphNodes$geneLogFC) &
                         graphNodes$type != "miRNA"] <- 0

  ## add integration status to nodes in the pathway
  graphNodes$integrated <- ifelse(graphNodes$name %in% pathTargs$mirna |
                                    graphNodes$name %in% pathTargs$targets,
                                  "integrated",
                                  "not-integrated")

  ## rebuild the tidygraph object
  tblGraph <- tidygraph::left_join(tblGraph, graphNodes, by = "name")

  ## add the direction of miRNA-target interaction
  graphEdges <- as.data.frame(tidygraph::activate(tblGraph, edges))

  graphEdges$edgeType[graphEdges$edgeType == "undefined"] <- mapply(
    function(fromNode, toNode) {
      origin <- graphNodes$name[fromNode]
      destination <- graphNodes$name[toNode]
      cType <- pathTargs$correlation[pathTargs$mirna == origin &
                                       pathTargs$targets == destination]
      if (cType == "negative") {
        "INHIBITION"
      } else if (cType == "positive") {
        "ACTIVATION"
      }
    },
    graphEdges$from[graphEdges$edgeType == "undefined"],
    graphEdges$to[graphEdges$edgeType == "undefined"])

  tblGraph <- tidygraph::mutate(tidygraph::activate(tblGraph, edges),
                                edgeType = graphEdges$edgeType)

  ## return tidygraph object
  return(tblGraph)

}





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
#' implemented in the [ggraph] package.
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
#' pathway)". Default is `NULL` not to include a plot title
#' @param edgesCol It must be an R color name that specifies the color of
#' interaction arrows in the network. Default is `gray`. All available colors
#' can be listed with [colors()]
#' @param bindingEdgesCol It must be an R color name that specifies the color
#' of chemical binding arrows in the network. Default is `gray`. All available
#' colors can be listed with [colors()]
#' @param stress.text.distance The distance between activation/inhibition
#' arrow caps and nodes for the `stress` design. Default is `0.2`
#' @param circle.text.distance The distance between activation/inhibition
#' arrow caps and nodes for the `circle` design. Default is `0.1`
#'
#' @returns
#' A `ggraph` object representing the network graph. For instructions on how
#' to handle and edit this object, please refer to the [ggraph] package.
#'
#' @examples
#' # explore a TCR-associated pathway present in Reactome
#' net <- mirnaPathway(obj,
#' pathway = "Phosphorylation of CD3 and TCR zeta chains", organism = "Homo
#' sapiens", database = "Reactome", onlyIntegrated = TRUE)
#'
#' # visualize miRNA-gene network
#' visualizeNetwork(net)
#'
#' # visualize network with a circular layout and with points
#' visualizeNetwork(net, layout = "circle", node = "point")
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
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
  if (!is.character(title) |
      length(title) != 1) {
    stop(paste("'title' must be the title of the plot (e.g. 'Apoptosis').",
               "For additional details see ?visualizeNetwork"),
         call. = FALSE)
  }
  if (!is.character(edgesCol) |
      length(edgesCol) != 1 |
      !edgesCol %in% colors()) {
    stop(paste("'edgesCol' must be an R color name. All available colors",
               "can be listed with 'colors()'."),
         call. = FALSE)
  }
  if (!is.character(bindingEdgesCol) |
      length(bindingEdgesCol) != 1 |
      !bindingEdgesCol %in% colors()) {
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

  ## set accepted pathway processes
  activation <- "ACTIVATION|Activation|activation|expression"
  inhibition <- "INHIBITION|Inhibition|inhibition"
  binding <- "Binding|binding"

  ## retain only integrated targets if wanted by the user
  if (onlyIntegrated == TRUE) {
    pathGraph <- tidygraph::to_subgraph(pathGraph,
                                        integrated == "integrated",
                                        subset_by = "nodes")$subgraph
  }

  ## create a ggraph network to visualize biological pathways with DE-miRNAs
  mirNet <- ggraph::ggraph(pathGraph, layout = layout) +

    ## use different edges for different biological interactions
    ggraph::geom_edge_link(
      ggplot2::aes(filter = grepl(binding, edgeType),
                   linetype = "binding"),
      color = bindingEdgesCol) +

    ggraph::geom_edge_link(
      ggplot2::aes(filter = grepl(inhibition, edgeType) &
                     !grepl(activation, edgeType),
                   linetype = "inhibition",
                   end_cap = {
                     ## set the arrow distance from nodes for boxed labels/points
                     if (node == "box") {
                       ggraph::label_rect(node2.name,
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
      ggplot2::aes(filter = grepl(activation, edgeType),
                   linetype = "activation",
                   end_cap = {
                     ## set the arrow distance from nodes for boxed labels/points
                     if (node == "box") {
                       ggraph::label_rect(node2.name,
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
                                         sep = "|"), edgeType),
                   linetype = "standard"),
      color = edgesCol)

  ## customize colors and labels with different scales for genes and miRNAs
  if (node == "box") {
    mirNet <- mirNet +
      ggraph::geom_node_label(ggplot2::aes(label = name, fill = geneLogFC), size = size) +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                    guide = ggplot2::guide_colorbar(order = 1)) +
      ggnewscale::new_scale_fill() +
      ggraph::geom_node_label(data = function(x) {tidygraph::filter(x, is.na(geneLogFC))},
                              ggplot2::aes(label = name, fill = mirnaLogFC), size = size) +
      ggplot2::scale_fill_gradient2(low = "green", mid = "grey", high = "yellow",
                                    guide = ggplot2::guide_colorbar(order = 2))
  } else if (node == "point") {
    mirNet <- mirNet +
      ggraph::geom_node_point(ggplot2::aes(fill = geneLogFC),
                              pch = 21, size = size) +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                    guide = ggplot2::guide_colorbar(order = 1)) +
      ggnewscale::new_scale_fill() +
      ggraph::geom_node_point(data = function(x) {tidygraph::filter(x, is.na(geneLogFC))},
                              ggplot2::aes(fill = mirnaLogFC),
                              pch = 21, size = size) +
      ggplot2::scale_fill_gradient2(low = "green", mid = "grey", high = "yellow",
                                    guide = ggplot2::guide_colorbar(order = 2))
  }

  ## add text elements for point nodes
  if (node == "point") {
    if (layout == "stress") {
      mirNet <- mirNet +
        ggraph::geom_node_text(ggplot2::aes(label = name),
                               size = size, nudge_y = - stress.text.distance)
    } else if (layout == "circle") {
      mirNet <- mirNet +
        ggraph::geom_node_text(ggplot2::aes(label = name), size = size,
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
                   colBy) {
            mirnaDotplot.MirnaEnrichment(object,
                                         showTerms = showTerms,
                                         splitDir = splitDir,
                                         ordBy = ordBy,
                                         sizeBy = sizeBy,
                                         colBy = colBy)
          })





#' @rdname mirnaDotplot
#' @export
setMethod("mirnaDotplot", "MirnaGsea",
          function(object,
                   showTerms,
                   splitDir,
                   ordBy,
                   sizeBy,
                   colBy) {
            mirnaDotplot.MirnaGsea(object,
                                   showTerms = showTerms,
                                   splitDir = splitDir,
                                   ordBy = ordBy,
                                   sizeBy = sizeBy,
                                   colBy = colBy)
          })





## create a dotplot for MirnaEnrichment objects
mirnaDotplot.MirnaEnrichment <- function(mirnaEnr,
                                         showTerms,
                                         splitDir,
                                         ordBy,
                                         sizeBy,
                                         colBy) {

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
                                         y = reorder(Subcategory,
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

  ## return ggplot2 graph
  return(dotRes)

}





## create a dotplot for MirnaGsea objects
mirnaDotplot.MirnaGsea <- function(mirnaGsea,
                                   showTerms,
                                   splitDir,
                                   ordBy,
                                   colBy,
                                   sizeBy) {

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
                                         y = reorder(Subcategory,
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
#' @param mirnaGsea A [`MirnaGsea`] object containing the results of gene set
#' enrichment analysis obtained with the [gseaMirnas()] function
#' @param showTerms It is the number of top enriched terms to show or,
#' alternatively, a character vector indicating the enriched terms to plot.
#' Default is `10`
#' @param colBy The parameter used to set the color scale. It must be one of
#' `P.adjusted` (default), `P.value` and `Observed`
#'
#' @returns
#' An object of class `ggplot` containing the ridgeplot of miRNA GSEA results.
#'
#' @examples
#' # perform miRNA GSEA analysis
#' gse_res <- gseaMirnas(obj, organism = "Homo sapiens",
#' category = "miRWalk_GO_mature")
#'
#' # plot results as a ridgeplot
#' mirnaRidgeplot(gse_res)
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
mirnaRidgeplot <- function(mirnaGsea, showTerms = 10, colBy = "P.adjusted") {

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
  ridgeRes <- ggplot2::ggplot(ridgeDf, ggplot2::aes(x = val,
                                                    y = reorder(term, meanLFC),
                                                    fill = !!ggplot2::sym(colBy))) +
    ggridges::geom_density_ridges() +
    ggplot2::scale_fill_continuous(low = "red",
                                   high = "blue",
                                   name = "P.adjusted",
                                   guide = ggplot2::guide_colorbar(reverse = TRUE)) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab("miRNAs logFC") +
    theme_enr()

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
                                         margin = ggplot2::margin(10, 5, 0, 0))
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
#' with [colors()]
#' @param mirFill It must be an R color name that specifies the fill color of
#' the miRNA locus. Default is `orange`. All available colors can be listed
#' with [colors()]
#' @param ... Other parameters that can be passed to [Gviz::plotTracks()]
#' function
#'
#' @note
#' This function retrieves genomic coordinates from the output of
#' [findMirnaSNPs()] function and then uses [Gviz] package to build
#' the trackplot.
#'
#' @references
#' Hahne, F., Ivanek, R. (2016). Visualizing Genomic Data Using Gviz and
#' Bioconductor. In: Mathé, E., Davis, S. (eds) Statistical Genomics. Methods
#' in Molecular Biology, vol 1418. Humana Press, New York, NY.
#' [https://doi.org/10.1007/978-1-4939-3578-9_16]
#'
#' @returns
#' A trackplot with information about chromosome, SNP and miRNA gene location.
#'
#' @examples
#' # retrieve associated SNPs
#' association <- findMirnaSNPs(obj, disId)
#'
#' # visualize association as a trackplot
#' mirVariantPlot(variantId = varId, snpAssociation = association)
#'
#' @author
#' Jacopo Ronchi, \email{j.ronchi2@@campus.unimib.it}
#'
#' @export
mirVariantPlot <- function(variantId,
                           snpAssociation,
                           showContext = FALSE,
                           showSequence = TRUE,
                           snpFill = "lightblue",
                           mirFill = "orange",
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
      !snpFill %in% colors()) {
    stop(paste("'snpFill' must be an R color name. All available colors",
               "can be listed with 'colors()'. Default is: 'lightblue'"),
         call. = FALSE)
  }
  if (!is.character(mirFill) |
      length(mirFill) != 1 |
      !mirFill %in% colors()) {
    stop(paste("'mirFill' must be an R color name. All available colors",
               "can be listed with 'colors()'. Default is: 'orange'"),
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
                                ...)


}

