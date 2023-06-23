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
#' This function uses the `graphite` package to retrieve biological pathways
#' from KEGG, Reactome and WikiPathways.
#'
#' @references
#' Sales, G., Calura, E., Cavalieri, D. et al. graphite - a Bioconductor
#' package to convert pathway topology to gene network.
#' BMC Bioinformatics 13, 20 (2012).
#' \url{https://doi.org/10.1186/1471-2105-13-20}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
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
#' refer to the `tidygraph` package.
#'
#' @examples
#' # load example MirnaExperiment object
#' obj <- loadExamples()
#' 
#' # perform the integration analysis between miRNAs and their targets
#' obj <- mirnaIntegration(obj)
#'
#' # explore Thyroid hormone synthesis pathway in KEGG
#' net <- mirnaPathway(obj,
#' pathway = "Thyroid hormone synthesis", organism = "Homo sapiens",
#' database = "KEGG")
#'
#' # visualize miRNA-gene network
#' visualizeNetwork(net)
#'
#' @note
#' This function uses the `graphite` package to retrieve biological pathways as
#' `graph` objects from KEGG, Reactome and WikiPathways.
#'
#' @references
#' Sales, G., Calura, E., Cavalieri, D. et al. graphite - a Bioconductor
#' package to convert pathway topology to gene network.
#' BMC Bioinformatics 13, 20 (2012).
#' \url{https://doi.org/10.1186/1471-2105-13-20}
#'
#' Minoru Kanehisa, Susumu Goto, KEGG: Kyoto Encyclopedia of Genes and Genomes,
#' Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Pages 27–30,
#' \url{https://doi.org/10.1093/nar/28.1.27}
#'
#' G. Joshi-Tope, M. Gillespie, I. Vastrik, P. D'Eustachio, E. Schmidt, B.
#' de Bono, B. Jassal, G.R. Gopinath, G.R. Wu, L. Matthews, S. Lewis, E.
#' Birney, L. Stein, Reactome: a knowledgebase of biological pathways, Nucleic
#' Acids Research, Volume 33, Issue suppl_1, 1 January 2005, Pages D428–D432,
#' \url{https://doi.org/10.1093/nar/gki072}
#'
#' Marvin Martens, Ammar Ammar, Anders Riutta, Andra Waagmeester, Denise N
#' Slenter, Kristina Hanspers, Ryan A. Miller, Daniela Digles, Elisson N Lopes,
#' Friederike Ehrhart, Lauren J Dupuis, Laurent A Winckers, Susan L Coort, Egon
#' L Willighagen, Chris T Evelo, Alexander R Pico, Martina Kutmon,
#' WikiPathways: connecting communities, Nucleic Acids Research, Volume 49,
#' Issue D1, 8 January 2021, Pages D613–D621,
#' \url{https://doi.org/10.1093/nar/gkaa1024}
#'
#' @author
#' Jacopo Ronchi, \email{jacopo.ronchi@@unimib.it}
#'
#' @importFrom rlang .data
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
  if (max(dim(integration(mirnaObj))) == 0) {
    stop(paste("Integration analysis is not detected in 'mirnaObj'!",
               "Before using this function, expression levels of miRNAs and",
               "genes must be integrated with the 'mirnaIntegration()'",
               "function. See '?mirnaIntegration' for the details."),
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
  
  ## listing the variables not present in the global environment
  edges <- nodes <- NULL
  
  ## define integration method used
  integration <- integration(mirnaObj)
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

