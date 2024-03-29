% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{enrichmentDatabase}
\alias{enrichmentDatabase}
\title{Access the database used for functional enrichment analyses}
\usage{
enrichmentDatabase(object)
}
\arguments{
\item{object}{An object of class
\code{\link[=FunctionalEnrichment-class]{FunctionalEnrichment}} containing
enrichment results}
}
\value{
A \code{character} containing the name of the database, such as \code{KEGG}.
}
\description{
This function accesses the \code{database} slot of a
\code{\link[=FunctionalEnrichment-class]{FunctionalEnrichment}} object and returns
a the name of the database used by the \code{\link[=enrichGenes]{enrichGenes()}} function to perform
the enrichment analysis.
}
\examples{
# load example FunctionalEnrichment object
obj <- loadExamples("FunctionalEnrichment")

# see the database
enrichmentDatabase(obj)

}
\author{
Jacopo Ronchi, \email{jacopo.ronchi@unimib.it}
}
