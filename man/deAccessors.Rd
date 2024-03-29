% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{deAccessors}
\alias{deAccessors}
\alias{mirnaDE}
\alias{geneDE}
\title{Extract differentially expressed miRNAs and genes from a
\code{\link[=MirnaExperiment-class]{MirnaExperiment}} object}
\usage{
mirnaDE(object, onlySignificant = TRUE, param = FALSE, returnObject = FALSE)

geneDE(object, onlySignificant = TRUE, param = FALSE, returnObject = FALSE)
}
\arguments{
\item{object}{A \code{\link[=MirnaExperiment-class]{MirnaExperiment}} object containing
miRNA and gene data}

\item{onlySignificant}{Logical, if \code{TRUE} differential expression results
will be returned just for statistically significant miRNAs/genes, if \code{FALSE}
the full table of miRNA/gene differential expression will be provided.
Default is \code{TRUE} to only report significant miRNAs/genes}

\item{param}{Logical, whether to return the complete \code{list} object with
the parameters used, or just the results stored in \code{data}. Default is FALSE}

\item{returnObject}{Logical, if \code{TRUE} this function will return the
limma/edgeR/DESeq2 object used for differential expression analysis}
}
\value{
A \code{data.frame} with miRNA/gene differential expression, or a \code{list} object
with the parameters used if \code{param = TRUE}.
}
\description{
The \code{mirnaDE()} and \code{geneDE()} are two accessor functions for the \code{mirnaDE}
and \code{geneDE} slots of \code{\link[=MirnaExperiment-class]{MirnaExperiment}} class,
respectively. Thus, they can be used to explore the results of miRNA and
gene differential expression analysis stored in a
\code{\link[=MirnaExperiment-class]{MirnaExperiment}} object.
}
\section{Functions}{
\itemize{
\item \code{mirnaDE()}: Extract miRNA differential expression results

\item \code{geneDE()}: Extract gene differential expression results

}}
\examples{
# load example MirnaExperiment object
obj <- loadExamples()

# access miRNA differential expression of a MirnaExperiment object
sig <- mirnaDE(obj)
all <- mirnaDE(obj, onlySignificant = FALSE)

# access gene differential expression of a MirnaExperiment object
sig <- geneDE(obj)
all <- geneDE(obj, onlySignificant = FALSE)

}
\author{
Jacopo Ronchi, \email{jacopo.ronchi@unimib.it}
}
