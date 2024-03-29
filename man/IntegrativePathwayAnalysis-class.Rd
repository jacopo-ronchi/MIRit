% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R, R/show-methods.R
\docType{class}
\name{IntegrativePathwayAnalysis-class}
\alias{IntegrativePathwayAnalysis-class}
\alias{IntegrativePathwayAnalysis}
\alias{integratedPathways,IntegrativePathwayAnalysis-method}
\alias{integrationDatabase,IntegrativePathwayAnalysis-method}
\alias{augmentedPathways,IntegrativePathwayAnalysis-method}
\alias{show,IntegrativePathwayAnalysis-method}
\title{The \code{IntegrativePathwayAnalysis} class}
\usage{
\S4method{integratedPathways}{IntegrativePathwayAnalysis}(object)

\S4method{integrationDatabase}{IntegrativePathwayAnalysis}(object)

\S4method{augmentedPathways}{IntegrativePathwayAnalysis}(object)

\S4method{show}{IntegrativePathwayAnalysis}(object)
}
\arguments{
\item{object}{An object of class
\code{\link[=IntegrativePathwayAnalysis-class]{IntegrativePathwayAnalysis}} containing
the results of a miRNA-mRNA pathway analysis}
}
\value{
\itemize{
\item the \code{integratedPathways()} function returns a \code{data.frame} with the
results of the integrative pathway analysis.
}

\itemize{
\item the \code{integrationDatabase()} function returns a \code{character} object
with the database used for the analysis.
}

\itemize{
\item the \code{augmentedPathways()} function returns a \code{list} object
containing the miRNA-augmented pathways considered for TAIPA.
}
}
\description{
This class stores the output of integrative multi-omic pathway analyses.
In particular, the slots of this class are suitable to represent the results
of topologically-aware integrative pathway analysis (TAIPA) returned from
the \code{\link[=topologicalAnalysis]{topologicalAnalysis()}} function.
}
\details{
\subsection{Analysis results}{

The \code{data} slot of this class consists in a \code{data.frame} object with six
columns, namely:
\itemize{
\item \code{pathway}, which indicates the name of the biological network;
\item \code{coverage}, which specifies the fraction of nodes with expression
measurement available;
\item \code{score}, which expresses the score of each individual pathway;
\item \code{normalized.score}, which indicates the pathway scores after standardizing
the values for the null distribution computed through permutations;
\item \code{P.Val}, the resulting p-value of each pathway;
\item \code{adj.P.Val}, the p-value adjusted for multiple testing.
}
}

\subsection{Organisms and databases}{

The \code{organism} and \code{database} slots specify the organism in study and the
database used for retrieving biological interactions, respectively. In
particular, the \code{\link[=topologicalAnalysis]{topologicalAnalysis()}} function supports \code{KEGG},
\code{WikiPathways}, and \code{Reactome} databases. Regarding organisms, the
\code{\link[=supportedOrganisms]{supportedOrganisms()}} function can be used to retrieve the available
species for each database.
}

\subsection{Statistical significance of the permutation test}{

\code{pCutoff} and \code{pAdjustment} slots refer to the cutoff used for the analysis.
\code{pCutoff} is the threshold used for defining statistically significant
pathways, whereas \code{pAdjustment} refers to the multiple testing correction
method used.

Furthermore, since the statistical significance of each pathway is defined
on the basis of a permutation test, the number of permutations is also
specified in the \code{nPerm} slot.
}

\subsection{Augmented pathways}{

The \code{pathways} slot contains a \code{list} with weighted \code{graph} objects, each
representing a biological pathway. These networks have been enlarged by
adding the observed miRNA-mRNA interactions. Each network has been
processed so that the weight of each edge is +1 for activation interactions,
and -1 for repression interactions, such as those occurring between miRNAs
and mRNAs.
}

\subsection{Differential expression results for both miRNAs and genes}{

The expression variation of all miRNAs and genes measured in the study is
stored in the \code{expression} slot. In particular, this slot consists of a
\code{data.frame} object with different information, including log2 fold changes,
node weights and p-values.
}

\subsection{Minimum percentage of measured features}{

The \code{minPc} slot indicates the minimum percentage of miRNAs/mRNAs above
which pathways have been considered for the integrative analysis. This is
needed because often, when differential expression analysis is performed,
lowly expressed features are removed. Therefore, some pathways might result
significantly affected even if only 1\% of nodes is perturbed.
}
}
\section{Functions}{
\itemize{
\item \code{integratedPathways(IntegrativePathwayAnalysis)}: Access the results of
integrative miRNA-mRNA pathway analysis

\item \code{integrationDatabase(IntegrativePathwayAnalysis)}: View the database used for
the integrative pathway analysis

\item \code{augmentedPathways(IntegrativePathwayAnalysis)}: Extract the list of biological
networks augmented with miRNA-mRNA interactions

\item \code{show(IntegrativePathwayAnalysis)}: Show method for objects of
class IntegrativePathwayAnalysis

}}
\section{Slots}{

\describe{
\item{\code{data}}{A \code{data.frame} object that contains the results of the
integrative pathway analysis. See the \emph{details} section for further details}

\item{\code{method}}{The method used for the analysis}

\item{\code{organism}}{The name of the organism under consideration (e.g.
\verb{Homo sapiens})}

\item{\code{database}}{The name of the database used for retrieving biological
pathways (e.g. \code{KEGG})}

\item{\code{pCutoff}}{A \code{numeric} value defining the threshold used for
statistical significance (e.g. \code{0.05})}

\item{\code{pAdjustment}}{A \code{character} indicating the method used to correct
p-values for multiple testing (e.g. \code{fdr})}

\item{\code{pathways}}{A \code{list} of \code{graph} objects containing the biological
networks retrieved from \code{database}, and augmented with
miRNA-mRNA interactions}

\item{\code{expression}}{A \code{data.frame} object containing differential expression
results for both miRNAs and genes}

\item{\code{minPc}}{The minimum percentage of measured features that a pathway must
have for being considered in the analysis}

\item{\code{nPerm}}{The number of permutation used for assessing the statistical
significance of each pathway}
}}

\author{
Jacopo Ronchi, \email{jacopo.ronchi@unimib.it}
}
