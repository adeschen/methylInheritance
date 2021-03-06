% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylInheritanceMethods.R
\name{loadConvergenceData}
\alias{loadConvergenceData}
\title{Load convergence information from RDS files}
\usage{
loadConvergenceData(
  analysisResultsDir,
  permutationResultsDir,
  type = c("sites", "tiles"),
  inter = c("i2", "iAll"),
  position,
  by = 100
)
}
\arguments{
\item{analysisResultsDir}{a \code{character} string, the path to the
directory that contains the analysis results. The path can be the same as
for the \code{permutatioNResultsDir} parameter.}

\item{permutationResultsDir}{a \code{character} string, the path to the
directory that contains the permutation results. The path can be the same
as for the \code{analysisResultsDir} parameter.}

\item{type}{One of the \code{"sites"} or \code{"tiles"} strings.
Specifies the type
of differentially methylated elements should be returned. For
retrieving differentially methylated bases \code{type} = \code{"sites"}; for
differentially methylated regions \code{type} = \code{"tiles"}.
Default: \code{"sites"}.}

\item{inter}{One of the \code{"i2"} or \code{"iAll"} strings. Specifies the
type of intersection should be returned. For
retrieving intersection results between two consecutive generations
\code{inter} = \code{"i2"}; for intersection results between three
generations or more \code{inter} = \code{"iAll"}.
Default: \code{"i2"}.}

\item{position}{a positive \code{integer}, the position in the \code{list}
where the information will be extracted.}

\item{by}{a \code{integer}, the increment of the number of permutations
where the significant level is tested. Default: 100.}
}
\value{
a graph showing the evolution of the significant level with the
number of permutations
}
\description{
Load convergence information from RDS files.
}
\examples{

## Get the name of the directory where files are stored
filesDir <- system.file("extdata", "TEST", package="methylInheritance")

## Load convergence information
results <- loadConvergenceData(analysisResultsDir = filesDir,
    permutationResultsDir = filesDir, type="sites", inter="i2", position=1,
    by=1)

}
\author{
Astrid Deschenes, Pascal Belleau
}
