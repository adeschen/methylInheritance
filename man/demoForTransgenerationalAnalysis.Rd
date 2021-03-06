% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylInheritance.R
\docType{data}
\name{demoForTransgenerationalAnalysis}
\alias{demoForTransgenerationalAnalysis}
\title{The methylation information from samples over three generations. Information
for each generation is stored in a
\code{methylRawList} format (for demo purpose).}
\format{
A \code{list} containing three \code{methylRawList} objects. Each
\code{methylRawList} contains the information for one generation
(first entry = first generation, etc..). Each sample information is
stored in a \code{methylRaw} object. There is \code{methylRaw} objects
(6 controls and 6 cases) in each generation.
}
\usage{
data(demoForTransgenerationalAnalysis)
}
\value{
A \code{list} containing three \code{methylRawList} objects. Each
\code{methylRawList} contains the information for one generation
(first entry = first generation, etc..). Each sample information is
stored in a \code{methylRaw} object. There is \code{methylRaw} objects
(6 controls and 6 cases) in each generation.
}
\description{
The object is a \code{list} with 3 entries. Each entry corresponds to the
information for one generation (first entry = first generation, etc..)
stored in a \code{methylRawList} object.
There are 12 samples (6 controls and 6 cases) for each generation. Each
sample information is stored in a \code{methylRaw} object.
}
\details{
This dataset can be used to test \code{runPermutation} and
\code{runObservation} functions.
}
\examples{

## Loading dataset
data(demoForTransgenerationalAnalysis)

## Run a permutation analysis
runObservation(methylKitData = demoForTransgenerationalAnalysis,
    outputDir = "test_demo", type = "tiles", vSeed = 2001)

## Get results
result <- loadAllRDSResults(analysisResultsDir = "test_demo",
    permutationResultsDir = NULL, doingSites = FALSE,
    doingTiles = TRUE)

## Remove result directory
if (dir.exists("test_demo")) {
    unlink("test_demo", recursive = TRUE)
}

}
\seealso{
\itemize{
    \item \code{\link{runPermutation}} {for running a
    permutation analysis, and optionally an observation analysis,
    using multi-generational dataset}
    \item \code{\link{runObservation}} {for running an
    observation analysis using methylKit info entry}
}
}
\keyword{datasets}
