#' methylInheritance: Permutation-Based Analysis associating Conserved
#' Differentially Methylated Elements from One Generation to the Next to
#' a Treatment Effect
#'
#' This package does a permutation analysis, based on Monte Carlo sampling,
#' for testing the hypothesis that the number of conserved differentially
#' methylated elements (sites or tiles), between
#' several generations, is associated to an effect inherited from a treatment
#' and that stochastic effect can be dismissed.
#'
#'
#' @name methylInheritance-package
#'
#' @aliases methylInheritance-package methylInheritance
#'
#' @author Astrid Deschênes,
#' Pascal Belleau and
#' Arnaud Droit
#'
#' Maintainer:
#' Astrid Deschenes <adeschen@hotmail.com>
#'
#' @seealso
#' \itemize{
#'     \item{\code{\link{runPermutation}}}{ for running a
#'     permutation analysis, and optionally an observation analysis, on a
#'     specified multi-generational dataset}
#'     \item{\code{\link{runObservation}}}{ for running an
#'     observation analysis on a specified multi-generational dataset}
#' }
#'
#' @return methylInheritance
#' @encoding UTF-8
#' @keywords package
#' @keywords internal
"_PACKAGE"

#' All samples information, formated by \code{methylKit}, in a
#' \code{methylRawList} format (for demo purpose).
#'
#' The object is a \code{list} with 3 entries. Each entry corresponds to the
#' information for one generation (first entry = first generation, etc..)
#' stored in a \code{methylRawList}.
#' There are 12 samples (6 controls and 6 cases) for each generation. Each
#' sample information is stored in a \code{methylRaw} object.
#'
#' This dataset can be
#' used to test the \code{runPermutation} function.
#'
#' @name samplesForTransgenerationalAnalysis
#'
#' @docType data
#'
#' @aliases samplesForTransgenerationalAnalysis
#'
#' @format A \code{list} containing three \code{methylRawList} objects. Each
#' \code{methylRawList} contains the information for one generation
#' (first entry = first generation, etc..). Each sample information is
#' stored in a \code{methylRaw} object. There is \code{methylRaw} objects
#' (6 controls and 6 cases) in each generation.
#'
#' @return A \code{list} containing three \code{methylRawList} objects. Each
#' \code{methylRawList} contains the information for one generation
#' (first entry = first generation, etc..). Each sample information is
#' stored in a \code{methylRaw} object. There is \code{methylRaw} objects
#' (6 controls and 6 cases) in each generation.
#'
#' @seealso
#' \itemize{
#'     \item{\code{\link{runPermutation}}}{ for running a
#'     permutation analysis, and optionally an observation analysis, using
#'     multi-generational dataset}
#' }
#'
#' @usage data(samplesForTransgenerationalAnalysis)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## Run a permutation analysis
#' runPermutation(methylKitData = samplesForTransgenerationalAnalysis,
#'     type = "tiles", nbrPermutations = 2, vSeed = 2332)
#'
NULL

#' The methylation information from samples over three generations. Information
#' for each generation is stored in a
#' \code{methylRawList} format (for demo purpose).
#'
#' The object is a \code{list} with 3 entries. Each entry corresponds to the
#' information for one generation (first entry = first generation, etc..)
#' stored in a \code{methylRawList} object.
#' There are 12 samples (6 controls and 6 cases) for each generation. Each
#' sample information is stored in a \code{methylRaw} object.
#'
#' This dataset can be used to test \code{runPermutation} and
#' \code{runObservation} functions.
#'
#' @name demoForTransgenerationalAnalysis
#'
#' @docType data
#'
#' @aliases demoForTransgenerationalAnalysis
#'
#' @format A \code{list} containing three \code{methylRawList} objects. Each
#' \code{methylRawList} contains the information for one generation
#' (first entry = first generation, etc..). Each sample information is
#' stored in a \code{methylRaw} object. There is \code{methylRaw} objects
#' (6 controls and 6 cases) in each generation.
#'
#' @return A \code{list} containing three \code{methylRawList} objects. Each
#' \code{methylRawList} contains the information for one generation
#' (first entry = first generation, etc..). Each sample information is
#' stored in a \code{methylRaw} object. There is \code{methylRaw} objects
#' (6 controls and 6 cases) in each generation.
#'
#' @seealso
#' \itemize{
#'     \item{\code{\link{runPermutation}}}{ for running a
#'     permutation analysis, and optionally an observation analysis,
#'     using multi-generational dataset}
#'     \item{\code{\link{runObservation}}}{ for running an
#'     observation analysis using methylKit info entry}
#' }
#'
#' @usage data(demoForTransgenerationalAnalysis)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(demoForTransgenerationalAnalysis)
#'
#' ## Run a permutation analysis
#' runObservation(methylKitData = demoForTransgenerationalAnalysis,
#'     outputDir = "test_demo", type = "tiles", vSeed = 2001)
#'
#' ## Get results
#' result <- loadAllRDSResults(analysisResultsDir = "test_demo",
#'     permutationResultsDir = NULL, doingSites = FALSE,
#'     doingTiles = TRUE)
#'
#' ## Remove result directory
#' if (dir.exists("test_demo")) {
#'     unlink("test_demo", recursive = TRUE)
#' }
#'
NULL


#' All observed and permutation results formatted in a
#' \code{methylInheritanceResults} class (for demo purpose).
#'
#' The object is a \code{list} with 2 entries: "OBSERVATION" and
#' "PERMUTATION".
#'
#' This dataset can be
#' used to test the \code{extractInfo} function.The extracted information can
#' be used to calculate the significant level or to create a graph.
#'
#' @name methylInheritanceResults
#'
#' @docType data
#'
#' @aliases methylInheritanceResults
#'
#' @format a \code{list} of class \code{methylInheritanceAllResults}
#' containing the following elements:
#' \itemize{
#' \item \code{OBSERVATION} a \code{list} containing:
#' \itemize{
#' \item \code{SITES} a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hyper differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hypo differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 1 entry,
#'the number of conserved
#' hyper differentially methylated sites between the three consecutive
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 1 entry, the number
#' of conserved hypo differentially methylated sites between the three
#' consecutive generations.
#' }
#' }
#' \item \code{TILES} a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hyper differentially methylated positions between two consecutive
#' generations. The first element represents the intersection of the
#' first and second generations; the second element, the intersection of
#' the second and third generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hypo differentially methylated positions between two consecutive
#' generations.The first element represents the intersection of the first and
#' second generations; the second element, the intersection of the second
#' and third generations.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hyper differentially methylated positions between the three consecutive
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hypo differentially methylated positions between the three consecutive
#' generations.
#' }
#' }
#' }
#' \item \code{PERMUTATION} a \code{list}
#' containing \code{nbrPermutations} entries. Each entry is
#' a \code{list} containing:
#' \itemize{
#' \item \code{SITES} a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hyper differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hypo differentially methylated sites between two consecutive generations.The
#' first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hyper differentially methylated sites between the three consecutive
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hypo differentially methylated sites between the three consecutive
#' generations.
#' }
#' }
#' \item \code{TILES} a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hyper differentially methylated positions between two consecutive
#' generations. The first element represents the intersection of the
#' first and second generations; the second element, the intersection of
#' the second and third generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hypo differentially methylated positions between two consecutive
#' generations.The first element represents the intersection of the first and
#' second generations; the second element, the intersection of the second
#' and third generations.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hyper differentially methylated positions between the three consecutive
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hypo differentially methylated positions between the three consecutive
#' generations.
#' }
#' }
#' }
#' }
#'
#' @return a \code{list} of class \code{methylInheritanceAllResults}
#' containing the following elements:
#' \itemize{
#' \item \code{OBSERVATION} a \code{list} containing:
#' \itemize{
#' \item \code{SITES} a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hyper differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hypo differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations.
#' }
#' \item \code{iAll} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 1 entry,
#'the number of conserved
#' hyper differentially methylated sites between the three consecutive
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 1 entry, the number
#' of conserved hypo differentially methylated sites between the three
#' consecutive generations.
#' }
#' }
#' \item \code{TILES} a \code{list} containing:
#' \itemize{
#' \item \code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hyper differentially methylated positions between two consecutive
#' generations. The first element represents the intersection of the
#' first and second generations; the second element, the intersection of
#' the second and third generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hypo differentially methylated positions between two consecutive
#' generations.The first element represents the intersection of the first and
#' second generations; the second element, the intersection of the second
#' and third generations.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hyper differentially methylated positions between the three consecutive
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hypo differentially methylated positions between the three consecutive
#' generations.
#' }
#' }
#' }
#' \item \code{PERMUTATION} a \code{list}
#' containing a number of entries corresponding to the number of permutations
#' that have been produced. Each entry is
#' a \code{list} containing:
#' \itemize{
#' \item \code{SITES} a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hyper differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hypo differentially methylated sites between two consecutive generations.The
#' first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hyper differentially methylated sites between the three consecutive
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hypo differentially methylated sites between the three consecutive
#' generations.
#' }
#' }
#' \item \code{TILES} a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hyper differentially methylated positions between two consecutive
#' generations. The first element represents the intersection of the
#' first and second generations; the second element, the intersection of
#' the second and third generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 2 entries,
#' the number of conserved
#' hypo differentially methylated positions between two consecutive
#' generations.The first element represents the intersection of the first and
#' second generations; the second element, the intersection of the second
#' and third generations.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hyper differentially methylated positions between the three consecutive
#' generations.
#' \item \code{HYPO} a \code{list} of \code{integer} with 1 entry,
#' the number of conserved
#' hypo differentially methylated positions between the three consecutive
#' generations.
#' }
#' }
#' }
#' }
#'
#' @seealso
#' \itemize{
#'     \item{\code{\link{extractInfo}}}{ for extracting the
#'     information specific to a subsection of the permutation analysis}
#' }
#'
#' @usage data(methylInheritanceResults)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset containing all results
#' data(methylInheritanceResults)
#'
#' ## Extract information for the intersection between conserved differentially
#' ## methylated sites (type = sites) between the intersection of 2
#' ## generations (inter = i2): F1 and F2 (position = 1)
#' extractInfo(allResults = methylInheritanceResults,
#'     type = "sites", inter="i2", 1)
#'
NULL
