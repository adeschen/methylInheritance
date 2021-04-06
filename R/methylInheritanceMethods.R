#' @title Run all permutations on the specified multi-generational dataset
#'
#' @description Run a permutation analysis, based on Monte Carlo sampling,
#' for testing the hypothesis that the number of conserved differentially
#' methylated elements (sites, tiles or both), between
#' several generations, is associated to an effect inherited from a treatment
#' and that stochastic effect can be dismissed.
#'
#' The multi-generational dataset or the name of the RDS file that contains the
#' dataset can be used as input.
#'
#' The observation analysis can also be run (optional). All permutation
#' results are saved in RDS files.
#'
#' @param methylKitData a \code{list} of \code{methylRawList} entries or the
#' name of the RDS file containing the \code{list}. Each
#' \code{methylRawList} entry must contain all the \code{methylRaw} entries
#' related to one generation (first entry = first generation, second
#' entry = second generation, etc..). The number of generations must
#' correspond to the number
#' of entries in the \code{methylKitData}. At least 2 generations
#' must be present to make a permutation analysis. More information can be
#' found in the methylKit package.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not
#' exist, it will be created. Default: \code{"output"}.
#'
#' @param runObservationAnalysis a \code{logical}, when
#' \code{runObservationAnalysis} = \code{TRUE}, a CpG analysis on the
#' observed dataset is done. Default: \code{TRUE}.
#'
#' @param nbrPermutations, a positive \code{integer}, the total number of
#' permutations that is going to be done. Default: \code{1000}.
#'
#' @param nbrCores a positive \code{integer}, the number of cores to use when
#' processing the analysis. Default: \code{1} and always \code{1} for Windows.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations. The parameter is
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the package
#' \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' corresponds to the \code{lo.count} parameter in the package
#' \code{methylKit}.
#'
#' @param minMethDiff a positive \code{double} between [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter corresponds to the \code{difference} parameter in
#' the methylKit package. Default: \code{10}.
#'
#' @param qvalue a positive \code{double} between [0,1], the cutoff
#' for qvalue of differential methylation statistics. Default: \code{0.01}.
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as an upper cutoff. Bases or regions
#' having higher
#' coverage than this percentile are discarded. The parameter is used for
#' both CpG sites and tiles analysis. The parameter
#' corresponds to the \code{hi.perc} parameter in the package \code{methylKit}.
#' Default: \code{99.9}.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. The parameter is used for both CpG
#' sites and tiles analysis.
#' Default: \code{FALSE}.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}.
#' Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#'
#' @param restartCalculation a \code{logical}, when \code{TRUE}, only
#' permutations that don't have an associated RDS result file are run. Useful
#' to restart a permutation analysis that has been interrupted. Beware that
#' the parameters have to be identical except for this one.
#'
#' @param saveInfoByGeneration a \code{logical}, when \code{TRUE}, the
#' information about differentially methylated sites and tiles for each
#' generation is saved in a RDS file. The information is saved in a different
#' file for each permutation. The files are saved in the directory specified
#' by the \code{outputDir} parameter.
#'
#' @return \code{0}.
#'
#' @seealso \code{\link{mergePermutationAndObservation}} for detail
#' description, in the Value section, of the
#' \code{methylInheritanceAllResults} object as
#' well as its \code{PERMUTATION} section.
#'
#' @examples
#'
#' ## Load methylKit information
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## Run a permutation analysis using the methylKit dataset
#' ## A real analysis would require a much higher number of permutations
#' runPermutation(methylKitData = samplesForTransgenerationalAnalysis,
#'     outputDir = "test_01", runObservationAnalysis = FALSE, type = "sites",
#'     nbrPermutations = 2, vSeed = 221)
#'
#' ## Get results
#' results_01 <- loadAllRDSResults(analysisResultsDir = NULL,
#'     permutationResultsDir = "test_01", doingSites = TRUE,
#'     doingTiles = FALSE)
#'
#' ## Remove results directory
#' if (dir.exists("test_01")) {
#'     unlink("test_01", recursive = TRUE, force = TRUE)
#' }
#'
#' ## Path to a methylKit RDS file
#' methylFile <- system.file("extdata", "methylObj_001.RDS",
#'     package = "methylInheritance")
#'
#' ## Run a permutation analysis using RDS file name
#' ## A real analysis would require a much higher number of permutations
#' runPermutation(methylKitData = methylFile, type = "tiles",
#'     outputDir = "test_02", nbrPermutations = 2, minCovBasesForTiles = 10,
#'     vSeed = 2001)
#'
#' ## Get results
#' results_02 <- loadAllRDSResults(analysisResultsDir = NULL,
#'     permutationResultsDir = "test_02", doingSites = FALSE,
#'     doingTiles = TRUE)
#'
#' ## Remove results directory
#' if (dir.exists("test_02")) {
#'     unlink("test_02", recursive = TRUE, force = TRUE)
#' }
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam bptry bpok
#' @importFrom parallel mclapply nextRNGSubStream
#' @importFrom methods new
#' @export
runPermutation <- function(methylKitData,
                            type=c("both", "sites", "tiles"),
                            outputDir="output",
                            runObservationAnalysis=TRUE,
                            nbrPermutations=1000,
                            nbrCores=1,
                            nbrCoresDiffMeth=1,
                            minReads=10,
                            minMethDiff=10,
                            qvalue=0.01,
                            maxPercReads=99.9,
                            destrand=FALSE,
                            minCovBasesForTiles=0,
                            tileSize=1000,
                            stepSize=1000,
                            vSeed=-1,
                            restartCalculation=FALSE,
                            saveInfoByGeneration=FALSE) {

    # Validate type value
    type <- match.arg(type)

    ## Parameters validation
    validateRunPermutation(methylKitData = methylKitData,
                    type = type, outputDir = outputDir,
                    runObservedAnalysis = runObservationAnalysis,
                    nbrPermutations = nbrPermutations, nbrCores = nbrCores,
                    nbrCoresDiffMeth = nbrCoresDiffMeth, minReads = minReads,
                    minMethDiff = minMethDiff, qvalue = qvalue,
                    maxPercReads = maxPercReads, destrand = destrand,
                    minCovBasesForTiles = minCovBasesForTiles,
                    tileSize = tileSize, stepSize = stepSize, vSeed = vSeed,
                    restartCalculation = restartCalculation,
                    saveInfoByGeneration = saveInfoByGeneration)

    ## Add last slash to path when absent
    if (!is.null(outputDir) &&
            (substr(outputDir, nchar(outputDir), nchar(outputDir)) != "/")) {
        outputDir <- paste0(outputDir, "/")
    }

    ## Load methylKit dataset when needed
    if (is.character(methylKitData)) {
        ## Extract information from RDS file
        methylKitData <- readRDS(methylKitData)
    }

    ## Set vSeed value when negative seed is given
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    RNGkind("L'Ecuyer-CMRG")
    set.seed(vSeed)

    ## Create directory for result files
    if (!is.null(outputDir)) {
        createOutputDir(outputDir,
                        doingSites = any(type %in% c("sites", "both")),
                        doingTiles = any(type %in% c("tiles", "both")),
                        saveInfoByGeneration = saveInfoByGeneration)
    }

    ## Call observation analysis
    if (runObservationAnalysis) {
        runObservation(methylKitData = methylKitData,
                                type = type,
                                outputDir = outputDir,
                                nbrCoresDiffMeth = nbrCoresDiffMeth,
                                minReads = minReads,
                                minMethDiff = minMethDiff,
                                qvalue = qvalue,
                                maxPercReads = maxPercReads,
                                destrand = destrand,
                                minCovBasesForTiles = minCovBasesForTiles,
                                tileSize = tileSize,
                                stepSize = stepSize,
                                vSeed = vSeed,
                                restartCalculation = restartCalculation,
                                saveInfoByGeneration = saveInfoByGeneration)
    }

    ## Upgrade seed
    .Random.seed <- nextRNGSubStream(.Random.seed)

    ## Call permutations in parallel mode
    if (nbrCores > 1) {

        mclapply(seq_len(nbrPermutations), FUN =
                                        runOnePermutationOnAllGenerations,
                            methylInfoForAllGenerations = methylKitData,
                            type = type,
                            outputDir = outputDir,
                            nbrCoresDiffMeth = nbrCoresDiffMeth,
                            minReads = minReads,
                            minMethDiff = minMethDiff,
                            qvalue = qvalue,
                            maxPercReads = maxPercReads,
                            destrand = destrand,
                            minCovBasesForTiles = minCovBasesForTiles,
                            tileSize = tileSize,
                            stepSize = stepSize,
                            restartCalculation = restartCalculation,
                            saveInfoByGeneration = saveInfoByGeneration,
                        mc.cores = nbrCores,
                        mc.preschedule = TRUE)
    } else {
        lapply(seq_len(nbrPermutations), FUN =
                                        runOnePermutationOnAllGenerations,
                                methylInfoForAllGenerations = methylKitData,
                                type = type,
                                outputDir = outputDir,
                                nbrCoresDiffMeth = nbrCoresDiffMeth,
                                minReads = minReads,
                                minMethDiff = minMethDiff,
                                qvalue = qvalue,
                                maxPercReads = maxPercReads,
                                destrand = destrand,
                                minCovBasesForTiles = minCovBasesForTiles,
                                tileSize = tileSize,
                                stepSize = stepSize,
                                restartCalculation = restartCalculation,
                                saveInfoByGeneration = saveInfoByGeneration)
    }

    return(0)
}


#' @title Run a differential methylation analysis on multi-generational
#' dataset
#'
#' @description Run a differential methylation analysis on each generation
#' present in a dataset. The number of conserved differentially
#' methylated elements (sites, tile or both) between generations is
#' them calculated. The
#' methylKit package is used to identify the differentially methylated
#' elements.
#'
#' The multi-generational dataset or the name of the RDS file that contains
#' the dataset can be used as input.
#'
#' The results can also be saved in RDS file (optional).
#'
#' @param methylKitData a \code{list} of \code{methylRawList} entries or the
#' name of the RDS file containing the list. Each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation (first entry = first generation, second entry = second
#' generation, etc..). The number of generations must correspond to the number
#' of entries in the \code{methylKitData}.At least 2 generations
#' must be present to calculate the conserved elements. More information can
#' be found in the methylKit package.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the analysis. If the directory does not
#' exist, it will be created. Default: \code{"output"}.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations.The parameter is
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the package
#' \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the package \code{methylKit}.
#'
#' @param minMethDiff a positive \code{double} between [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter corresponds to the \code{difference} parameter in
#' the methylKit package. Default: \code{10}.
#'
#' @param qvalue a positive \code{double} between [0,1], the cutoff
#' for qvalue of differential methylation statistics. Default: \code{0.01}.
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as an upper cutoff. Bases or regions
#' having higher
#' coverage than this percentile are discarded. The parameter is used for
#' both CpG sites and tiles analysis. The parameter
#' corresponds to the \code{hi.perc} parameter in the package \code{methylKit}.
#' Default: \code{99.9}.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both CpG
#' sites and tiles analysis.
#' Default: \code{FALSE}.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}.
#' Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#'
#' @param restartCalculation a \code{logical}, when \code{TRUE}, only
#' permutations that don't have a RDS result final are run. Useful
#' to restart a permutation analysis that has been interrupted. Beware that
#' the parameters have to be identical except for this one.
#'
#' @param saveInfoByGeneration a \code{logical}, when \code{TRUE}, the
#' information about differentially methylated sites and tiles for each
#' generation is saved in a RDS file. The files are saved in the directory
#' specified by the \code{outputDir} parameter.
#'
#' @return \code{0}.
#'
#' @seealso \code{\link{mergePermutationAndObservation}}  for detail
#' description, in the Value section, of the \code{OBSERVATION} section of the
#' \code{methylInheritanceAllResults} object.
#'
#' @examples
#'
#' ## Load methylation information
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## Run an observation analysis
#' runObservation(methylKitData = samplesForTransgenerationalAnalysis,
#'     outputDir = "test", type = "sites", vSeed = 221)
#'
#' ## Load the results
#' results <- loadAllRDSResults(analysisResultsDir = "test",
#'     permutationResultsDir = NULL, doingSites = TRUE, doingTiles = FALSE)
#'
#' ## Print the results
#' results
#'
#' ## Remove directory
#' if (dir.exists("test")) {
#'     unlink("test", recursive = TRUE, force = FALSE)
#' }
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
runObservation <- function(methylKitData,
                                    type=c("both", "sites", "tiles"),
                                    outputDir="output",
                                    nbrCoresDiffMeth=1,
                                    minReads=10,
                                    minMethDiff=10,
                                    qvalue=0.01,
                                    maxPercReads=99.9,
                                    destrand=FALSE,
                                    minCovBasesForTiles=0,
                                    tileSize=1000,
                                    stepSize=1000,
                                    vSeed=-1,
                                    restartCalculation=FALSE,
                                    saveInfoByGeneration=FALSE) {

    # Validate type value
    type <- match.arg(type)

    ## Parameters validation
    validateRunObservation(methylKitData = methylKitData,
                            type = type, outputDir = outputDir,
                            nbrCoresDiffMeth = nbrCoresDiffMeth,
                            minReads = minReads, minMethDiff = minMethDiff,
                            qvalue = qvalue,
                            maxPercReads = maxPercReads, destrand = destrand,
                            minCovBasesForTiles = minCovBasesForTiles,
                            tileSize = tileSize,
                            stepSize = stepSize, vSeed = vSeed,
                            restartCalculation = restartCalculation,
                            saveInfoByGeneration = saveInfoByGeneration)

    ## Add last slash to path when absent
    if (!is.null(outputDir) &&
        (substr(outputDir, nchar(outputDir), nchar(outputDir)) != "/")) {
        outputDir <- paste0(outputDir, "/")
    }

    ## Load methylKit dataset when needed
    if (is.character(methylKitData)) {
        ## Extract information from RDS file
        methylKitData <- readRDS(methylKitData)
    }

    ## Set vSeed value when negative seed is given
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    RNGkind("L'Ecuyer-CMRG")
    set.seed(vSeed)

    methylInfo <- list(sample = methylKitData, id = 0)

    if (!is.null(outputDir)) {
        doTiles <- any(type %in% c("tiles", "both"))
        doSites <- any(type %in% c("sites", "both"))
        createOutputDir(outputDir, doingSites = doSites, doingTiles = doTiles,
                        saveInfoByGeneration = saveInfoByGeneration)
    }

    ## Extract information
    runOnePermutationOnAllGenerations(id = 0,
                                methylInfoForAllGenerations = methylKitData,
                                type = type, outputDir = outputDir,
                                nbrCoresDiffMeth = nbrCoresDiffMeth,
                                minReads = minReads,
                                minMethDiff = minMethDiff,
                                qvalue = qvalue,
                                maxPercReads = maxPercReads,
                                destrand = destrand,
                                minCovBasesForTiles = minCovBasesForTiles,
                                tileSize = tileSize,
                                stepSize = stepSize,
                                restartCalculation = restartCalculation,
                                saveInfoByGeneration = saveInfoByGeneration)

    return(0)
}


#' @title Load all RDS files created by the permutation  and observation
#' analysis
#'
#' @description  Load all RDS files created by the permutation and
#' observation analysis. The function
#' returns an object of \code{class} "methylInheritanceAllResults" that holds
#' all the pertinent information.
#'
#' @param analysisResultsDir a \code{character} string, the path to the
#' directory that contains the analysis results. The path can be the same as
#' for the \code{permutationResultsDir} parameter. When \code{NULL}, the
#' observation results are not loaded. Default = \code{NULL}.
#'
#' @param permutationResultsDir a \code{character} string, the path to the
#' directory that contains the permutation results. The path can be the same
#' as for the \code{analysisResultsDir} parameter. When \code{NULL}, the
#' permutation results are not loaded. Default = \code{NULL}.
#'
#' @param doingSites a \code{logical}, the data related to differentially
#' methylated sites are loaded when
#' \code{doingSites} = \code{TRUE}. Default: \code{TRUE}.
#'
#' @param doingTiles a \code{logical}, the data related to differentially
#' methylated tiles are loaded when
#' \code{doingTiles} = \code{TRUE}. Default: \code{TRUE}.
#'
#' @param maxID \code{NA} or a positive \code{integer}, the maximum
#' identification number of the permutation files to be loaded. When \code{NA},
#' all files present in the directory are loaded. Default: \code{NA}.
#'
#' @return a \code{list} of class \code{methylInheritanceAllResults}
#' containing the result of the observation analysis as well as the results
#' of all the permutations.
#'
#' @seealso \code{\link{mergePermutationAndObservation}} for detail
#' description, in the Value section, of the
#' \code{methylInheritanceAllResults} object.
#'
#' @examples
#'
#' ## Get the name of the directory where files are stored
#' filesDir <- system.file("extdata", "TEST", package="methylInheritance")
#'
#' ## Load information from files
#' results <- loadAllRDSResults(analysisResultsDir = filesDir,
#'     permutationResultsDir = filesDir, doingSites = TRUE, doingTiles = TRUE)
#'
#' ## Print the observation results
#' results
#'
#' ## Access the results for the first permutation only for sites
#' results$PERMUTATION[[1]]$SITES
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom rebus number_range
#' @export
loadAllRDSResults <- function(analysisResultsDir,
                                    permutationResultsDir,
                                    doingSites=TRUE, doingTiles=FALSE,
                                    maxID = NA) {

    ## Validate parameters
    if (!is.na(maxID) && !is.numeric(maxID) && maxID > 0) {
        stop("maxID must be NA or a positive integer")
    }

    ## Add last slash to analysisResultsDIR when absent
    if (!is.null(analysisResultsDir) &&
        (substr(analysisResultsDir, nchar(analysisResultsDir),
                    nchar(analysisResultsDir)) != "/")) {
        analysisResultsDir <- paste0(analysisResultsDir, "/")
    }

    ## Add last slash to permutationResultsDIR when absent
    if (!is.null(permutationResultsDir) &&
        (substr(permutationResultsDir, nchar(permutationResultsDir),
                    nchar(permutationResultsDir)) != "/")) {
        permutationResultsDir <- paste0(permutationResultsDir, "/")
    }

    result<-list()

    if (!is.na(maxID)) {
        rxRange <- suppressWarnings(number_range(1, maxID))
        newRange <- gsub(pattern = "\\?\\:", replacement = "",
                            toString(rxRange))
        filePattern <- paste0("[^[:digit:]]", newRange,  ".RDS")
    } else {
        filePattern <- "[[:digit:]].RDS"
    }

    ## SITES
    if (doingSites) {
        if (!is.null(analysisResultsDir)) {
            analysisResults <- readRDS(file = paste0(analysisResultsDir,
                                        "SITES/SITES_observed_results.RDS"))
            analysisStruct <- createDataStructure(interGenerationGR =
                                                    analysisResults)
            result[["OBSERVATION"]][["SITES"]] <- analysisStruct
        }

        if (!is.null(permutationResultsDir)) {
            filesInDir <- list.files(path = paste0(permutationResultsDir,
                                                                "SITES/"),
                                pattern = filePattern, all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE,
                                no.. = FALSE)

            sitesPerm <- lapply(filesInDir, FUN = function(x) {
                                                        readRDS(file = x)})

            t <- lapply(sitesPerm, FUN = function(x) {
                    struct <- createDataStructure(interGenerationGR = x)
                    res <- list("SITES" = struct)
                    return(res)})

            result[["PERMUTATION"]] <- t
        }
    }

    ## TILES
    if (doingTiles) {
        if (!is.null(analysisResultsDir)) {
            analysisResults <- readRDS(file = paste0(analysisResultsDir,
                                        "TILES/TILES_observed_results.RDS"))
            analysisStruct <- createDataStructure(interGenerationGR =
                                                    analysisResults)
            result[["OBSERVATION"]][["TILES"]] <- analysisStruct
        }

        if (!is.null(permutationResultsDir)) {
            filesInDir <- list.files(path = paste0(permutationResultsDir,
                                                            "TILES/"),
                                pattern = filePattern, all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE,
                                no.. = FALSE)

            tilesPerm <- lapply(filesInDir, FUN = function(x) {
                                                    readRDS(file = x)})

            t <- lapply(tilesPerm, FUN = function(x) {
                    struct <- createDataStructure(interGenerationGR = x)
                    res <- list("TILES" = struct)
                    return(res)})

            if (!doingSites) {
                result[["PERMUTATION"]] <- t
            } else {
                for (i in seq_len(length(result[["PERMUTATION"]]))) {
                    result[["PERMUTATION"]][[i]]$TILES <- t[[i]]$TILES
                }
            }
        }
    }

    class(result)<-"methylInheritanceAllResults"

    return(result)
}


#' @title Merge the permutation results with the observation results.
#'
#' @description  Merge the permutation results with the observation results.
#' The merging is only needed when permutation and observation have been
#' processed separately. The returned value is a
#' \code{methylInheritanceAllResults} object that can be used by
#' the \code{extractInfo} function.
#'
#' @param permutationResults a \code{list} with 1 entry called
#' \code{PERMUTATION}. The  \code{PERMUTATION} entry is a \code{list} with
#' a number of entries corresponding
#' to the number of permutations that have been processed. Each entry contains
#' the result of one permutation.
#'
#' @param observationResults a \code{list} with 1 entry called
#' \code{OBSERVATION}. The \code{OBSERVATION} entry is a \code{list} containing
#' the result obtained
#' with the observed dataset (not shuffled).
#'
#' @return a \code{list} of class \code{methylInheritanceAllResults} with
#' 2 entries. The 2 entries are:
#' \itemize{
#' \item \code{PERMUTATION} \code{list} with a number of entries corresponding
#' to the number of permutations that have been processed. Each entry contains
#' the result of one permutation.The elements in each entry are:
#' \itemize{
#' \item \code{SITES} Only present when a sites analysis has been achieved,
#' a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between two consecutive generations.The
#' first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc.The number of entries depends on the number
#' of generations.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc. The number of entries depends on the number of
#' generations.
#' }
#' }
#' \item \code{TILES} Only present when a tiles analysis has been achieved,
#' a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated positions between two consecutive
#' generations. The first element represents the intersection of the
#' first and second generations; the second element, the intersection of
#' the second and third generations; etc.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated positions between two consecutive
#' generations.The first element represents the intersection of the first and
#' second generations; the second element, the intersection of the second
#' and third generations; etc.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated positions between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc. The number of entries depends on the number
#' of generations.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated positions between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc.The number of entries depends on the number of
#' generations.
#' }
#' }
#' }
#' \item \code{OBSERVATION} a \code{list} containing
#' the result obtained with the observed dataset (not shuffled). The
#' elements are:
#' \itemize{
#' \item \code{SITES} Only present when a sites analysis has been achieved,
#' a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between two consecutive generations.The
#' first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc.The number of entries depends on the number
#' of generations.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc. The number of entries depends on the number of
#' generations.
#' }
#' }
#' \item \code{TILES} Only present when a tiles analysis has been achieved,
#' a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated positions between two consecutive
#' generations. The first element represents the intersection of the
#' first and second generations; the second element, the intersection of
#' the second and third generations; etc.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated positions between two consecutive
#' generations.The first element represents the intersection of the first and
#' second generations; the second element, the intersection of the second
#' and third generations; etc.
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated positions between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc. The number of entries depends on the number
#' of generations.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated positions between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc.The number of entries depends on the number of
#' generations.
#' }
#' }
#' }
#' }
#'
#' @examples
#'
#' ## Create a observation result
#' observed <- list()
#' observed[["OBSERVATION"]] <- list()
#' observed[["OBSERVATION"]][["SITES"]] <- list()
#' observed[["OBSERVATION"]][["SITES"]][["i2"]] <- list(HYPER = list(11, 10),
#'     HYPO = list(13, 12))
#' observed[["OBSERVATION"]][["SITES"]][["iAll"]] <- list(HYPER = list(1),
#'     HYPO = list(3))
#'
#' ## Create a permutation result containing only 1 permutation result
#' ## Real perumtations results would have more entries
#' permutated <- list()
#' permutated[["PERMUTATION"]] <- list()
#' permutated[["PERMUTATION"]][[1]] <- list()
#' permutated[["PERMUTATION"]][[1]][["SITES"]] <- list()
#' permutated[["PERMUTATION"]][[1]][["SITES"]][["i2"]] <- list(HYPER =
#'     list(11, 12), HYPO = list(8, 11))
#' permutated[["PERMUTATION"]][[1]][["SITES"]][["iAll"]] <- list(HYPER =
#'     list(0), HYPO = list(1))
#'
#' ## Merge permutation and observation results
#' mergePermutationAndObservation(permutationResults = permutated,
#'     observationResults = observed)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
mergePermutationAndObservation <- function(permutationResults,
                                                observationResults) {

    ## Validate parameters
    validateMergePermutationAndObservation(permutationResults,
                                            observationResults)

    mergedData <- list()
    mergedData[["PERMUTATION"]] <- permutationResults[["PERMUTATION"]]
    mergedData[["OBSERVATION"]] <- observationResults[["OBSERVATION"]]

    class(mergedData)<-"methylInheritanceAllResults"

    return(mergedData)
}


#' @title Extract the information specific to a subsection of the permutation
#' analysis
#'
#' @description  Extract the information specific to a subsection of the
#' permutation analysis. The extracted information will be specific to one
#' type of differential methylation analysis (tiles or sites), to one type
#' of intersection (two consecutive generation or more) and to one specific
#' group of generations.
#'
#' @param allResults a \code{list} of class \code{methylInheritanceAllResults}
#' as created by the
#' \code{runPermutation} function. The \code{list} must contain
#' two entries : \code{"PERMUTATION"} and \code{"OBSERVATION"}. The
#' \code{"PERMUTATION"} \code{list} must contain all results from all
#' permutations while
#' the \code{"OBSERVATION"} \code{list} must contain the result obtained with
#' the observed dataset (not shuffled).
#'
#' @param type One of the \code{"sites"} or \code{"tiles"} strings.
#' Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases \code{type} = \code{"sites"}; for
#' differentially methylated regions \code{type} = \code{"tiles"}.
#' Default: \code{"sites"}.
#'
#' @param inter One of the \code{"i2"} or \code{"iAll"} strings. Specifies the
#' type of intersection should be returned. For
#' retrieving intersection results between two consecutive generations
#' \code{inter} = \code{"i2"}; for intersection results between three
#' generations or more \code{inter} = \code{"iAll"}.
#' Default: \code{"i2"}.
#'
#' @param position a positive \code{integer}, the position in the \code{list}
#' where the information will be extracted. Default=\code{1}.
#'
#' @return a \code{data.frame}
#' containing the observation results (using real
#' data) and the permutation results (using shuffled data). Both hyper and
#' hypo differentially conserved methylation results are present.
#'
#' @examples
#'
#' ## Get the name of the directory where files are stored
#' filesDir <- system.file("extdata", "TEST", package="methylInheritance")
#'
#' ## Load information from files
#' results <- loadAllRDSResults(analysisResultsDir = filesDir,
#'     permutationResultsDir = filesDir, doingSites = TRUE, doingTiles = TRUE)
#'
#' ## Extract information for the intersection between conserved differentially
#' ## methylated sites (type = sites) between the intersection of 2
#' ## generations (inter = i2): F1 and F2 (position = 1)
#' info <- extractInfo(allResults = results, type = "sites", inter="i2", 1)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
extractInfo <- function(allResults, type=c("sites", "tiles"),
                            inter=c("i2", "iAll"), position=1) {

    ## Validate type value
    type <- match.arg(type)

    ## Validate type value
    inter <- match.arg(inter)

    ## Validate other parameters
    validateExtractInfo(allResults = allResults, type, inter, position)

    type <- toupper(type)

    real <- allResults[["OBSERVATION"]][[type]][[inter]]

    dataConserved <- data.frame(TYPE=c("HYPO", "HYPER"),
                                RESULT=c(real[["HYPO"]][[position]],
                                            real[["HYPER"]][[position]]),
                                SOURCE=c("OBSERVATION", "OBSERVATION"))

    for (i in seq_len(length(allResults[["PERMUTATION"]]))) {
        permutation <- allResults[["PERMUTATION"]][[i]][[type]][[inter]]
        dataConserved <- rbind(dataConserved,
                                data.frame(TYPE=c("HYPO", "HYPER"),
                        RESULT=c(permutation[["HYPO"]][[position]],
                                    permutation[["HYPER"]][[position]]),
                        SOURCE=c("PERMUTATION", "PERMUTATION")))
    }

    return(dataConserved)
}


#' @title Generate a graph for a permutation analysis
#'
#' @description  Generate a graph for a permutation analysis using observed
#' and shuffled results.
#'
#' @param formatForGraphDataFrame a \code{data.frame} containing the
#' observation results (using real
#' data) and the permutation results (using shuffled data). Both hyper and
#' hypo differentially conserved methylation results must be present. The
#' \code{data.frame} must have 3 columns : "TYPE", "RESULT" and "SOURCE".
#' The "TYPE" can be either "HYPER" or "HYPO". The "RESULT" is the number
#' of conserved differentially elements. The "SOURCE" can be either
#' "OBSERVATION" or "PERMUTATION".
#'
#' @return a graph showing the permutation analysis results
#'
#' @examples
#'
#' ## Loading dataset containing all results
#' data(methylInheritanceResults)
#'
#' ## Extract information for the intersection between conserved differentially
#' ## methylated sites (type = sites) between the intersection of 2
#' ## generations (inter = i2): F2 and F3 (position = 2)
#' info <- extractInfo(allResults = methylInheritanceResults,
#'     type = "sites", inter="i2", 2)
#'
#' ## Create graph
#' plotGraph(info)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom ggplot2 ggplot geom_text facet_grid theme geom_vline
#' geom_histogram labs aes scale_color_manual
#' @importFrom gridExtra grid.arrange tableGrob
#' @export
plotGraph <- function(formatForGraphDataFrame) {

    # Basic graph using data.frame
    # Columns names : TYPE (HYPER or HYPO), RESULT (nbr conseved sites),
    # SOURCE (OBSERVED or PERMUTATION)
    p <- ggplot(data=formatForGraphDataFrame,
                    aes(x=formatForGraphDataFrame$RESULT)) +
                    geom_histogram(col="blue", fill="lightblue",
                                    binwidth=2, alpha = .2) +
        labs(title = "") +
        labs(x = "Number of conserved differentially methylated sites",
                y = "Frequency")

    # Split to have one section for HYPER and one for HYPO
    p <- p + facet_grid(.~TYPE)

    # Add vertical line corresponding to the number of conserved elements
    # in the observed results (real results)

    RESULT <- NULL
    interceptFrame <- subset(formatForGraphDataFrame,
                            formatForGraphDataFrame$SOURCE == "OBSERVATION")
    p <- p + geom_vline(data = interceptFrame,
                        aes(xintercept=RESULT, color="observed"),
                        linetype="longdash", show.legend=TRUE)

    p <- p + scale_color_manual(name = "", values = c(observed = "red")) +
        theme(legend.position="bottom")

    # Calculate the significant level for HYPER AND HYPO
    signif <- calculateSignificantLevel(formatForGraphDataFrame)

    signifLevelHypo  <- signif[["HYPO"]]
    signifLevelHyper <- signif[["HYPER"]]

    hyperNumber <- interceptFrame[interceptFrame$TYPE == "HYPER",]$RESULT
    hypoNumber  <- interceptFrame[interceptFrame$TYPE == "HYPO",]$RESULT

    # Number of observed conserved elements as annotated text
    info <- data.frame(type = c("HYPER", "HYPO"),
                            lab = c(hyperNumber, hypoNumber),
                            signif = c(signifLevelHyper, signifLevelHypo))
    colnames(info)<-c("Type", "Observed Value", "Significant Level")

    ## Put graph and table in grid
    g <- grid.arrange(p, tableGrob(info, rows = NULL), nrow = 2,
                        heights = c(2, 1), clip = FALSE)

    return(g)
}


#' @title Load convergence information from RDS files
#'
#' @description  Load convergence information from RDS files.
#'
#' @param analysisResultsDir a \code{character} string, the path to the
#' directory that contains the analysis results. The path can be the same as
#' for the \code{permutatioNResultsDir} parameter.
#'
#' @param permutationResultsDir a \code{character} string, the path to the
#' directory that contains the permutation results. The path can be the same
#' as for the \code{analysisResultsDir} parameter.
#'
#' @param type One of the \code{"sites"} or \code{"tiles"} strings.
#' Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases \code{type} = \code{"sites"}; for
#' differentially methylated regions \code{type} = \code{"tiles"}.
#' Default: \code{"sites"}.
#'
#' @param inter One of the \code{"i2"} or \code{"iAll"} strings. Specifies the
#' type of intersection should be returned. For
#' retrieving intersection results between two consecutive generations
#' \code{inter} = \code{"i2"}; for intersection results between three
#' generations or more \code{inter} = \code{"iAll"}.
#' Default: \code{"i2"}.
#'
#' @param position a positive \code{integer}, the position in the \code{list}
#' where the information will be extracted.
#'
#' @param by a \code{integer}, the increment of the number of permutations
#' where the significant level is tested. Default: 100.
#'
#' @return a graph showing the evolution of the significant level with the
#' number of permutations
#'
#' @examples
#'
#' ## Get the name of the directory where files are stored
#' filesDir <- system.file("extdata", "TEST", package="methylInheritance")
#'
#' ## Load convergence information
#' results <- loadConvergenceData(analysisResultsDir = filesDir,
#'     permutationResultsDir = filesDir, type="sites", inter="i2", position=1,
#'     by=1)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
loadConvergenceData <- function(analysisResultsDir,
                            permutationResultsDir, type = c("sites", "tiles"),
                            inter = c("i2", "iAll"), position, by = 100) {
    ## Validate type value
    type <- match.arg(type)

    ## Validate type value
    inter <- match.arg(inter)

    ## Validate other parameters
    validateLoadConvergenceData(analysisResultsDir, permutationResultsDir,
                                    position, by)

    ## Add last slash to analysisResultsDir when absent
    if (substr(analysisResultsDir, nchar(analysisResultsDir),
                nchar(analysisResultsDir)) != "/") {
        analysisResultsDir <- paste0(analysisResultsDir, "/")
    }

    ## Add last slash to permutationResultsDir when absent
    if (substr(permutationResultsDir, nchar(permutationResultsDir),
                nchar(permutationResultsDir)) != "/") {
        permutationResultsDir <- paste0(permutationResultsDir, "/")
    }

    if (type == "sites") {
        doingSites = TRUE
        doingTiles = FALSE
    } else {
        doingSites = FALSE
        doingTiles = TRUE
    }

    ## Get maximum number of files and create a incrementing sequence using
    ## the "by" parameter
    filesInDir <- list.files(path = paste0(permutationResultsDir, toupper(type),
                        "/"), pattern = "[[:digit:]].RDS", all.files = FALSE,
                        full.names = TRUE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

    nbFiles <- length(filesInDir)
    seqFiles <- unique(c(seq(by, nbFiles, by), nbFiles))

    ## The returned data.frame
    final <- data.frame(NBR_PERMUTATIONS = integer(), ELEMENT = character(),
                            TYPE = character(), ANALYSIS = character(),
                            POSITION = integer(),
                            SIGNIFICANT_LEVEL = numeric(),
                            stringsAsFactors=FALSE)

    for (i in seqFiles) {
        ## Load data associated to the selected number of files
        data <- loadAllRDSResults(analysisResultsDir = analysisResultsDir,
                    permutationResultsDir = permutationResultsDir,
                    doingSites = doingSites, doingTiles = doingTiles,
                    maxID = i)

        ## Extract info
        info <- extractInfo(allResults = data, type = type,
                            inter = inter, position)

        ## Calculate significant level
        result <- calculateSignificantLevel(info)

        ## Create data.frame using extracted information
        temp <- data.frame(NBR_PERMUTATIONS = rep(i, 2),
                    ELEMENT = rep(toupper(type), 2), ANALYSIS = rep(inter, 2),
                    POSITION = rep(position, 2), TYPE = c("HYPER", "HYPO"),
                    SIGNIFICANT_LEVEL = c(result$HYPER, result$HYPO),
                    stringsAsFactors = FALSE)

        ## Add the information to the returned data.frame
        final <- rbind(final, temp)
    }

    ## The final data.frame containing significant levels for all entries
    ## of the incrementing sequence of number of permutations
    return(final)
}


#' @title Generate a graph showing the convergence for a permutation analysis
#'
#' @description  Generate a graph showing the convergence for a permutation
#' analysis using observed and permuted results.
#'
#' @param dataFrameConvergence a \code{data.frame} containing the
#' significant levels at different number of cycles (total number of
#' permuted data analysed).  The
#' \code{data.frame} must have 6 columns : "NBR_PERMUTATIONS", "ELEMENT".
#' "ANALYSIS", "POSITION", "TYPE" and "SIGNIFICANT_LEVEL". The "ELEMENT" can
#' be either "SITES" or "TILES". The "TYPE" can be either "HYPER" or "HYPO".
#'
#' @return a \code{ggplot} object.
#'
#' @examples
#'
#' ## Get the name of the directory where files are stored
#' filesDir <- system.file("extdata", "TEST", package="methylInheritance")
#'
#' ## Extract convergenc information for F1 and F2 and F3
#' data <- loadConvergenceData(analysisResultsDir = filesDir,
#'     permutationResultsDir = filesDir, type = "sites", inter = "iAll",
#'     position = 1, by = 1)
#'
#' ## Create convergence graph
#' plotConvergenceGraph(data)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom ggplot2 ggplot geom_point geom_line facet_grid aes xlab ylab
#' @export
plotConvergenceGraph <- function(dataFrameConvergence) {

    ## Create graph
    NBR_PERMUTATIONS <- NULL
    SIGNIFICANT_LEVEL <- NULL
    graph <- ggplot(data=dataFrameConvergence, aes(x=NBR_PERMUTATIONS,
        y=SIGNIFICANT_LEVEL)) + geom_point(color = 'blue', size = 2) +
        geom_line(color='blue', linetype = "dashed", size = 1) +
        facet_grid(TYPE ~ .) + ylab("Significant Level") +
        xlab("Number of permutations")

    ## Return graph
    return(graph)
}
