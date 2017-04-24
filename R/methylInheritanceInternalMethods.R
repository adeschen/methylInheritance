#' @title Parameters validation for the \code{\link{runPermutation}} function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{runPermutation}} function.
#'
#' @param methylKitData a \code{list} of \code{methylRawList} entries or the
#' name of the RDS file containing the \code{list}. Each
#' \code{methylRawList} entry must contain all the \code{methylRaw} entries
#' related to one generation (first entry = first generation, second
#' entry = second generation, etc..). The number of generations must
#' correspond to the number
#' of entries in the \code{methylKitData}. At least 2 generations
#' must be present to do a permutation analysis. More information can be found
#' in the methylKit package.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param runObservedAnalysis a \code{logical}, when \code{runObservedAnalysis}
#' = \code{TRUE}, a CpG analysis on the observed dataset is done.
#'
#' @param nbrPermutations, a positive \code{integer}, the total number of
#' permutations that is going to be done.
#'
#' @param nbrCores a positive \code{integer}, the number of cores to use when
#' processing the analysis.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations.Parameter
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in
#' the \code{methylKit} package.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' corresponds to the \code{lo.count} parameter in the  \code{methylKit}
#' package.
#'
#' @param minMethDiff a positive \code{double} betwwen [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter corresponds to the \code{difference} parameter in
#' the \code{methylKit} package.
#'
#' @param qvalue a positive \code{double} betwwen [0,1], the cutoff
#' for qvalue of differential methylation statistic. TODO
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both CpG
#' sites and tiles analysis.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}. Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used.
#'
#' @param restartCalculation a \code{logical}, when \code{TRUE}, only
#' permutations that don't have an associated RDS result file are run. Useful
#' to restart a permutation analysis that has been interrupted.
#'
#' @param saveInfoByGeneration a \code{logical}, when \code{TRUE}, the
#' information about differentially methylated sites and tiles for each
#' generation is saved in a RDS file. The information is saved in a different
#' file for each permutation. The files are only saved when the
#' \code{outputDir} is not \code{NULL}.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load dataset
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## The function returns 0 when all paramaters are valid
#' methylInheritance:::validateRunPermutation(
#'     methylKitData = samplesForTransgenerationalAnalysis, type = "sites",
#'     outputDir = NULL, runObservedAnalysis = TRUE,
#'     nbrPermutations = 10000, nbrCores = 1,
#'     nbrCoresDiffMeth = 1, minReads = 10, minMethDiff = 25, qvalue = 0.01,
#'     maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
#'     tileSize = 1000, stepSize = 500, vSeed = 12, restartCalculation = FALSE,
#'     saveInfoByGeneration = FALSE)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{methylInheritance:::validateRunPermutation(
#'     methylKitData = "HI",type = "tiles", outputDir = NULL,
#'     runObservedAnalysis = FALSE, nbrPermutations = 10000, nbrCores = 1,
#'     nbrCoresDiffMeth = 1, minReads = 10, minMethDiff = 25, qvalue = 0.01,
#'     maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
#'     tileSize = 1000, stepSize = 500, vSeed = 12, restartCalculation = FALSE,
#'     saveInfoByGeneration = FALSE)}
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunPermutation <- function(methylKitData,
                            type, outputDir, runObservedAnalysis,
                            nbrPermutations, nbrCores, nbrCoresDiffMeth,
                            minReads, minMethDiff, qvalue, maxPercReads,
                            destrand, minCovBasesForTiles, tileSize,
                            stepSize, vSeed, restartCalculation,
                            saveInfoByGeneration) {

    ## Validate methylKitData, outputDir, nbrCoresDiffMeth
    ## minReads, minMethDiff, qvalue, maxPercReads, destrand,
    ## minCovBasesForTiles, tileSize, stepSize, vSeed
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

    ## Validate that the runObservedAnalysis is a logical
    if (!is.logical(runObservedAnalysis)) {
        stop("runObservedAnalysis must be a logical")
    }

    ## Validate that nbrCores is an positive integer
    if (!(isSingleInteger(nbrCores) || isSingleNumber(nbrCores)) ||
            as.integer(nbrCores) < 1) {
        stop("nbrCores must be a positive integer or numeric")
    }

    ## Validate that nbrCores is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" && as.integer(nbrCores) != 1) {
        stop("nbrCores must be 1 on a Windows system.")
    }

    ## Validate that nbrPermutations is an positive integer
    if (!(isSingleInteger(nbrPermutations) ||
            isSingleNumber(nbrPermutations)) ||
            as.integer(nbrPermutations) < 1) {
        stop("nbrPermutations must be a positive integer or numeric")
    }

    return(0)
}


#' @title Validation of some parameters of the
#' \code{\link{runObservation}} function
#'
#' @description Validation of some parameters needed by the public
#' \code{\link{runObservation}} function.
#'
#' @param methylKitData a \code{list} of \code{methylRawList} entries or the
#' name of the RDS file containing the list. Each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation (first entry = first generation, second entry = second
#' generation, etc..). The number of generations must correspond to the number
#' of entries in the \code{methylKitData}. At least 2 generations
#' must be present to calculate the conserved elements. More information can
#' be found in the methylKit package.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations.Parameter
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in
#' the \code{methylKit} package.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the  \code{methylKit}
#' package.
#'
#' @param minMethDiff a positive \code{double} betwwen [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter correspond to the \code{difference} parameter in
#' the \code{methylKit} package.
#'
#' @param qvalue a positive \code{double} betwwen [0,1], the cutoff
#' for qvalue of differential methylation statistic.
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both CpG
#' sites and tiles analysis.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}. Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used.
#'
#' @param restartCalculation a \code{logical}, when \code{TRUE}, only
#' permutations that don't have an associated RDS result file are run. Useful
#' to restart a permutation analysis that has been interrupted.
#'
#' @param saveInfoByGeneration a \code{logical}, when \code{TRUE}, the
#' information about differentially methylated sites and tiles for each
#' generation is saved in a RDS file. The information is saved in a different
#' file for each permutation. The files are only saved when the
#' \code{outputDir} is not \code{NULL}.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load dataset
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## The function returns 0 when all paramaters are valid
#' methylInheritance:::validateRunObservation(
#'     methylKitData = samplesForTransgenerationalAnalysis, type = "sites",
#'     outputDir = NULL, nbrCoresDiffMeth = 1, minReads = 10,
#'     minMethDiff = 25, qvalue = 0.01,
#'     maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
#'     tileSize = 1000, stepSize = 500, vSeed = 12, restartCalculation = TRUE,
#'     saveInfoByGeneration = FALSE)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{methylInheritance:::validateRunObservation(
#'     methylKitData = samplesForTransgenerationalAnalysis,
#'     type = "tiles", outputDir = NULL, nbrCoresDiffMeth = 1, minReads = "HI",
#'     minMethDiff = 25, qvalue = 0.01,
#'     maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
#'     tileSize = 1000, stepSize = 500, vSeed = 12, restartCalculation = FALSE,
#'     saveInfoByGeneration = FALSE)}
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunObservation <- function(methylKitData,
                                    type, outputDir,
                                    nbrCoresDiffMeth,
                                    minReads, minMethDiff, qvalue,
                                    maxPercReads, destrand,
                                    minCovBasesForTiles, tileSize,
                                    stepSize, vSeed, restartCalculation,
                                    saveInfoByGeneration) {

    ## Validate that methylKitData is a valid RDS file when string is passed
    if (is.character(methylKitData)) {
        if (!file.exists(methylKitData)) {
            stop(paste0("The file \"", methylKitData, "\" does not exist."))
        } else {
            methylKitData <- readRDS(methylKitData)
        }
    }

    ## Validate that methylKitData is a list of methylRawList
    if (class(methylKitData) != "list" ||
            !all(sapply(methylKitData, class) == "methylRawList")) {
        stop(paste0("methylKitData must be a list containing ",
                    "\"methylRawList\" entries; each entry must contain ",
                    "all \"methylRaw\" objects related to one generation"))
    }

    ## Validate that the output_dir is an not empty string
    if (!is.null(outputDir) && !is.character(outputDir)) {
        stop("output_dir must be a character string or NULL")
    }

    ## Validate that nbrCoresDiffMeth is an positive integer
    if (!(isSingleInteger(nbrCoresDiffMeth) ||
            isSingleNumber(nbrCoresDiffMeth)) ||
        as.integer(nbrCoresDiffMeth) < 1) {
        stop("nbrCoresDiffMeth must be a positive integer or numeric")
    }

    ## Validate that nbrCoresDiffMeth is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" &&
            as.integer(nbrCoresDiffMeth) != 1) {
        stop("nbrCoresDiffMeth must be 1 on a Windows system.")
    }

    ## Validate that minReads is an positive integer
    if (!(isSingleInteger(minReads) || isSingleNumber(minReads)) ||
        as.integer(minReads) < 1) {
        stop("minReads must be a positive integer or numeric")
    }

    ## Validate that minMethDiff is an positive double between [0,100]
    if (!(isSingleNumber(minMethDiff)) ||
            minMethDiff < 0.00 || minMethDiff > 100.00) {
        stop("minMethDiff must be a positive double between [0,100]")
    }

    ## Validate that qvalue is an positive double between [0,1]
    if (!(isSingleNumber(qvalue)) ||
            qvalue < 0.00 || qvalue > 1.00) {
        stop("qvalue must be a positive double between [0,1]")
    }

    ## Validate that maxPercReads is an positive double between [0,100]
    if (!(isSingleNumber(maxPercReads)) ||
            maxPercReads < 0.00 || maxPercReads > 100.00) {
        stop("maxPercReads must be a positive double between [0,100]")
    }

    ## Validate that destrand is a logical
    if (!is.logical(destrand)) {
        stop("destrand must be a logical")
    }

    if (any(type %in% c("both", "tiles"))) {
        ## Validate that minCovBasesForTiles is an positive integer
        if (!(isSingleInteger(minCovBasesForTiles) ||
                isSingleNumber(minCovBasesForTiles)) ||
                    as.integer(minCovBasesForTiles) < 0) {
            stop("minCovBasesForTiles must be a positive integer or numeric")
        }

        ## Validate that tileSize is an positive integer
        if (!(isSingleInteger(tileSize) || isSingleNumber(tileSize)) ||
                as.integer(tileSize) < 1) {
            stop("tileSize must be a positive integer or numeric")
        }

        ## Validate that stepSize is an positive integer
        if (!(isSingleInteger(stepSize) || isSingleNumber(stepSize)) ||
                as.integer(stepSize) < 1) {
            stop("stepSize must be a positive integer or numeric")
        }
    }

    ## Validate that vSeed is an integer
    if (!(isSingleInteger(vSeed) || isSingleNumber(vSeed))) {
        stop("vSeed must be an integer or numeric")
    }
    ## Validate that restartCalculation is a logical
    if (!is.logical(restartCalculation)) {
        stop("restartCalculation must be a logical")
    }

    ## Validate that saveInfoByGeneration is a logical
    if (!is.logical(saveInfoByGeneration)) {
        stop("saveInfoByGeneration must be a logical")
    }

    return(0)
}


#' @title Validation of some parameters of the
#' \code{\link{extractInfo}} function
#'
#' @description Validation of some parameters needed by the public
#' \code{\link{extractInfo}} function.
#'
#' @param allResults a \code{list} as created by the
#' \code{runPermutation} or the \code{loadAllRDSResults} functions.
#'
#' @param type One of the \code{"sites"} or \code{"tiles"} strings.
#' Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases \code{type} = \code{"sites"}; for
#' differentially methylated regions \code{type} = \code{"tiles"}.
#'
#' @param inter One of the \code{"i2"} or \code{"iAll"} strings. Specifies the
#' type of intersection should be returned. For
#' retrieving intersection results between two consecutive generations
#' \code{inter} = \code{"i2"}; for intersection results between three
#' generations or more \code{inter} = \code{"iAll"}.
#'
#' @param position a positive \code{integer}, the position in the \code{list}
#' where the information will be extracted. The position must be an existing
#' position inside \code{allResults}
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load dataset
#' data(methylInheritanceResults)
#'
#' ## The function returns 0 when all paramaters are valid
#' methylInheritance:::validateExtractInfo(
#'     allResults = methylInheritanceResults, type = "sites",
#'     inter = "i2", 2)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{methylInheritance:::validateExtractInfo(
#'     allResults = methylInheritanceResults, type = "sites",
#'     inter = "i2", 12)}
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateExtractInfo <- function(allResults, type, inter, position) {

    if (position < 1) {
        stop("position must be a positive integer")
    }

    if (!"methylInheritanceAllResults" %in% class(allResults)) {
        stop("allResults must be of class \"methylInheritanceAllResults\"")
    }

    if (!is.list(allResults)) {
        stop("allResults must be a list")
    }

    if (is.null(allResults$OBSERVATION)) {
        stop("allResults must have an element called \"OBSERVATION\"")
    }

    if (is.null(allResults$OBSERVATION[[toupper(type)]])) {
        stop("allResults must have an element called \"", toupper(type),
                "\" in its \"OBSERVATION\" list")
    }

    if (is.null(allResults$OBSERVATION[[toupper(type)]][[inter]])) {
        stop("allResults must have an element called \"", inter,
                "\" in the \"", toupper(type),
                "\" list present in its \"OBSERVATION\" list")
    }

    if (position > length(allResults$OBSERVATION[[toupper(type)]][[inter]])) {
        stop(paste0("position must correspond to a valid entry in the \"",
            "allResults$OBSERVATION[[", toupper(type), "]][[", inter, "]]"))
    }

    if (is.null(allResults$PERMUTATION)) {
        stop("allResults must have an element called \"PERMUTATION\"")

    }

    if (!is.list(allResults$PERMUTATION)) {
        stop(paste0("allResults must have an element called \"PERMUTATION\". ",
                "The \"PERMUTATION\" must be a list"))
    }

    if (length(allResults$PERMUTATION) < 1) {
        stop(paste0("allResults must have an element called \"PERMUTATION\". ",
                "The \"PERMUTATION\" must be a list with at least one entry"))
    }

    ## Validate that all entries in allResults$PERMUTATION must contain
    ## a entry corresponding to "type"
    results <- sapply(allResults$PERMUTATION, function(x) {
                is.null(x[[toupper(type)]])})
    if (any(results)) {
        stop(paste0("all entries in allResults$PERMUTATION must have an ",
                "element called \"", toupper(type)))
    }

    return(0)
}


#' @title Validation of some parameters of the
#' \code{\link{mergePermutationAndObservation}} function
#'
#' @description Validation of some parameters needed by the public
#' \code{\link{mergePermutationAndObservation}} function.
#'
#' @param @param permutationResults a \code{list} with 1 entry called
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
#' @return \code{0} indicating that all parameters validations have been
#' successful.
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
#' methylInheritance:::validateMergePermutationAndObservation(
#'     permutationResults = permutated, observationResults = observed)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{methylInheritance:::validateMergePermutationAndObservation(
#'     permutationResults = permutated, observationResults = NULL)}
#'
#' @author Astrid Deschenes
#' @keywords internal
validateMergePermutationAndObservation <- function(permutationResults,
                                                    observationResults) {

    if (!is.list(permutationResults)) {
        stop("permutationResults must be a list")
    }

    if (is.null(permutationResults$PERMUTATION)) {
        stop("permutationResults must have an element called \"PERMUTATION\"")
    }

    if (!is.list(observationResults)) {
        stop("observationResults must be a list")
    }

    if (is.null(observationResults$OBSERVATION)) {
        stop("observationResults must have an element called \"OBSERVATION\"")
    }

    return(0)
}


#' @title Transform results from a CpG site or region analysis done on mutliple
#' generations into a \code{list} of \code{GRanges} objects
#'
#' @description Transform a \code{list} of \code{methylDiff} objects into
#' a \code{list} of \code{GRanges} objects. Each \code{methylDiff} object
#' represent a CpG site or region analysis done on one generation.
#'
#' @param methDiff a \code{list} of S4 \code{methylDiff} class objects, each
#' entry of the \code{list} represents the differentially methylated results
#' for one generation (first entry = first genertation, second entry =
#' second generation, etc..). Each \code{methylDiff} object holds statistics
#' and locations
#' for differentially methylated regions/bases.
#'
#' @param pDiff a positive \code{double} between \code{0} and \code{100},
#' the cutoff for absolute value of methylation percentage change
#' between test and control.
#'
#' @param qvalue a positive \code{double} inferior to \code{1}, the cutoff
#' for qvalue of differential methylation statistic.
#'
#' @param type One of the \code{"hyper"},\code{"hypo"} or \code{"all"} strings,
#' the string specifies what type of differentially methylated bases/tiles
#' should be treated  For
#' retrieving hyper-methylated tiles/sites \code{type} = \code{"hyper"}; for
#' hypo-methylated \code{type} = \code{"hypo"}. Default: \code{"all"}.
#'
#' @return a \code{list} of \code{GRanges} objects, each
#' entry of the \code{list} represents the differentially methylated results
#' for one generation (first entry = first genertation, second entry =
#' second generation, etc..). Each \code{GRanges} object holds statistics
#' for differentially methylated regions/bases.
#'
#' @examples
#'
#' ## Load permutation results on sites
#' permutationResultsFile <- system.file("extdata",
#'     "permutationResultsForSites.RDS", package="methylInheritance")
#' permutationResults <- readRDS(permutationResultsFile)
#'
#' ## Transform result to GRanges
#' resultsGR <- methylInheritance:::getGRangesFromMethylDiff(methDiff =
#'     permutationResults, pDiff = 10, qvalue = 0.01, type = "hyper")
#'
#' @author Pascal Belleau
#' @importFrom methylKit getMethylDiff
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @keywords internal
getGRangesFromMethylDiff <- function(methDiff, pDiff, qvalue,
                                        type = c("all", "hyper", "hypo")) {
    # Validate type value
    type <- match.arg(type)

    ## Transform each methylDiff object present in the list to a
    ## GRanges object
    methDiffK <- lapply(1:length(methDiff), FUN = function(i, methDiff,
                                                    pDiff, qCut, typeD) {
        methK <- getMethylDiff(methDiff[[i]], difference = pDiff,
                                qvalue = qCut, type = typeD)
        GRanges(seqnames = methK$chr, ranges = IRanges(start = methK$start,
                                                        end = methK$end),
                strand = methK$strand, pvalue = methK$pvalue,
                qvalue = methK$qvalue, meth.diff = methK$meth.diff)
    }, methDiff = methDiff, pDiff = pDiff, qCut = qvalue, typeD = type)

    return(methDiffK)
}


#' @title Calculate the intersection of the differentially methylated
#' results for two
#' or more consercutive generations
#'
#' @description Calculate the intersection of the differentially methylated
#' results for two
#' or more consercutive generations using a \code{list} of \code{GRanges} where
#' each entry represents the results for one generation.
#'
#' @param resultAllGenGR a \code{list} of \code{GRanges} as created by the
#' \code{getGRangesFromMethylDiff} function. Each
#' entry of the \code{list} represents the differentially methylated results
#' for one generation (first entry = first genertation, second entry =
#' second generation, etc..). Each \code{GRanges} object holds statistics
#' for differentially methylated regions/bases.
#'
#' @return a \code{list} containing the following elements:
#' \itemize{
#' \item\code{i2} a \code{list} of \code{GRanges} Each
#' \code{GRanges} represents the intersection of analysis results between two
#' consecutive generations. The first element represents the intersection
#' of the
#' first and second generations; the second element, the intersection of
#' the second and third generations; etc.. The number of entries depends
#' of the number of generations.
#' \item\code{iAll} a \code{list} of \code{GRanges}. Each \code{GRanges}
#' represents the intersection fo the analysis results between three or more
#' consecutive generations. The first element represents the
#' intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number
#' of generations.
#' }
#'
#' @examples
#'
#' ## Load permutation results on sites
#' permutationResultsFile <- system.file("extdata",
#'     "permutationResultsForSites.RDS", package="methylInheritance")
#' permutationResults <- readRDS(permutationResultsFile)
#'
#' ## Transform result to GRanges
#' resultsGR <- methylInheritance:::getGRangesFromMethylDiff(methDiff =
#'     permutationResults, pDiff = 10, qvalue = 0.01, type = "hyper")
#'
#' ## Extract inter generational conserved sites
#' conservedSitesGR <- methylInheritance:::interGeneration(resultsGR)
#'
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom GenomicRanges intersect GRanges
#' @importFrom S4Vectors DataFrame values<- values
#' @keywords internal
interGeneration <- function(resultAllGenGR) {

    lInter <- list("i2" = list(), "iAll" = list())

    # Calculate intersection of two consecutive generations
    lInter$i2 <- lapply(2:length(resultAllGenGR), FUN = function(i,b){
        upM <- intersect(b[[i-1]][b[[i-1]]$meth.diff > 0],
                            b[[i]][b[[i]]$meth.diff > 0])
        downM <- intersect(b[[i-1]][b[[i-1]]$meth.diff < 0],
                            b[[i]][b[[i]]$meth.diff < 0])
        typeDiff <- DataFrame(typeDiff=rep(1,length(upM)))
        values(upM) <- cbind(values(upM), typeDiff)
        typeDiff <- DataFrame(typeDiff=rep(-1,length(downM)))
        values(downM) <- cbind(values(downM), typeDiff)
        c(upM, downM)
    }, b = resultAllGenGR)

    # Calculate intersection of three or more consercutive generations
    cur <- lInter$i2[[1]]
    for(i in 3:length(resultAllGenGR)){
        upM <- intersect(cur[cur$typeDiff > 0],
                        resultAllGenGR[[i]][resultAllGenGR[[i]]$meth.diff > 0])
        downM <- intersect(cur[cur$typeDiff < 0],
                            resultAllGenGR[[i]][
                            resultAllGenGR[[i]]$meth.diff < 0])
        typeDiff <- DataFrame(typeDiff=rep(1,length(upM)))
        values(upM) <- cbind(values(upM), typeDiff)
        typeDiff <- DataFrame(typeDiff=rep(-1,length(downM)))
        values(downM) <- cbind(values(downM), typeDiff)

        lInter$iAll[[i-2]] <- c(upM,downM)
        cur <- lInter$iAll[[i-2]]
    }

    return(lInter)
}


#' @title Create directories that will contained the results of the
#' permutations in RDS format
#'
#' @description Create directories that will contained the results of the
#' permutations in RDS format.
#'
#' @param outputDir a string of \code{character}, the name of the main
#' directory to be created.
#'
#' @param doingSites a \code{logical}, a directory consecrated to contain the
#' results of the permutation analysis for sites is created when
#' \code{doingSites} = \code{TRUE}. Default: \code{TRUE}.
#'
#' @param doingTiles a \code{logical}, a directory consecrated to contain the
#' results of the permutation analysis for tiles is created when
#' \code{doingTiles} = \code{TRUE}. Default: \code{FALSE}.
#'
#' @param saveInfoByGeneration a \code{logical}, when \code{TRUE}, the
#' information about differentially methylated sites and tiles for each
#' generation is saved in a RDS file. The information is saved in a different
#' file for each permutation.
#'
#' @return \code{0} when all directories are created without problem.
#'
#' @examples
#'
#' ## Create an output directory for SITES only
#' methylInheritance:::createOutputDir(outputDir = "testSites",
#'     doingSites = TRUE, doingTiles = FALSE, saveInfoByGeneration = TRUE)
#'
#' @author Astrid Deschenes
#' @keywords internal
createOutputDir <- function(outputDir, doingSites = TRUE,
                                doingTiles = FALSE,
                            saveInfoByGeneration) {

    # Create directories for output files
    if (!dir.exists(outputDir)) {
        dir.create(outputDir, showWarnings = TRUE)
    }

    if (doingSites) {
        type <-  "SITES"
        dirName <- paste0(outputDir, type)
        if (!dir.exists(dirName)) {
            dir.create(dirName, showWarnings = TRUE)
        }
    }

    if (doingTiles) {
        type <-  "TILES"
        dirName <- paste0(outputDir, type)
        if (!dir.exists(dirName)) {
            dir.create(dirName, showWarnings = TRUE)
        }
    }

    # Create directory for the information for each generation
    if (saveInfoByGeneration) {
        dirName <- paste0(outputDir, "InfoByGeneration")
        if (!dir.exists(dirName)) {
            dir.create(dirName, showWarnings = TRUE)
        }
    }
    return(0)
}


#' @title Run the analysis on one permutation dataset, including all
#' generations, using \code{methylKit} package
#'
#' @description Run CpG site or region analysis using the \code{methylKit}
#' package for each generation present in the dataset. The intersection of
#' conserved elements is obtained for each group of two consecutive
#' generations, as well as, for larger group subset. The output of the
#' analysis is saved in a RDS file when an directory is
#' specified.
#'
#' @param methylInfoForAllGenerations a \code{list} containing the
#' following elements:
#' \itemize{
#' \item \code{sample} a \code{list} of \code{methylRawList} entries, each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation. The number of generations must correspond to the number
#' of entries in the \code{sample}. At least 2 generations
#' must be present to do a permutation analysis.
#' \item \code{id} an integer, the permutation id.
#' }
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation or \code{NULL}. If the directory does not
#' exist, it will be created. When \code{NULL}, the results of the permutation
#' are not saved. Default: \code{NULL}.
#'
#' @param nbrCoresDiffMeth a positive integer, the number of cores to use for
#' parallel differential methylation calculations.Parameter used for both
#' sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the
#' package \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the  \code{methylKit} package.
#'
#' @param minMethDiff a positive integer betwwen [0,100], the absolute value
#' of methylation percentage change between cases and controls. The parameter
#' correspond to the \code{difference} parameter in the
#' package \code{methylKit}. Default: \code{10}.
#'
#' @param qvalue a positive \code{double} inferior to \code{1}, the cutoff
#' for qvalue of differential methylation statistic. Default: \code{0.01}.
#'
#' @param maxPercReads a double between [0-100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#' Default: \code{99.9}.
#'
#' @param destrand a logical, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both
#' sites and tiles analysis. Default: \code{FALSE}.
#'
#' @param minCovBasesForTiles a non-negative integer, the minimum number of
#' bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the
#' package \code{methylKit}.
#' Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive integer, the size of the tiling window. The
#' parameter corresponds to the \code{win.size} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param stepSize a positive integer, the step size of tiling windows. The
#' parameter corresponds to the \code{stepSize} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param restartCalculation a \code{logical}, when \code{TRUE}, only
#' permutations that don't have a RDS result final are run.
#'
#' @param saveInfoByGeneration a \code{logical}, when \code{TRUE}, the
#' information about differentially methylated sites and tiles for each
#' generation is saved in a RDS file. The information is saved in a different
#' file for each permutation. The files are only saved when the
#' \code{outputDir} is not \code{NULL}.
#'
#' @return a \code{list} containing the following elements:
#' \itemize{
#' \item \code{SITES} Only present when \code{type} = \code{"sites"} or
#' \code{"both"}, a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc..
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between two consecutive generations.The
#' first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc..
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number
#' of generations.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number of
#' generations.
#' }
#' }
#' \item \code{TILES} Only present when \code{type} = \code{"tiles"} or
#' \code{"both"}, a \code{list} containing:
#' itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated positions between two consecutive
#' generations. The first element represents the intersection of the
#' first and second generations; the second element, the intersection of
#' the second and third generations; etc..
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated positions between two consecutive
#' generations.The first element represents the intersection of the first and
#' second generations; the second element, the intersection of the second
#' and third generations; etc..
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated positions between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number
#' of generations.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated positions between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number of
#' generations.
#' }
#' }
#' }
#'
#' @examples
#'
#' ## Load methyl information
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## Run a permutation analysis
#' methylInheritance:::runOnePermutationOnAllGenerations(id = 2,
#'     methylInfoForAllGenerations = samplesForTransgenerationalAnalysis,
#'     type = "tiles", outputDir = NULL,
#'     nbrCoresDiffMeth = 1, minReads = 10, minMethDiff = 10, qvalue = 0.01,
#'     maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 0,
#'     tileSize = 1000, stepSize = 1000, restartCalculation = FALSE)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom methylKit filterByCoverage normalizeCoverage unite
#' calculateDiffMeth getMethylDiff getData tileMethylCounts methRead
#' @importFrom GenomicRanges width
#' @keywords internal
runOnePermutationOnAllGenerations <- function(id,
                        methylInfoForAllGenerations,
                        type = c("both", "sites", "tiles"),
                        outputDir = NULL,
                        nbrCoresDiffMeth = 1,
                        minReads = 10, minMethDiff = 10,
                        qvalue = 0.01, maxPercReads = 99.9,
                        destrand = FALSE, minCovBasesForTiles = 0,
                        tileSize = 1000, stepSize = 1000,
                        restartCalculation, saveInfoByGeneration) {

    # Validate type value
    type <- match.arg(type)

    doTiles <- any(type %in% c("tiles", "both"))
    doSites <- any(type %in% c("sites", "both"))

    ## Extract info from input list
    if (id > 0) {
        methylRawForAllGenerations <- formatInputMethylData(methylKitData =
                                            methylInfoForAllGenerations)
    } else {
        methylRawForAllGenerations <- methylInfoForAllGenerations
    }
    nbrGenerations <- length(methylRawForAllGenerations)

    ## Preparing list that will receive final results
    permutationList <- list()
    if (doTiles) {
        permutationList[["TILES"]] <- list()
    }
    if (doSites) {
        permutationList[["SITES"]] <- list()
    }

    readyTiles <- FALSE
    if (doTiles && restartCalculation) {
        readyTiles <- isInterGenerationResults(outputDir, id, "tiles")
    }

    readySites <- FALSE
    if (doSites && restartCalculation) {
        readySites <- isInterGenerationResults(outputDir, id, "sites")
    }

    for (i in 1:nbrGenerations) {

        allSamplesForOneGeneration <- methylRawForAllGenerations[[i]]

        ## SITES
        if (doSites && !readySites) {

            ## Filter sites by coverage
            filtered.sites <- filterByCoverage(allSamplesForOneGeneration,
                                                lo.count = minReads,
                                                lo.perc = NULL,
                                                hi.count = NULL,
                                                hi.perc = maxPercReads)

            ## Normalize coverage
            filtered.sites <- normalizeCoverage(filtered.sites, "median")

            ## Merge all samples to one table
            meth.sites <- unite(filtered.sites, destrand = destrand)

            if (length(meth.sites@.Data[[1]]) == 0) {
                stop("meth.sites IS EMPTY")
            }

            ## Get differentially methylated sites
            allSites <- suppressWarnings(
                calculateDiffMeth(meth.sites, mc.cores = nbrCoresDiffMeth))

            permutationList[["SITES"]][[i]] <- suppressWarnings(
                getMethylDiff(allSites, difference = minMethDiff,
                                qvalue = qvalue))
        }

        ## TILES
        if (doTiles && !readyTiles) {

            ## Summarize methylated base counts over tilling windows
            tiles <- tileMethylCounts(allSamplesForOneGeneration,
                                        win.size = tileSize,
                                        step.size = stepSize,
                                        cov.bases = minCovBasesForTiles)

            ## Filter tiles by coverage
            filtered.tiles <- filterByCoverage(tiles,
                                                lo.count = minReads,
                                                lo.perc = NULL,
                                                hi.count = NULL,
                                                hi.perc = maxPercReads)

            ## Normalize coverage
            filtered.tiles <- normalizeCoverage(filtered.tiles, "median")

            ## Merge all samples to one table
            meth.tiles <- unite(filtered.tiles, destrand = destrand)

            ## Get diff methylated tiles
            allTiles <- suppressWarnings(
                calculateDiffMeth(meth.tiles, mc.cores = nbrCoresDiffMeth))

            permutationList[["TILES"]][[i]] <- suppressWarnings(
                getMethylDiff(allTiles, difference = minMethDiff,
                              qvalue = qvalue))
        }
    }

    ## Save all results per generation in RDS file when specified
    if (!is.null(outputDir) && saveInfoByGeneration) {
        saveRDS(object = permutationList, file = paste0(outputDir,
                        "InfoByGeneration/DMEByGeneration_", id, ".RDS"))
    }

    permutationFinal <- list()

    ## Calculate the number of SITES in the intersection
    if (doSites) {
        if (!readySites) {
            ## Transform initial results to GRanges
            resultGR <- getGRangesFromMethylDiff(permutationList[["SITES"]],
                                        minMethDiff, qvalue, type = "all")

            ## Extract inter generational conserved sites
            result <- interGeneration(resultGR)

            ## Save results in RDS file when specified
            if (!is.null(outputDir)) {
                saveInterGenerationResults(outputDir, id, type = "sites",
                                            result)
            }
        } else {
            result<- readInterGenerationResults(outputDir, id, type = "sites")
        }

        ## Create list that will contain final results
        # permutationFinal[["SITES"]] <- list()
        # permutationFinal[["SITES"]][["i2"]] <- list()
        # permutationFinal[["SITES"]][["i2"]][["HYPER"]] <- list()
        # permutationFinal[["SITES"]][["i2"]][["HYPO"]]  <- list()
        # permutationFinal[["SITES"]][["iAll"]][["HYPER"]]  <- list()
        # permutationFinal[["SITES"]][["iAll"]][["HYPO"]]   <- list()
        #
        # permutationFinal[["SITES"]][["i2"]][["HYPER"]] <- lapply(result$i2,
        #                 FUN = function(x) {sum(width(x[x$typeDiff > 0]))})
        #
        # permutationFinal[["SITES"]][["i2"]][["HYPO"]]  <- lapply(result$i2,
        #                 FUN = function(x) {sum(width(x[x$typeDiff < 0]))})
        #
        # permutationFinal[["SITES"]][["iAll"]][["HYPER"]] <- lapply(
        #                 result$iAll,
        #                 FUN = function(x) {sum(width(x[x$typeDiff > 0]))})
        #
        # permutationFinal[["SITES"]][["iAll"]][["HYPO"]]  <- lapply(
        #                 result$iAll,
        #                 FUN = function(x) {sum(width(x[x$typeDiff < 0]))})
    }

    ## Calculate the number of TILES in the intersection
    if (doTiles) {
        if (!readyTiles) {
            ## Transform initial results to GRanges
            resultGR <- getGRangesFromMethylDiff(permutationList[["TILES"]],
                                minMethDiff, qvalue, type = "all")

            ## Extract inter generational conserved tiles
            result <- interGeneration(resultGR)

            ## Save results in RDS file when specified
            if (!is.null(outputDir)) {
                saveInterGenerationResults(outputDir, id, type = "tiles",
                                                result)
            }
        } else {
            result <- readInterGenerationResults(outputDir, id, type = "tiles")
        }

        ## Create list that will contain final results
        # permutationFinal[["TILES"]] <- list()
        # permutationFinal[["TILES"]][["i2"]] <- list()
        # permutationFinal[["TILES"]][["i2"]][["HYPER"]] <- list()
        # permutationFinal[["TILES"]][["i2"]][["HYPO"]]  <- list()
        # permutationFinal[["TILES"]][["iAll"]][["HYPER"]]  <- list()
        # permutationFinal[["TILES"]][["iAll"]][["HYPO"]]   <- list()
        #
        # permutationFinal[["TILES"]][["i2"]][["HYPER"]] <- lapply(result$i2,
        #                 FUN = function(x) {sum(width(x[x$typeDiff > 0]))})
        #
        # permutationFinal[["TILES"]][["i2"]][["HYPO"]]  <- lapply(result$i2,
        #                 FUN = function(x) {sum(width(x[x$typeDiff < 0]))})
        #
        # permutationFinal[["TILES"]][["iAll"]][["HYPER"]] <- lapply(
        #                 result$iAll,
        #                 FUN = function(x) {sum(width(x[x$typeDiff > 0]))})
        #
        # permutationFinal[["TILES"]][["iAll"]][["HYPO"]]  <- lapply(
        #                 result$iAll,
        #                 FUN = function(x) {sum(width(x[x$typeDiff < 0]))})
    }

    return(permutationFinal)
}


#' @title Save the result of on CpG site or tile analysis on all generations.
#' The anaysis can come from observed or shuffled dataset. Each case is
#' saved with a different extension.
#'
#' @description Save the result of on CpG site or tile analysis on all
#' generations. The results are saved in a RDS file. The anaysis can have been
#' done on the observed or shuffled dataset.
#' Each permutation is saved using its identifiant in the file name.
#'
#' @param outputDir a string of \code{character}, the name of the directory
#' that will contain
#' the results of the permutation. The name should end with a slash. The
#' directory should already exists.
#'
#' @param permutationID an \code{integer}, the identifiant of the permutation.
#' When the \code{permutationID} = \code{0}, the results are considered as the
#' observed results and are saved in a file with the "_observed_results.RDS"
#' extension. When the \code{permutationID} != \code{0}, the results are
#' considered as permutation results and are saved in a file with the
#' "_permutation_{permutationID}.RDS" extension.
#'
#' @param type One of the \code{"sites"} or \code{"tiles"} strings. Specifies
#' the type of differentially methylated elements should be saved.
#' Default: \code{"sites"}.
#'
#' @param interGenerationResult a \code{list} that corresponds to the output
#' of the \code{interGeneration} function, the result of on CpG site or tile
#' analysis on all generations.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load permutation results on sites
#'
#' permutationResultsFile <- system.file("extdata",
#'     "permutationResultsForSites.RDS", package="methylInheritance")
#' permutationResults <- readRDS(permutationResultsFile)
#'
#' ## Transform result to GRanges
#' resultsGR <- methylInheritance:::getGRangesFromMethylDiff(methDiff =
#'     permutationResults, pDiff = 10, qvalue = 0.01, type = "hyper")
#'
#' ## Extract inter-generationally conserved sites
#' interGenerationResult <- methylInheritance:::interGeneration(resultsGR)
#'
#' ## Create directories
#' dir.create("TEST", showWarnings = TRUE)
#' dir.create("TEST/SITES", showWarnings = TRUE)
#'
#' ## Save results
#' methylInheritance:::saveInterGenerationResults(
#'     outputDir = "TEST/", permutationID=100, type = "sites",
#'     interGenerationResult = interGenerationResult)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @keywords internal
saveInterGenerationResults <- function(outputDir, permutationID,
                                        type = c("sites", "tiles"),
                                        interGenerationResult) {

    if (permutationID != 0) {
        ## Save the permutation results
        saveRDS(object = interGenerationResult,
                file = paste0(outputDir,  toupper(type), "/",
                        toupper(type), "_permutation_", permutationID, ".RDS"))
    } else {
        ## Save the observed results
        saveRDS(object = interGenerationResult,
                    file = paste0(outputDir, toupper(type), "/",
                        toupper(type), "_observed_results.RDS"))
    }

    return(0)
}


#' @title Verify if a specific file containing intergenerational results
#' exists or not.
#'
#' @description Verify if a specific file containing intergenerational results
#' exists or not.
#'
#' @param outputDir a string of \code{character}, the name of the directory
#' that will contain
#' the results of the permutation. The name should end with a slash. The
#' directory should already exists.
#'
#' @param permutationID an \code{integer}, the identifiant of the permutation.
#' When the \code{permutationID} = \code{0}, the results are considered as the
#' observed results and are saved in a file with the "_observed_results.RDS"
#' extension. When the \code{permutationID} != \code{0}, the results are
#' considered as permutation results and are saved in a file with the
#' "_permutation_{permutationID}.RDS" extension.
#'
#' @param type One of the \code{"sites"} or \code{"tiles"} strings. Specifies
#' the type of differentially methylated elements should be saved.
#' Default: \code{"sites"}.
#'
#' @return \code{TRUE} when file present; otherwise \code{FALSE}.
#'
#' @examples
#'
#' ## Get the name of the directory where the file is stored
#' filesDir <- system.file("extdata", "TEST", package="methylInheritance")
#'
#' ## Verify that DMS intergenerational results for the observed data exists
#' methylInheritance:::isInterGenerationResults(outputDir =
#'     paste0(filesDir, "/"), 0, "sites")
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @keywords internal
isInterGenerationResults <- function(outputDir, permutationID,
                                     type = c("sites", "tiles")) {

    if (permutationID != 0) {
        result <- file.exists(paste0(outputDir,  toupper(type), "/",
                        toupper(type), "_permutation_", permutationID, ".RDS"))
    } else {
        result <- file.exists(paste0(outputDir, toupper(type), "/",
                                     toupper(type), "_observed_results.RDS"))
    }

    return(result)
}

#' @title Read and return intergenerational results contained in a
#' RDS file
#'
#' @description Read and return intergenerational results contained in a
#' RDS file
#'
#' @param outputDir a string of \code{character}, the name of the directory
#' that will contain
#' the results of the permutation. The name should end with a slash. The
#' directory should already exists.
#'
#' @param permutationID an \code{integer}, the identifiant of the permutation.
#' When the \code{permutationID} = \code{0}, the results are considered as the
#' observed results and are saved in a file with the "_observed_results.RDS"
#' extension. When the \code{permutationID} != \code{0}, the results are
#' considered as permutation results and are saved in a file with the
#' "_permutation_{permutationID}.RDS" extension.
#'
#' @param type One of the \code{"sites"} or \code{"tiles"} strings. Specifies
#' the type of differentially methylated elements should be saved.
#' Default: \code{"sites"}.
#'
#' @return a \code{list} containing the intergenerational results for the
#' specified permutation.
#'
#' @examples
#'
#' ## Get the name of the directory where the file is stored
#' filesDir <- system.file("extdata", "TEST", package="methylInheritance")
#'
#' ## Read DMS intergenerational results for the observed data
#' methylInheritance:::readInterGenerationResults(outputDir =
#'     paste0(filesDir, "/"), 0, "sites")
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @keywords internal
readInterGenerationResults <- function(outputDir, permutationID,
                                     type = c("sites", "tiles")) {

    if (permutationID != 0) {
        result <- readRDS(file = paste0(outputDir,  toupper(type), "/",
                    toupper(type), "_permutation_", permutationID, ".RDS"))
    } else {
        result <- readRDS(file = paste0(outputDir, toupper(type), "/",
                                     toupper(type), "_observed_results.RDS"))
    }

    return(result)
}


#' @title Extract the number of conserved differentially methylated
#' elements in \code{GRanges}.
#'
#' @description Extract the number of conserved differentially methylated
#' elements in \code{GRanges}. Each \code{GRanges}
#' is the result of one intersection between two or more consecutive
#' generations for one analysis done on all generations.
#' The hypo and hyper differentially methylated elements are counted
#' separatly.
#'
#' @param interGenerationGR a \code{list} that contains the information for
#' all differentially methylated analysis done on each generation present in
#' the initial dataset. The \code{list} must contain the following elements:
#' \itemize{
#' \item\code{i2} a \code{list} of \code{GRanges} Each
#' \code{GRanges} represents the intersection of analysis results between two
#' consecutive generations. The first element represents the intersection
#' of the
#' first and second generations; the second element, the intersection of
#' the second and third generations; etc.. The number of entries depends
#' of the number of generations.
#' \item\code{iAll} a \code{list} of \code{GRanges}. Each \code{GRanges}
#' represents the intersection fo the analysis results between three or more
#' consecutive generations. The first element represents the
#' intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number
#' of generations.
#' }
#'
#' @return a \code{list} containing the following elements:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc..
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between two consecutive generations.The
#' first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc..
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number
#' of generations.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number of
#' generations.
#' }
#' }
#'
#' @examples
#'
#' ## Get the name of the directory where the file is stored
#' filesDir <- system.file("extdata", "TEST", package="methylInheritance")
#'
#' ## Load file containing results from a observation analysis
#' obsResults <- readRDS(file = paste0(filesDir,
#'     "/SITES/SITES_observed_results.RDS"))
#'
#' ## Create data structure using information form the observation analysis
#' formatedResults <- methylInheritance:::createDataStructure(obsResults)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom GenomicRanges width
#' @keywords internal
createDataStructure <- function(interGenerationGR) {

    result <- list()
    result[["i2"]] <- list()
    result[["i2"]][["HYPER"]] <- list()
    result[["i2"]][["HYPO"]]  <- list()
    result[["iAll"]][["HYPER"]]  <- list()
    result[["iAll"]][["HYPO"]]   <- list()
    result[["i2"]][["HYPER"]] <- lapply(interGenerationGR$i2,
                        FUN = function(x) {sum(width(x[x$typeDiff > 0]))})
    result[["i2"]][["HYPO"]]  <- lapply(interGenerationGR$i2,
                        FUN = function(x) {sum(width(x[x$typeDiff < 0]))})
    result[["iAll"]][["HYPER"]] <- lapply(interGenerationGR$iAll,
                        FUN = function(x) {sum(width(x[x$typeDiff > 0]))})
    result[["iAll"]][["HYPO"]]  <- lapply(interGenerationGR$iAll,
                        FUN = function(x) {sum(width(x[x$typeDiff < 0]))})

    return(result)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param methylKitData TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom methods new
#' @keywords internal
formatInputMethylData <- function(methylKitData) {

    ## Extract information
    nbGenerations <- length(methylKitData)
    nbSamplesByGeneration <- sapply(methylKitData, length)
    nbSamples  <- sum(nbSamplesByGeneration)
    allSamples <- unlist(methylKitData, recursive = FALSE)

    ## Random sample
    permutationSample <- sample(seq_len(nbSamples))

    ## Create list that will contain information for all generations
    ## related to the same permutation analysis
    permutationList <- list()
    start <- 1
    for (j in 1:nbGenerations) {
        end <- start + nbSamplesByGeneration[j] - 1
        samplePos <- permutationSample[start:end]
        treatment <- methylKitData[[j]]@treatment
        newSampleList <- new("methylRawList", allSamples[samplePos],
                                treatment = treatment)
        permutationList[[j]] <- newSampleList
        start <- end + 1
    }

    return(permutationList)
}
