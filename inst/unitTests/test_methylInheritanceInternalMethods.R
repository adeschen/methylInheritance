###################################################
# Created by Astrid Deschenes
# 2017-01-09
###################################################

###################################################
## Test the methylInheritanceInternalMethods functions
###################################################

METHYL_OBJ_FILE <- system.file("extdata", "methylObj_001.RDS",
                                package = "methylInheritance")

METHYL_OBJ <- readRDS(METHYL_OBJ_FILE)

data(methylInheritanceResults)

###########################################################
## runOnePermutationOnAllGenerations() function
###########################################################

## Test sites when all parameters are valid
# test.validateRunPermutationUsingMethylKitInfo_sites_good_01 <- function() {
#     ## Extract information
#     set.seed(111)
#     allSamples <- sample(unlist(METHYL_OBJ, recursive = FALSE), 36, replace = F)
#     treatment <- c(0,0,0,0,0,0,1,1,1,1,1,1)
#     sampleList01 <- new("methylRawList", allSamples[1:12],
#                         treatment = treatment)
#     sampleList02 <- new("methylRawList", allSamples[13:24],
#                         treatment = treatment)
#     sampleList03 <- new("methylRawList", allSamples[25:36],
#                         treatment = treatment)
#     input <- list(sampleList01, sampleList02, sampleList03)
#
#     obs <- tryCatch(methylInheritance:::runOnePermutationOnAllGenerations(
#         id = 1, methylInfoForAllGenerations = input, outputDir = NULL, type = "sites",
#         nbrCoresDiffMeth = 1,
#         minReads = 10, minMethDiff = 10, qvalue = 0.05,
#         maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
#         tileSize = 1000, stepSize = 100, restartCalculation = FALSE),
#         error=conditionMessage)
#
#     exp <- list()
#     exp[["SITES"]] <- list()
#     exp[["SITES"]][["i2"]] <- list()
#     exp[["SITES"]][["i2"]][["HYPER"]] <- list(1,0)
#     exp[["SITES"]][["i2"]][["HYPO"]]  <- list(0,0)
#     exp[["SITES"]][["iAll"]][["HYPER"]]  <- list(0)
#     exp[["SITES"]][["iAll"]][["HYPO"]]   <- list(0)
#
#     message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_sites_good_01() ",
#                     "- Valid parameters did not generated expected results.")
#     checkEquals(obs, exp, msg = message)
# }

## Test tiles when all parameters are valid
# test.validateRunPermutationUsingMethylKitInfo_tiles_good_01 <- function() {
#     ## Extract information
#     set.seed(11)
#     allSamples <- sample(unlist(METHYL_OBJ, recursive = FALSE), 36,
#                             replace = F)
#     treatment <- c(0,0,0,0,0,0,1,1,1,1,1,1)
#     sampleList01 <- new("methylRawList", allSamples[1:12],
#                             treatment = treatment)
#     sampleList02 <- new("methylRawList", allSamples[13:24],
#                             treatment = treatment)
#     sampleList03 <- new("methylRawList", allSamples[25:36],
#                             treatment = treatment)
#     input <- list(sampleList01, sampleList02, sampleList03)
#
#     obs <- tryCatch(methylInheritance:::runOnePermutationOnAllGenerations(
#         id = 2,
#         methylInfoForAllGenerations = input, outputDir = NULL, type = "tiles",
#         nbrCoresDiffMeth = 1,
#         minReads = 5, minMethDiff = 5, qvalue = 0.05,
#         maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
#         tileSize = 1000, stepSize = 100, restartCalculation = FALSE),
#         error=conditionMessage)
#
#     exp <- list()
#     exp[["TILES"]] <- list()
#     exp[["TILES"]][["i2"]] <- list()
#     exp[["TILES"]][["i2"]][["HYPER"]] <- list(0, 0)
#     exp[["TILES"]][["i2"]][["HYPO"]]  <- list(1900, 0)
#     exp[["TILES"]][["iAll"]][["HYPER"]]  <- list(0)
#     exp[["TILES"]][["iAll"]][["HYPO"]]   <- list(0)
#
#     message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_tiles_good_01() ",
#                 "- Valid parameters did not generated expected results.")
#
#     checkEquals(obs, exp, msg = message)
# }

###########################################################
## isInterGenerationResults() function
###########################################################

test.isInterGenerationResults_true <- function() {

    filesDir <- system.file("extdata", "TEST", package="methylInheritance")

    obs <- methylInheritance:::isInterGenerationResults(outputDir =
                paste0(filesDir, "/"), 0, "sites")

    message <- paste0("test.isInterGenerationResults_true() ",
                      "- Function should return TRUE")

    checkTrue(obs, msg = message)
}


###########################################################
## validateExtractInfo() function
###########################################################

test.validateExtractInfo_position_zero <- function() {
    obs <- tryCatch(methylInheritance:::validateExtractInfo(
        allResults = methylInheritanceResults, type = "sites",
        inter = "i2", position = 0),
        error=conditionMessage)

    exp <- "position must be a positive integer"

    message <- paste0("test.validateExtractInfo_position_zero() - ",
                      "Zero position value did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateExtractInfo_position_string <- function() {
    obs <- tryCatch(methylInheritance:::validateExtractInfo(
        allResults = methylInheritanceResults, type = "sites",
        inter = "i2", position = "hi"),
        error=conditionMessage)

    exp <- "position must be a positive integer"

    message <- paste0("test.validateExtractInfo_position_string() - ",
                      "Zero position value did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateExtractInfo_allResults_vector <- function() {
    obs <- tryCatch(methylInheritance:::validateExtractInfo(
        allResults = c(1,2,3), type = "sites",
        inter = "i2", position = 1),
        error=conditionMessage)

    exp <- "allResults must be of class \"methylInheritanceAllResults\""

    message <- paste0("test.validateExtractInfo_allResults_vector() - ",
                      "allResults vector did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateExtractInfo_type_wrong <- function() {
    obs <- tryCatch(methylInheritance:::validateExtractInfo(
        allResults = methylInheritanceResults, type = "toto",
        inter = "i2", position = 1),
        error=conditionMessage)

    exp <- "allResults must have an element called \"TOTO\" in its \"OBSERVATION\" list"

    message <- paste0("test.validateExtractInfo_type_wrong() - ",
                      "Wrong type did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateExtractInfo_type_wrong <- function() {
    obs <- tryCatch(methylInheritance:::validateExtractInfo(
        allResults = methylInheritanceResults, type = "sites",
        inter = "hi", position = 1),
        error=conditionMessage)

    exp <- "allResults must have an element called \"hi\" in the \"SITES\" list present in its \"OBSERVATION\" list"

    message <- paste0("test.validateExtractInfo_type_wrong() - ",
                      "Wrong inter did not generated expected message.")

    checkEquals(obs, exp, message)
}

###########################################################
## calculateSignificantLevel() function
###########################################################

test.calculateSignificantLevel_true <- function() {

    filesDir <- system.file("extdata", "TEST", package="methylInheritance")

    results <- loadAllRDSResults(analysisResultsDir = filesDir,
                      permutationResultsDir = filesDir, doingSites = TRUE,
                      doingTiles = FALSE)

    iAll <- extractInfo(results, type = "sites", inter = "iAll", 1)

    obs <- methylInheritance:::calculateSignificantLevel(iAll)

    message <- paste0("test.calculateSignificantLevel_true() ",
                      "- Function did not return expected values")

    exp <- list(HYPER=1.0, HYPO=(2.0/4.0))

    checkEquals(obs, exp, msg = message)
}
