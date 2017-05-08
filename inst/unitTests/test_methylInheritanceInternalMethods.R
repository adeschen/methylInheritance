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
