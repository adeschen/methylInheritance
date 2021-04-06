###################################################
# Created by Astrid Deschenes
# 2017-01-05
###################################################

###################################################
## Test the methylInheritanceMethods functions
###################################################

METHYL_OBJ_FILE_01 <- system.file("extdata", "methylObj_001.RDS",
                                    package = "methylInheritance")

METHYL_OBJ_01 <- readRDS(METHYL_OBJ_FILE_01)

TEST_DIR <- system.file("extdata", "TEST", package = "methylInheritance")

data("methylInheritanceResults")

.tearDown <- function() {
    if (dir.exists("test_002")) {
        unlink("test_002", recursive = TRUE, force = TRUE)
    }
}

###########################################################
## runPermutation() function
###########################################################

## Test when methylKitData is not a valid RDS file name
test.runPermutation_methylKitData_not_valid_RDS <- function() {
    obs <- tryCatch(runPermutation(
        methylKitData = "HI",  outputDir = "test_002",
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "The file \"HI\" does not exist."

    message <- paste0(" test.runPermutation_methylKitData_not_valid_RDS() ",
                      "- Not valid file for methylKitData did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when all parameters valid
test.runPermutation_good_001 <- function() {

    if (!dir.exists("test_002")) {
        dir.create("test_002/SITES/", recursive = TRUE)
        dir.create("test_002/TILES/", recursive = TRUE)
    }

    obs <- runPermutation(methylKitData = METHYL_OBJ_FILE_01, runObservationAnalysis = FALSE,
                    type = "both", nbrPermutations = 2, minReads = 5, minMethDiff = 5,
                    outputDir = "test_002", vSeed = 2000)

    obsV <- methylInheritance::loadAllRDSResults(permutationResultsDir = "test_002",
                                                analysisResultsDir = NULL, doingSites = TRUE,
                                                doingTiles = TRUE)

    exp <- list()
    exp[["PERMUTATION"]] <- list()
    exp[["PERMUTATION"]][[1]] <- list()
    exp[["PERMUTATION"]][[1]][["SITES"]] <- list()
    exp[["PERMUTATION"]][[1]][["SITES"]][["i2"]] <- list()
    exp[["PERMUTATION"]][[1]][["SITES"]][["i2"]][["HYPER"]] <- list(0,0)
    exp[["PERMUTATION"]][[1]][["SITES"]][["i2"]][["HYPO"]] <- list(1,3)
    exp[["PERMUTATION"]][[1]][["SITES"]][["iAll"]] <- list()
    exp[["PERMUTATION"]][[1]][["SITES"]][["iAll"]][["HYPER"]] <- list(0)
    exp[["PERMUTATION"]][[1]][["SITES"]][["iAll"]][["HYPO"]] <- list(1)
    exp[["PERMUTATION"]][[1]][["TILES"]] <- list()
    exp[["PERMUTATION"]][[1]][["TILES"]][["i2"]] <- list()
    exp[["PERMUTATION"]][[1]][["TILES"]][["i2"]][["HYPER"]] <- list(0,0)
    exp[["PERMUTATION"]][[1]][["TILES"]][["i2"]][["HYPO"]] <- list(1000,0)
    exp[["PERMUTATION"]][[1]][["TILES"]][["iAll"]] <- list()
    exp[["PERMUTATION"]][[1]][["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    exp[["PERMUTATION"]][[1]][["TILES"]][["iAll"]][["HYPO"]] <- list(0)
    exp[["PERMUTATION"]][[2]] <- list()
    exp[["PERMUTATION"]][[2]][["SITES"]] <- list()
    exp[["PERMUTATION"]][[2]][["SITES"]][["i2"]] <- list()
    exp[["PERMUTATION"]][[2]][["SITES"]][["i2"]][["HYPER"]] <- list(0,0)
    exp[["PERMUTATION"]][[2]][["SITES"]][["i2"]][["HYPO"]] <- list(0,1)
    exp[["PERMUTATION"]][[2]][["SITES"]][["iAll"]] <- list()
    exp[["PERMUTATION"]][[2]][["SITES"]][["iAll"]][["HYPER"]] <- list(0)
    exp[["PERMUTATION"]][[2]][["SITES"]][["iAll"]][["HYPO"]] <- list(0)
    exp[["PERMUTATION"]][[2]][["TILES"]] <- list()
    exp[["PERMUTATION"]][[2]][["TILES"]][["i2"]] <- list()
    exp[["PERMUTATION"]][[2]][["TILES"]][["i2"]][["HYPER"]] <- list(0,0)
    exp[["PERMUTATION"]][[2]][["TILES"]][["i2"]][["HYPO"]] <- list(0,0)
    exp[["PERMUTATION"]][[2]][["TILES"]][["iAll"]] <- list()
    exp[["PERMUTATION"]][[2]][["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    exp[["PERMUTATION"]][[2]][["TILES"]][["iAll"]][["HYPO"]] <- list(0)

    message <- paste0(" test.runPermutation_good_001() ",
                      "- Valid parameters did not generated expected message.")

    checkEquals(obs, 0, msg = message)
    checkEquals(obsV$PERMUTATION, exp$PERMUTATION, msg = message)
}


###########################################################
# runObservation() function
###########################################################

## Test when methylKitData is not a valid RDS file name
test.runObservation_methylKitData_not_valid <- function() {
    obs <- tryCatch(runObservation(
        methylKitData = "ALLO",  outputDir = NULL, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "The file \"ALLO\" does not exist."

    message <- paste0(" test.runObservation_methylKitData_not_valid() ",
                      "- Not valid file for methylKitData did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when all parameters valid
test.runObservation_good_001 <- function() {

    if (!dir.exists("test_002")) {
        dir.create("test_002/SITES/", recursive = TRUE)
    }

    obs <- runObservation(
        methylKitData = METHYL_OBJ_FILE_01, type = "sites",
        outputDir = "test_002", nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 5, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 200, saveInfoByGeneration = FALSE)

    obsV <- methylInheritance::loadAllRDSResults(permutationResultsDir = NULL,
                                                 analysisResultsDir = "test_002")

    exp <- list()
    exp[["OBSERVATION"]] <- list()
    exp[["OBSERVATION"]][["SITES"]] <- list()
    exp[["OBSERVATION"]][["SITES"]][["i2"]] <- list()
    exp[["OBSERVATION"]][["SITES"]][["i2"]][["HYPER"]] <- list(0,0)
    exp[["OBSERVATION"]][["SITES"]][["i2"]][["HYPO"]] <- list(3,3)
    exp[["OBSERVATION"]][["SITES"]][["iAll"]] <- list()
    exp[["OBSERVATION"]][["SITES"]][["iAll"]][["HYPER"]] <- list(0)
    exp[["OBSERVATION"]][["SITES"]][["iAll"]][["HYPO"]] <- list(2)

    message <- paste0(" test.runObservation_good_001() ",
                      "- All valid parameters did not generated expected result.")

    checkEquals(obs, 0, msg = message)
    checkEquals(obsV$OBSERVATION, exp$OBSERVATION, msg = message)
}


###########################################################
# extractInfo() function
###########################################################

# Test result when all parameters are good
test.extractInfo_good_01 <- function() {
    obs <- tryCatch(extractInfo(allResults = methylInheritanceResults,
                    type = "sites", inter="i2", 1),
                    error=conditionMessage)

    exp <- data.frame(TYPE = rep(c("HYPO","HYPER"), 21),
                      RESULT = c(2,4,2,4,4,3,1,5,3,3,4,2,0,0,0,1,1,0,6,2,2,5,1,
                            3,2,4,222,67,6,4,183,53,1,6,34,102,2,2,4,3,2,2),
                      SOURCE = c("OBSERVATION", "OBSERVATION",
                                    rep("PERMUTATION", 40)))

    message <- paste0(" test.extractInfo_good_01() ",
                      "- Valid parameters for formatForGraph did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
## loadAllRDSResults() function
###########################################################

## Test result when all parameters are good
test.loadAllRDSResults_good_01 <- function() {
    obs <- tryCatch(loadAllRDSResults(analysisResultsDir = TEST_DIR,
                        permutationResultsDir = TEST_DIR, doingSites = TRUE,
                        doingTiles = TRUE),
                    error=conditionMessage)

    exp <- list()

    exp[["OBSERVATION"]] <- list("SITES" = list(), "TILES" = list())
    exp[["OBSERVATION"]][["SITES"]][["i2"]] <- list(HYPER=list(21, 10), HYPO=list(15, 12))
    exp[["OBSERVATION"]][["SITES"]][["iAll"]] <- list(HYPER=list(1), HYPO=list(3))
    i2 <- list(HYPER=list(2000, 3000), HYPO=list(2000, 3000))
    iAll <- list(HYPER=list(0), HYPO=list(0))
    exp[["OBSERVATION"]][["TILES"]][["i2"]] <- list(HYPER=list(2000, 3000), HYPO=list(2000, 3000))
    exp[["OBSERVATION"]][["TILES"]][["iAll"]] <- list(HYPER=list(0), HYPO=list(0))

    exp[["PERMUTATION"]] <- list()
    cas_01 <- list("SITES" = list(), "TILES" = list())
    cas_01[["SITES"]][["i2"]] <- list(HYPER=list(5, 7), HYPO=list(10, 11))
    cas_01[["SITES"]][["iAll"]] <- list(HYPER=list(2), HYPO=list(3))
    cas_01[["TILES"]][["i2"]] <- list(HYPER=list(0, 0), HYPO=list(1000, 3000))
    cas_01[["TILES"]][["iAll"]] <- list(HYPER=list(0), HYPO=list(1000))
    exp[["PERMUTATION"]][[1]] <- cas_01
    cas_02 <- list()
    cas_02[["SITES"]] <- list()
    cas_02[["SITES"]][["i2"]] <- list(HYPER=list(8, 9), HYPO=list(4, 7))
    cas_02[["SITES"]][["iAll"]] <- list(HYPER=list(3), HYPO=list(0))
    cas_02[["TILES"]] <- list()
    cas_02[["TILES"]][["i2"]] <- list(HYPER=list(1000, 1000), HYPO=list(1000, 1000))
    cas_02[["TILES"]][["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
    exp[["PERMUTATION"]][[2]] <- cas_02
    cas_03 <- list()
    cas_03[["SITES"]] <- list()
    cas_03[["SITES"]][["i2"]] <- list(HYPER=list(10, 7), HYPO=list(11, 7))
    cas_03[["SITES"]][["iAll"]] <- list(HYPER=list(1), HYPO=list(2))
    cas_03[["TILES"]] <- list()
    cas_03[["TILES"]][["i2"]] <- list(HYPER=list(0, 3000), HYPO=list(2000, 0))
    cas_03[["TILES"]][["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
    exp[["PERMUTATION"]][[3]] <- cas_03

    class(exp) <- "methylInheritanceAllResults"

    message <- paste0(" test.loadAllRDSResults_good_01() ",
                        "- Valid parameters for loadAllRDSResults() ",
                        "did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
## mergePermutationAndObservation() function
###########################################################

## Test when observationResults is not a list
test.mergePermutationAndObservation_observation_not_list <- function() {
    perm <- list()
    perm[["PERMUTATION"]] <- methylInheritanceResults$PERMUTATION

    obs <- tryCatch(mergePermutationAndObservation(permutationResults = perm,
                                observationResults = "33"),
                    error=conditionMessage)

    exp <- "observationResults must be a list"

    message <- paste0(" test.mergePermutationAndObservation_observation_not_list() ",
                      "- Not a list for observationResults did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}

## Test when permutationResults is not a list
test.mergePermutationAndObservation_permutation_not_list <- function() {
    res <- list()
    res[["OBSERVATION"]] <- methylInheritanceResults$OBSERVATION

    obs <- tryCatch(mergePermutationAndObservation(permutationResults = "allo",
                    observationResults = res),
                    error=conditionMessage)

    exp <- "permutationResults must be a list"

    message <- paste0(" test.mergePermutationAndObservation_permutation_not_list() ",
                      "- Not a list for permutationResults did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}

## Test result when all parameters are good
test.mergePermutationAndObservation_good_01 <- function() {

    observed <- list()
    observed[["OBSERVATION"]] <- list()
    observed[["OBSERVATION"]][["SITES"]] <- list()
    observed[["OBSERVATION"]][["SITES"]][["i2"]] <- list(HYPER=list(21, 10), HYPO=list(15, 12))
    observed[["OBSERVATION"]][["SITES"]][["iAll"]] <- list(HYPER=list(1), HYPO=list(3))
    observed[["OBSERVATION"]][["TILES"]] <- list()
    observed[["OBSERVATION"]][["TILES"]][["i2"]] <- list(HYPER=list(2000, 3000), HYPO=list(2000, 3000))
    observed[["OBSERVATION"]][["TILES"]][["iAll"]] <- list(HYPER=list(0), HYPO=list(0))

    permutated <- list()
    permutated[["PERMUTATION"]] <- list()
    cas_01 <- list()
    cas_01[["i2"]] <- list(HYPER=list(5, 7), HYPO=list(10, 11))
    cas_01[["iAll"]] <- list(HYPER=list(2), HYPO=list(3))
    cas_02 <- list()
    cas_02[["i2"]] <- list(HYPER=list(8, 9), HYPO=list(4, 7))
    cas_02[["iAll"]] <- list(HYPER=list(3), HYPO=list(0))
    cas_03 <- list()
    cas_03[["i2"]] <- list(HYPER=list(10, 7), HYPO=list(11, 7))
    cas_03[["iAll"]] <- list(HYPER=list(1), HYPO=list(2))
    permutated[["PERMUTATION"]][["SITES"]] <- list(cas_01, cas_02, cas_03)
    cas_01 <- list()
    cas_01[["i2"]] <- list(HYPER=list(0, 0), HYPO=list(1000, 3000))
    cas_01[["iAll"]] <- list(HYPER=list(0), HYPO=list(1000))
    cas_02 <- list()
    cas_02[["i2"]] <- list(HYPER=list(1000, 1000), HYPO=list(1000, 1000))
    cas_02[["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
    cas_03 <- list()
    cas_03[["i2"]] <- list(HYPER=list(0, 3000), HYPO=list(2000, 0))
    cas_03[["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
    permutated[["PERMUTATION"]][["TILES"]] <- list(cas_01, cas_02, cas_03)


    obs <- tryCatch(mergePermutationAndObservation(permutationResults = permutated,
                            observationResults = observed),
                    error=conditionMessage)

    exp <- list()
    exp[["PERMUTATION"]] <- permutated[["PERMUTATION"]]
    exp[["OBSERVATION"]] <- observed[["OBSERVATION"]]
    class(exp) <- "methylInheritanceAllResults"

    message <- paste0(" test.mergePermutationAndObservation_good_01() ",
                      "- Valid parameters for mergePermutationAndObservation() ",
                      "did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}

###########################################################
## plotGraph() function
###########################################################

# Test result when all parameters are good
test.plotGraph_good_01 <- function() {

    g <- extractInfo(allResults = methylInheritanceResults,
                                type = "sites", inter="i2", 1)

    obs <- plotGraph(g)

    message <- paste0(" test.plotGraph_good_01() ",
                      "- Valid parameters for plotGraph did not generated expected results.")

    checkTrue("gtable" %in% class(obs), msg = message)
    checkEquals(class(obs[[1]]), "list", msg = message)
    checkEquals(class(obs[[2]]), "data.frame", msg = message)
}


###########################################################
## plotConvergenceGraph() function
###########################################################

# Test result when all parameters are good
test.plotConvergenceGraph_good_01 <- function() {

    filesDir <- system.file("extdata", "TEST", package="methylInheritance")

    ##Extract convergenc information for F1 and F2 and F3
    data <- loadConvergenceData(analysisResultsDir = filesDir,
            permutationResultsDir = filesDir, type = "sites", inter = "iAll",
            position = 1, by = 1)

    ## Create convergence graph
    obs <- plotConvergenceGraph(data)

    message <- paste0(" test.plotConvergenceGraph_good_01() ",
                      "- Valid parameters for plotGraph did not generated expected results.")

    checkTrue("ggplot" %in% class(obs), msg = message)
    checkTrue(is(obs[[1]], "data.frame"), msg = message)
    checkTrue(is(obs[[2]], "list"), msg = message)

}


###########################################################
## loadConvergenceData() function
###########################################################

# Test result when all parameters are good
test.loadConvergenceData_good_01 <- function() {

    filesDir <- system.file("extdata", "TEST", package="methylInheritance")
    ##Extract convergence information for F1 and F2 and F3
    data <- loadConvergenceData(analysisResultsDir = filesDir,
                                permutationResultsDir = filesDir,
                                type = "sites", inter = "iAll",
                                position = 1, by = 1)

    expected <- data.frame(NBR_PERMUTATIONS=c(1,1,2,2,3,3),
                    ELEMENT=rep("SITES", 6), ANALYSIS=rep("iAll", 6),
                    POSITION=rep(1, 6), TYPE=rep(c("HYPER", "HYPO"), 3),
                    SIGNIFICANT_LEVEL=c(1.00000000, 1.0000000, 1.0000000000,
                        0.666666666666666666, 1.000000000, 0.50000000000000))

    message <- paste0(" test.loadConvergenceData_good_01() ",
                      "- Valid parameters for loadConvergenceData did not generated expected results.")

    checkEquals(class(data), "data.frame", msg = message)
    checkEquals(data , expected, msg = message)
}


# Test result when all parameters are good
test.loadConvergenceData_two_different_directories <- function() {

    permutationDir <- system.file("extdata", "TEST_01/permutations",
                                  package="methylInheritance")

    observationDir <- system.file("extdata", "TEST_01/observations",
                                  package="methylInheritance")

    ##Extract convergence information for F1 and F2 and F3
    data <- loadConvergenceData(analysisResultsDir = observationDir,
                                permutationResultsDir = permutationDir,
                                type = "sites", inter = "iAll",
                                position = 1, by = 1)

    expected <- data.frame(NBR_PERMUTATIONS=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10),
                           ELEMENT=rep("SITES", 20), ANALYSIS=rep("iAll", 20),
                           POSITION=rep(1, 20), TYPE=rep(c("HYPER", "HYPO"), 10),
                           SIGNIFICANT_LEVEL=c(1.00000000000000000000, 0.50000000000000000000,
                                               0.66666666666666662966, 0.66666666666666662966,
                                               0.75000000000000000000, 0.50000000000000000000,
                                               0.80000000000000004441, 0.40000000000000002220,
                                               0.83333333333333337034, 0.33333333333333331483,
                                               0.71428571428571430157, 0.42857142857142854764,
                                               0.75000000000000000000, 0.37500000000000000000,
                                               0.77777777777777779011, 0.33333333333333331483,
                                               0.80000000000000004441, 0.29999999999999998890,
                                               0.72727272727272729291, 0.27272727272727270709))

    message <- paste0(" test.loadConvergenceData_two_different_directories() ",
                      "- Using two different directories for observation and permutation with loadConvergenceData did not generated expected results.")

    checkTrue(is(data, "data.frame"), msg = message)
    checkEquals(data , expected, msg = message)
}


# Test result when all parameters are good
test.loadConvergenceData_two_different_directories_by_5 <- function() {

    permutationDir <- system.file("extdata", "TEST_01/permutations",
                                  package="methylInheritance")

    observationDir <- system.file("extdata", "TEST_01/observations",
                                  package="methylInheritance")

    ##Extract convergence information for F1 and F2 and F3
    data <- loadConvergenceData(analysisResultsDir = observationDir,
                                permutationResultsDir = permutationDir,
                                type = "sites", inter = "iAll",
                                position = 1, by = 5)

    expected <- data.frame(NBR_PERMUTATIONS=c(5,5,10,10),
                           ELEMENT=rep("SITES", 4), ANALYSIS=rep("iAll", 4),
                           POSITION=rep(1, 4), TYPE=rep(c("HYPER", "HYPO"), 2),
                           SIGNIFICANT_LEVEL=c(0.83333333333333337034, 0.33333333333333331483,
                                               0.72727272727272729291, 0.27272727272727270709))

    message <- paste0(" test.loadConvergenceData_two_different_directories() ",
                      "- Using two different directories for observation and permutation with loadConvergenceData did not generated expected results.")

    checkEquals(class(data), "data.frame", msg = message)
    checkEquals(data , expected, msg = message)
}
