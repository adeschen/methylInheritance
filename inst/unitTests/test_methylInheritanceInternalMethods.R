###################################################
# Created by Astrid Deschenes
# 2017-01-09
###################################################

###################################################
## Test the methylInheritanceInternalMethods functions
###################################################


library("GenomicRanges")

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

test.validateExtractInfo_position_too_high <- function() {
    obs <- tryCatch(methylInheritance:::validateExtractInfo(
        allResults = methylInheritanceResults, type = "sites",
        inter = "i2", position = 4),
        error=conditionMessage)

    exp <- "position must correspond to a valid entry in the \"allResults$OBSERVATION[[SITES]][[i2]]\""

    message <- paste0("test.validateExtractInfo_position_too_high() - ",
                      "Too high position value did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateExtractInfo_allResults_no_permutation <- function() {
    g<-list()
    g[["OBSERVATION"]] <- methylInheritanceResults$OBSERVATION
    class(g)<-"methylInheritanceAllResults"

    obs <- tryCatch(methylInheritance:::validateExtractInfo(
        allResults = g, type = "sites",
        inter = "i2", position = 1),
        error=conditionMessage)

    exp <- "allResults must have an element called \"PERMUTATION\""

    message <- paste0("test.validateExtractInfo_allResults_no_permutation() - ",
                      "allResult without permutation value did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateExtractInfo_allResults_no_observation <- function() {
    g<-list()
    g[["PERMUTATION"]] <- methylInheritanceResults$OBSERVATION
    class(g)<-"methylInheritanceAllResults"

    obs <- tryCatch(methylInheritance:::validateExtractInfo(
        allResults = g, type = "sites",
        inter = "i2", position = 1),
        error=conditionMessage)

    exp <- "allResults must have an element called \"OBSERVATION\""

    message <- paste0("test.validateExtractInfo_allResults_no_observation() - ",
                      "allResult without observation value did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateExtractInfo_allResults_not_list <- function() {
    g<-vector()
    class(g)<-"methylInheritanceAllResults"

    obs <- tryCatch(methylInheritance:::validateExtractInfo(
        allResults = g, type = "sites",
        inter = "i2", position = 1),
        error=conditionMessage)

    exp <- "allResults must be a list"

    message <- paste0("test.validateExtractInfo_allResults_not_list() - ",
                      "allResult not list did not generated expected message.")

    checkEquals(obs, exp, message)
}

###########################################################
## validateLoadConvergenceData() function
###########################################################

test.validateLoadConvergenceData_integer_analysisDir <- function() {
    obs <- tryCatch(methylInheritance:::validateLoadConvergenceData(analysisResultsDir = 33,
            permutationResultsDir = "./", position = 1, by = 2),
            error=conditionMessage)

    exp <- "analysisResultsDir must be a character string"

    message <- paste0("test.validateLoadConvergenceData_integer_analysisDir() - ",
                      "Integer as analysis directory did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateLoadConvergenceData_integer_permutationDir <- function() {
    obs <- tryCatch(methylInheritance:::validateLoadConvergenceData(analysisResultsDir = "./",
                    permutationResultsDir = 44, position = 1, by = 2),
                    error=conditionMessage)

    exp <- "permutationResultsDir must be a character string"

    message <- paste0("test.validateLoadConvergenceData_integer_permutationDir() - ",
                      "Integer as permutation directory did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateLoadConvergenceData_position_zero <- function() {
    obs <- tryCatch(methylInheritance:::validateLoadConvergenceData(analysisResultsDir = "./",
                        permutationResultsDir = "./", position = 0, by = 2),
                    error=conditionMessage)

    exp <- "position must be a positive integer"

    message <- paste0("test.validateLoadConvergenceData_position_zero() - ",
                      "Zero as position did not generated expected message.")

    checkEquals(obs, exp, message)
}

test.validateLoadConvergenceData_by_zero <- function() {
    obs <- tryCatch(methylInheritance:::validateLoadConvergenceData(analysisResultsDir = "./",
                permutationResultsDir = "./", position = 1, by = 0),
                error=conditionMessage)

    exp <- "by must be a positive integer"

    message <- paste0("test.validateLoadConvergenceData_by_zero() - ",
                      "Zero as by did not generated expected message.")

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

###########################################################
## readInterGenerationResults() function
###########################################################

test.readInterGenerationResults_good_01 <- function() {

    filesDir <- system.file("extdata", "TEST", package="methylInheritance")

    ## Read DMS intergenerational results for the observed data
    obs <- methylInheritance:::readInterGenerationResults(outputDir =
            paste0(filesDir, "/"), 0, "sites")

    message <- paste0("test.readInterGenerationResults_good_01() ",
                      "- Function did not return expected values")

    iAll_1 <- GenomicRanges::GRanges(seqnames = rep("S", 4),
                           ranges = IRanges::IRanges(start = c(30222185, 15048832, 22963194, 23499048),
                                                     end = c(30222185, 15048832, 22963194, 23499048)),
                           strand = rep("+", 4),
                           typeDiff = c(1, rep(-1, 3)))

    i2_2 <- GenomicRanges::GRanges(seqnames = rep("S", 22),
                                   ranges = IRanges::IRanges(start = c(97481, 572272, 3281006, 11121503, 19260516,
                                                                       19445653, 22874019, 27232572, 30222185, 35929511,
                                                                       6085769,  8045209, 10355001, 11147625, 15048832,
                                                                       15438496, 22745004, 22899924, 22963194, 23499048,
                                                                       28622167, 34611139),
                                                             end = c(97481, 572272, 3281006, 11121503, 19260516,
                                                                     19445653, 22874019, 27232572, 30222185, 35929511,
                                                                     6085769,  8045209, 10355001, 11147625, 15048832,
                                                                     15438496, 22745004, 22899924, 22963194, 23499048,
                                                                     28622167, 34611139)),
                                   strand = rep("+", 22),
                                   typeDiff = c(rep(1, 10), rep(-1, 12)))

    i2_1 <- GenomicRanges::GRanges(seqnames = rep("S", 36),
                                   ranges = IRanges::IRanges(start = c(570115, 2573229, 5063112, 8247138, 8765627, 8791494,
                                                                       9955639, 11667875, 19095767, 26225126, 26798489, 27089337,
                                                                       27188724, 27236909, 27421271, 30222185, 30786437, 31364173,
                                                                       31396094, 33611091, 33886929,  3139258, 14391040, 15048832,
                                                                       15438613, 16630377, 17795264, 18396852, 22963194, 23499048,
                                                                       23499106, 23499111, 27019812, 30204193, 30746773, 35827911),
                                                             end = c(570115, 2573229, 5063112, 8247138, 8765627, 8791494,
                                                                     9955639, 11667875, 19095767, 26225126, 26798489, 27089337,
                                                                     27188724, 27236909, 27421271, 30222185, 30786437, 31364173,
                                                                     31396094, 33611091, 33886929,  3139258, 14391040, 15048832,
                                                                     15438613, 16630377, 17795264, 18396852, 22963194, 23499048,
                                                                     23499106, 23499111, 27019812, 30204193, 30746773, 35827911)),
                                   strand = rep("+", 36),
                                   typeDiff = c(rep(1, 21), rep(-1, 15)))

    exp <- list("i2" = list(i2_1, i2_2), "iAll" = list(iAll_1))

    checkEquals(obs, exp, msg = message)
}

test.readInterGenerationResults_good_02 <- function() {

    filesDir <- system.file("extdata", "TEST", package="methylInheritance")

    ## Read DMS intergenerational results for the observed data
    obs <- methylInheritance:::readInterGenerationResults(outputDir =
                                                              paste0(filesDir, "/"), 1, "sites")

    message <- paste0("test.readInterGenerationResults_good_02() ",
                      "- Function did not return expected values")

    iAll_1 <- GenomicRanges::GRanges(seqnames = rep("S", 5),
                                     ranges = IRanges::IRanges(start = c(17191066, 17424070, 1130005, 22786110, 26615081),
                                                               end = c(17191066, 17424070, 1130005, 22786110, 26615081)),
                                     strand = rep("+", 5),
                                     typeDiff = c(1, 1, rep(-1, 3)))

    i2_2 <- GenomicRanges::GRanges(seqnames = rep("S", 18),
                                   ranges = IRanges::IRanges(start = c(5507460,  8045221, 17191066, 17424070, 32510128,
                                                                       32911751, 33071323,  1130005,  5826332,  8045207,
                                                                       17804285, 18396726, 19869520, 21672484, 22786110,
                                                                       26615081, 32216298, 32216478),
                                                             end = c(5507460,  8045221, 17191066, 17424070, 32510128,
                                                                     32911751, 33071323,  1130005,  5826332,  8045207,
                                                                     17804285, 18396726, 19869520, 21672484, 22786110,
                                                                     26615081, 32216298, 32216478)),
                                   strand = rep("+", 18),
                                   typeDiff = c(rep(1, 7), rep(-1, 11)))

    i2_1 <- GenomicRanges::GRanges(seqnames = rep("S", 15),
                                   ranges = IRanges::IRanges(start = c(3401344, 17191066, 17424070, 24135743,
                                                                       27019812,  1130005,  1345775,  2573229,
                                                                       6717075, 18926407, 19260516, 22786110,
                                                                       23655774, 26615081, 33278578),
                                                             end = c(3401344, 17191066, 17424070, 24135743,
                                                                     27019812,  1130005,  1345775,  2573229,
                                                                     6717075, 18926407, 19260516, 22786110,
                                                                     23655774, 26615081, 33278578)),
                                   strand = rep("+", 15),
                                   typeDiff = c(rep(1, 5), rep(-1, 10)))

    exp <- list("i2" = list(i2_1, i2_2), "iAll" = list(iAll_1))

    checkEquals(obs, exp, msg = message)
}


###########################################################
## formatInputMethylData() function
###########################################################

test.formatInputMethylData_good_01 <- function() {

    initGR_01 <- list()
    initGR_01[[1]] <- new("methylRaw", data.frame(chr = rep("S", 3),
                                                start = c(1005, 1011, 1017), end = c(1005, 1011, 1017),
                                                strand = strand(rep("+", 3)), coverage = c(71, 90, 95),
                                                numCs = c(0, 1, 1), numTs = c(71, 89, 94)), sample.id = "F1_1_C",
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    initGR_01[[2]] <- new("methylRaw", data.frame(chr = rep("S", 3),
                                                start = c(1006, 1011, 1017), end = c(1006, 1011, 1017),
                                                strand = strand(rep("+", 3)), coverage = c(93, 92, 93),
                                                numCs = c(1, 4, 0), numTs = c(92, 88, 93)), sample.id = "F1_2_C",
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')

    initGR_02 <- list()
    initGR_02[[1]] <- new("methylRaw", data.frame(chr = rep("S", 3),
                                                start = c(1007, 1011, 1017), end = c(1007, 1011, 1017),
                                                strand = strand(rep("+", 3)), coverage = c(90, 85, 79),
                                                numCs = c(0, 0, 1), numTs = c(90, 85, 78)), sample.id = "F2_1_C",
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    initGR_02[[2]] <- new("methylRaw", data.frame(chr = rep("S", 3),
                                                start = c(1008, 1011, 1017), end = c(1008, 1011, 1017),
                                                strand = strand(rep("+", 3)), coverage = c(73, 93, 78),
                                                numCs = c(0, 2, 0), numTs = c(73, 91, 78)), sample.id = "F2_2_C",
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')

    initGR_03 <- list()
    initGR_03[[1]] <- new("methylRaw", data.frame(chr = rep("S", 3),
                                                start = c(1009, 1011, 1017), end = c(1009, 1011, 1017),
                                                strand = strand(rep("+", 3)), coverage = c(80, 73, 84),
                                                numCs = c(0, 0, 1), numTs = c(80, 73, 83)), sample.id = "F3_1_C",
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    initGR_03[[2]] <- new("methylRaw", data.frame(chr = rep("S", 3),
                                                start = c(1010, 1011, 1017), end = c(1010, 1011, 1017),
                                                strand = strand(rep("+", 3)), coverage = c(77, 80, 94),
                                                numCs = c(0, 2, 2), numTs = c(77, 78, 92)), sample.id = "F3_2_C",
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')

    initGR <- list()
    initGR[[1]] <- new("methylRawList", initGR_01, treatment = c(0, 1))
    initGR[[2]] <- new("methylRawList", initGR_02, treatment = c(0, 1))
    initGR[[3]] <- new("methylRawList", initGR_03, treatment = c(0, 1))

    set.seed(20011)
    obs <- methylInheritance:::formatInputMethylData(initGR)

    expGR <- list()
    expGR[[1]] <- new("methylRawList", list(initGR_03[[2]], initGR_02[[2]]), treatment = c(0, 1))
    expGR[[2]] <- new("methylRawList", list(initGR_01[[1]], initGR_03[[1]]), treatment = c(0, 1))
    expGR[[3]] <- new("methylRawList", list(initGR_01[[2]], initGR_02[[1]]), treatment = c(0, 1))

    message <- paste0("test.formatInputMethylData_good_01() ",
                      "- Function did not return expected values")

    checkEquals(obs, expGR, message)
}
