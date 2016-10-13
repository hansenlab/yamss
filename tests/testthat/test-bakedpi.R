library(yamss)
context("Testing Background Correction")
load(file.path(path.package("yamss"), "scripts", "digests.rda"))
data(cmsRawExample)
background <- yamss:::backgroundCorrection(cmsRawExample)
test_that("Background correction", {
    expect_equal(yamss:::.digestDataTableRaw(yamss:::.rawDT(background)), digests$backgroundRaw)
    expect_equal(yamss:::.digestDataTableBG(yamss:::.bgcorrDTDT(background)), digests$backgroundCorr)
    })
