library(yamss)
load(file.path(path.package("yamss"), "scripts", "digests.rda"))
data(cmsRawExample)

context("Testing bakedpi")
background <- yamss:::backgroundCorrection(cmsRawExample)
test_that("Background correction", {
    expect_equal(yamss:::.digestDataTableRaw(yamss:::.rawDT(background)), digests$backgroundRaw)
    expect_equal(yamss:::.digestDataTableBG(yamss:::.bgcorrDT(background)), digests$backgroundCorr)
    })

rtAlign <- yamss:::rtAlignment(background)
test_that("Retention time alignment", {
    expect_equal(yamss:::.digestDataTableRaw(yamss:::.rawDT(rtAlign)), digests$rtRaw)
    expect_equal(yamss:::.digestDataTableBG(yamss:::.bgcorrDT(rtAlign)), digests$rtCorr)
})

baked <- bakedpi(cmsRawExample, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1),
                 dortalign = TRUE, mzsubset = c(500, 510))
test_that("bakedpi", {
    expect_equal(yamss:::.digestDataTableRaw(yamss:::.rawDT(baked)), digests$bakedRaw)
    expect_equal(yamss:::.digestDataTableBG(yamss:::.bgcorrDT(baked)), digests$bakedCorr)
})

cutoff <- densityQuantiles(baked)["99.9%"]
sliced <- slicepi(baked, cutoff = cutoff, verbose = TRUE)
test_that("slicepi", {
    expect_equal(yamss:::.digestPeakBounds(peakBounds(sliced)), digests$slicedPeakBounds)
    expect_equal(yamss:::.digestPeakQuants(peakQuants(sliced)), digests$slicedPeakQuants)
})
