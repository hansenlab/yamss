## Creating test digest objects
library(yamss)
data(cmsRawExample)
digests <- list()

background <- yamss:::backgroundCorrection(cmsRawExample)
digests$backgroundRaw <- yamss:::.digestDataTableRaw(yamss:::.rawDT(background))
digests$backgroundCorr <- yamss:::.digestDataTableBG(yamss:::.bgcorrDT(background))

rtAlign <- yamss:::rtAlignment(background)
digests$rtRaw <- yamss:::.digestDataTableRaw(yamss:::.rawDT(rtAlign))
digests$rtCorr <- yamss:::.digestDataTableBG(yamss:::.bgcorrDT(rtAlign))

baked <- bakedpi(cmsRawExample, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1),
                 dortalign = TRUE, mzsubset = c(500, 510))
digests$bakedRaw <- yamss:::.digestDataTableRaw(yamss:::.rawDT(baked))
digests$bakedCorr <- yamss:::.digestDataTableBG(yamss:::.bgcorrDT(baked))

cutoff <- densityQuantiles(baked)["99.9%"]
sliced <- slicepi(baked, cutoff = cutoff, verbose = TRUE)
digests$slicedRowData <- yamss:::.digestRowData(rowData(sliced))
digests$slicedPeakQuants <- yamss:::.digestPeakQuants(peakQuants(sliced))

save(digests, file = "digests.rda")                                                  
sessionInfo()

