## Creating test digest objects
library(yamss)
data(cmsRawExample)
digests <- list()
background <- yamss:::backgroundCorrection(cmsRawExample)
digests$backgroundRaw <- yamss:::.digestDataTableRaw(yamss:::.rawDT(background))
digests$backgroundCorr <- yamss:::.digestDataTableBG(yamss:::.bgcorrDT(background))
save(digests, file = "digests.rda")                                                  
sessionInfo()

