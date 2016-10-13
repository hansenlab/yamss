library(yamss)
library(mtbls2)

filepath <- file.path(find.package("mtbls2"), "mzData")
files <- list.files(filepath, pattern = "MSpos-Ex1.*Ag-[12]", recursive = TRUE, full.names = TRUE)
classes <- rep(c("wild-type", "mutant"), each = 2)

colData <- DataFrame(sampClasses = classes, ionmode = "pos")
cmsRawExample <- readMSdata(files = files, colData = colData, mzsubset = c(500,502), verbose = TRUE)
format(object.size(cmsRawExample), units = "Mb")
save(cmsRawExample, file = "../../data/cmsRawExample.rda")


cmsProcExample <- bakedpi(cmsRawExample, dbandwidth = c(0.01, 10), dgridstep = c(0.01, 1),
                   outfileDens = NULL, dortalign = TRUE, verbose = TRUE)
format(object.size(cmsProcExample), units = "Mb")


cutoff <- densityQuantiles(cmsProcExample)["99.9%"]
cmsSliceExample <- slicepi(cmsProcExample, cutoff = cutoff, verbose = TRUE)