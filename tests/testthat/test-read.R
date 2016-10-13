library(yamss)
context("Reading raw data")
if (require(mtbls2)) {
	filepath <- file.path(find.package("mtbls2"), "mzData")
	files <- list.files(filepath, pattern = "MSpos-Ex1.*Ag-[12]", recursive = TRUE, full.names = TRUE)
	classes <- rep(c("wild-type", "mutant"), each = 2)
	colData <- DataFrame(sampClasses = classes, ionmode = "pos")
	cmsRaw <- readMSdata(files = files, colData = colData, mzsubset = c(500,520), verbose = FALSE)
	rawDT <- .rawDT(cmsRaw)
	test_that("rawDT has correct properties", {
		expect_equal(length(unique(rawDT[,sample])), length(files))
		expect_true(all(is.integer(rawDT[,mz])))
		expect_true(all(rawDT[,intensity] > 0))
		expect_true(all(rawDT[,scan] > 0))
	})
}