




test_that("TracerExperiment operations work", {

  n <- 10
  m <- 15
  C13 <- exampleTracerExperiment(nsamples = n, nmets = m)

  expect_equal(length(colnames(C13)), n)
  expect_equal(nrow(colData(C13)), n)
  expect_equal(length(C13$group), n)

  expect_equal(length(metnames(C13)), m)
  expect_setequal(unique(gsub("_m\\d", "", rownames(C13))), metnames(C13))
  expect_setequal(rownames(metData(C13)), metnames(C13))
  expect_setequal(rownames(isoData(C13)), rownames(C13))

  expect_s4_class(C13, "TracerExperiment")
  expect_s4_class(subset(C13, group == "A"), "TracerExperiment")

  expect_length(split(C13, by = ~ group), 2)

})






test_that("TracerExperiment operations work", {

  C13 <- exampleTracerExperiment()
  assay(C13, "mid") <- MID(C13, assay = "raw")
  assay(C13, "frac", type = "met") <- isoEnrichment(C13, assay = "raw")

  expect_gte(min(assay(C13, "mid")), 0)
  expect_lte(max(assay(C13, "mid")), 1)

  expect_gte(min(assay(C13, "frac", type = "met")), 0)
  expect_lte(max(assay(C13, "frac", type = "met")), 1)

  expect_equal(sum(sumMets(C13, assay = "mid") - 1), 0)

})
















