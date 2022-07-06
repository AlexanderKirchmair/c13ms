


test_that("correctIso works", {

  suppressMessages({
    res <- replicate(5, correctIso(exampleTracerExperiment(seed = sample(1:10^3, 1)), assay = "raw"))
  })

  corr <- sapply(res, function(x) "corr" %in% names(x@isoAssays) )
  expect_true(all(corr))
})




test_that("correctIso functions", {


  C13 <- exampleTracerExperiment(nsamples = 10, nmets = 100, seed = sample(1:10^3, 1))
  assaydata <- .getAssays(C13, assay = 1, type = "iso")
  molecules <- metData(C13)[,c("Molecule", "MSion")]


  MolTable <- .makeMoleculeTable(molecules = molecules, Molecule = rlang::sym("Molecule"), MSion = rlang::sym("MSion"), mode = "negative", tracer = "C")
  ElemTable <- .makeElementsTable(df = MolTable, isocomp = isocomp, tracer = "C", tracer_purity = NA)

  res <- isoCorr(assaydata, molecules, Molecule = Molecule, MSion = MSion, highres = TRUE, tmpdir = NULL, isocomp = isocomp, allresults = FALSE)


  expect_setequal(MolTable$Molecule, metnames(C13))
  expect_true(all(dim(res) == dim(assaydata)))







})




test_that("Element compositions download", {

  isocomp2 <- getElementCompositions()
  expect_equal(isocomp2, c13ms::isocomp)

})

















