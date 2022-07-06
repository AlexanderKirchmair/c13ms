




test_that("differential testing recovers expected metabolites", {

  # The first metabolite has some differences added between the two sample groups:
  C13 <- exampleTracerExperiment(nsamples = 10, nmets = 12, seed = sample(1:10^3, 1))
  met <- metnames(C13)[1]
  iso <- isoData(C13) %>% subset(metabolite == met) %>% rownames()

  C13 %<>% impute(assay = "raw")
  C13 %<>% correctIso(assay = "imp")
  C13 %<>% normalize(method = ~ COLSUM, assay = "corr")

  assay(C13, "mid") <- MID(C13)
  assay(C13, "frac", type = "met") <- isoEnrichment(C13)
  assay(C13, "norm", type = "met") <- sumMets(C13)

  contrasts <- list(groupBvsA = list("group" = c("B", "A")))

  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "met", assay = "norm", method = "ttest")
  res_met_norm <- results(C13, "met", "norm", "ttest") %>% subset(padj <= 0.05) %>% rownames()

  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "met", assay = "frac", method = "beta")
  res_met_frac <- results(C13, "met", "frac", "beta") %>% subset(padj <= 0.05) %>% rownames()

  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "iso",  assay = "mid", method = "beta")
  res_iso_mid <- results(C13, "iso", "mid", "beta") %>% subset(padj <= 0.05) %>% rownames()

  # Test if the expected metabolites are recovered:
  expect_true(met %in% res_met_norm[1:min(length(res_met_norm), 3)])
  expect_true(met %in% res_met_norm[1:min(length(res_met_frac), 3)])
  expect_gte(mean(iso %in% res_iso_mid), 0.8)
})














