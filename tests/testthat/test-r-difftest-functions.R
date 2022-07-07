




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




test_that("differential testing output works", {

  C13 <- exampleTracerExperiment(nsamples = 10, nmets = 5, seed = sample(1:10^3, 1))
  assay(C13, "mid") <- MID(C13, assay = "raw")
  assay(C13, "frac", type = "met") <- isoEnrichment(C13, assay = "raw")
  assay(C13, "norm", type = "met") <- sumMets(C13, assay = "raw")

  colData(C13)$group <- factor(colData(C13)$group)
  colData(C13)$donor <- factor(gsub("A|B", "", colData(C13)$name))

  contrasts <- list(groupAvsB = list("group" = c("A", "B")),
                    groupBvsA = list("group" = c("B", "A")))

  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "met", assay = "norm", method = "ttest")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group + donor, type = "met", assay = "norm", method = "lm")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "met", assay = "norm", method = "lm")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group + (1|donor), type = "met", assay = "norm", method = "limma")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "met", assay = "norm", method = "limma")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, random = ~ 1 | donor, type = "met", assay = "norm", method = "lmm")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group + (1|donor), type = "met", assay = "norm", method = "dream")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group + donor, type = "met", assay = "frac", method = "beta")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group + donor, type = "met", assay = "frac", method = "betareg")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group + (1|donor), type = "iso",  assay = "mid", method = "beta")
  C13 %<>% diffTest(contrasts = contrasts, formula = ~ group + (1|donor), type = "iso",  assay = "mid", method = "limma")

  diffres <- results(C13)

  expect_true(all(abs(diffres$met$norm$ttest$groupAvsB$lfc + diffres$met$norm$ttest$groupBvsA$lfc) < 10^-10)) # ttest
  expect_true(all(abs(diffres$met$norm$lm$groupAvsB$lfc + diffres$met$norm$lm$groupBvsA$lfc) < 10^-10)) # lm
  expect_setequal(diffres$met$norm$limma$groupAvsB$lfc, -diffres$met$norm$limma$groupBvsA$lfc) # limma
  expect_setequal(diffres$met$norm$lmm$groupAvsB$lfc, -diffres$met$norm$lmm$groupBvsA$lfc) # lmm
  expect_setequal(diffres$met$norm$dream$groupAvsB$lfc, -diffres$met$norm$dream$groupBvsA$lfc) # dream
  expect_setequal(diffres$iso$mid$beta$groupAvsB$diff, -diffres$iso$mid$beta$groupBvsA$diff) # beta

  expect_gt(cor(diffres$met$norm$ttest$groupAvsB$lfc, diffres$met$norm$lm$groupAvsB$lfc), 0.8)
  expect_gt(cor(diffres$met$norm$ttest$groupAvsB$lfc, diffres$met$norm$limma$groupAvsB$lfc), 0.8)
  expect_gt(cor(diffres$met$norm$ttest$groupAvsB$lfc, diffres$met$norm$lmm$groupAvsB$lfc), 0.8)
  expect_gt(cor(diffres$met$norm$ttest$groupAvsB$lfc, diffres$met$norm$dream$groupAvsB$lfc), 0.8)
  expect_gt(cor(diffres$met$frac$beta$groupAvsB$diff, diffres$met$frac$betareg$groupAvsB$diff, use = "pairwise.complete.obs"), 0.8)
  expect_gt(cor(diffres$iso$mid$beta$groupAvsB$lfc, diffres$iso$mid$limma$groupAvsB$lfc, use = "pairwise.complete.obs"), 0)

})


test_that("contrasts work", {

  C13 <- exampleTracerExperiment(nsamples = 15, nmets = 5)
  colData(C13)$group <- factor(rep(LETTERS[1:3], each = 5))
  colData(C13)$name <- paste0(colData(C13)$group, 1:5)
  colData(C13)$donor <- gsub("A|B|C", "", colData(C13)$name)
  colData(C13)$class <- ifelse(colData(C13)$donor == "5", "5", "not5")
  assay(C13, "norm", type = "met") <- sumMets(C13, assay = "raw")
  data <- assay(C13, "norm", type = "met")
  data <- as.data.frame(t(data))
  design <- colData(C13)

  contrasts <- list(groupBvsA = list(group = c("B", "A"), class = "not5"),
                    groupCvsA = list(group = c("C", "A"), class = "not5"),
                    groupBvsC = list(group = c("B", "C"), class = "5"))

  subsets <- getSubsets(contrasts)
  subset_data <- lapply(subsets, function(x) getDataSubset(data, design, contrast = x) )
  samples <- getContrastSamples(design, contrasts$groupBvsA)
  expect_setequal(names(subsets), c("class==5", "class==not5"))
  expect_setequal(rownames(subset_data$`class==5`), rownames(subset(design, class == "5")))
  expect_setequal(rownames(subset_data$`class==not5`), rownames(subset(design, class == "not5")))
  expect_setequal(names(samples), rownames(subset(design, class == "not5" & group != "C")))

})



