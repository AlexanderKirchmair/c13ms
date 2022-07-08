
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Analysis of 13C-labelled metabolomics data

<!-- badges: start -->
<!-- badges: end -->

R package for the statistical and functional analysis of 13C labelling
data.

## Installation

``` r
devtools::install_github("AlexanderKirchmair/c13ms")
```

## Introduction

Isotopologue abundances are required as a main input. The workflow is
centered around the ‘TracerExperiment’ class object, which acts as an
integrative storage of isotopologue abundances, isotopologue and
metabolite annotations, sample annotations and other metadata.

Here are some basic functionalities:

``` r
C13 <- exampleTracerExperiment(nsamples = 6, nmets = 5)

isoData(C13)
metData(C13)

colnames(C13)
metnames(C13)
rownames(C13)

colData(C13)
C13$group

C13 <-  makeTracerExperiment(cbind(isoData(C13), assay(C13, "raw")), metData = metData(C13), colData = colData(C13))

C13 %>% subset(group == "A")
C13groups <- split(C13, by = ~ group)
C13 <- with(C13groups, A + B)

assay(C13, "raw")
sumMets(C13, assay = "raw")
```

## Workflow

A typical analysis workflow may look like as demonstrated below.

Pre-processing (imputation, natural isotope abundance correction,
normalization, …):

``` r
C13 %<>% impute(assay = "raw")
C13 %<>% correctIso(assay = "imp")
C13 %<>% normalize(method = ~ COLSUM, assay = "corr")
```

Calculation of relative mass isotopomer distributions (MID), fractional
enrichment and summarization of isotopologues to metabolite levels:

``` r
assay(C13, "mid") <- MID(C13)
assay(C13, "frac", type = "met") <- isoEnrichment(C13)
assay(C13, "norm", type = "met") <- sumMets(C13)
```

Statistical testing for differences in abundances and labelling:

``` r
contrasts <- list(groupBvsA = list("group" = c("B", "A")))

C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "met", assay = "norm", method = "ttest")
C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "met", assay = "frac", method = "beta")
C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "iso",  assay = "mid", method = "beta")

results(C13, "iso", "mid", "beta") %>% head(10)
```

## Visualization

``` r
isoplot(C13, mets = metnames(C13)[1], cumulative = T)
```

<img src="man/figures/README-unnamed-chunk-7-1.svg" width="100%" />
