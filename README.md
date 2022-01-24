
<!-- README.md is generated from README.Rmd. Please edit that file -->

# c13ms

<!-- badges: start -->
<!-- badges: end -->

Analysis of isotope-labelled metabolomics data

``` r
library(c13ms)

Attaching package: 'c13ms'
The following object is masked from 'package:base':

    array
library(magrittr)
```

### Introduction

``` r
C13 <- exampleTracerExperiment()

colnames(C13)
metnames(C13)
rownames(C13)

colData(C13)
C13$group

# C13 %>% subset(group == "A")
# C13groups <- split(C13, by = ~ group)
# C13 <- with(C13groups, A + B)

metData(C13)
isoData(C13)

assay(C13, "raw")
```

## Workflow

A typical analysis workflow may look like as demonstrated below:

``` r
C13 %<>% impute(assay = "raw")
C13 %<>% correctIso()
C13 %<>% normalize(method = ~ COLSUM)
Normalization using normCOLSUM 

assay(C13, "mid") <- MID(C13)
assay(C13, "frac", type = "met") <- isoEnrichment(C13)
assay(C13, "norm", type = "met") <- sumMets(C13)

contrasts <- list(groupBvsA = list("group" = c("B", "A")))

C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "met", method = "limma", assay = "norm")
C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "met", method = "beta", assay = "frac")
C13 %<>% diffTest(contrasts = contrasts, formula = ~ group, type = "iso", method = "beta", assay = "mid")

results(C13, "iso", "beta", "mid") %>% head()
NULL
```

## Visualization

``` r
isoplot(C13, mets = metnames(C13)[1])
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
