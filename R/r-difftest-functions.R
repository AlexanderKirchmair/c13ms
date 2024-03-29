
# DIFFERENTIAL ABUNDANCE TESTING FUNCTIONS

#' Differential abundance testing
#'
#' Various methods for differential abundance testing of isotope labelling data
#'
#' Metabolite and isotopologue abundances (values from 0 to Inf):
#' Fractional labelling of metabolites and isotopologues (values between 0 and 1):
#'
#' @param TE TE object
#' @param contrasts List of contrasts in the format list(name = c("B", "A"), subclass = "celltype1")
#' @param formula Formula
#' @param method
#' @param type "iso" or "met"
#' @param assay assay name
#' @param logged are data already log-transformed?
#' @param p.adj.method
#' @param conf LOQ confidence values
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' exampleTracerExperiment() %>% diffTest(contrasts = list(BvsA = list(group = c("B", "A"))), formula = ~ group, assay = "raw")
diffTest <- function(TE, contrasts, formula = NULL, method = "ttest", type = "iso", assay = "norm", assay_ref = "clean", logged = FALSE, p.adj.method = "fdr", conf = TRUE, ...){


  ### Input ----

  design <- TE@colData


  data <- c13ms:::.getAssays(TE, assay = assay, type = type)


  ### Pre-filtering ----
  # fterms <- labels(terms(formula)) |> sub(pattern = " ", replacement = "") |> sub(pattern = "1|", replacement = "", fixed = TRUE) |> trimws()
  # if (length(fterms) > 0){
  #   split_by <- as.formula(paste0("~", paste(fterms, collapse = " + ")))
  # } else {
  #   split_by <- ~1
  # }
  #
  # data <- clean(TE, assay = assay, type = type, qc_LOQ = qc_LOQ, thres_LOQ = thres_LOQ, min_rep_per_group = 0, split_by = split_by)


  ### Call method functions ----
  # Any function defined as 'test"METHOD"' can be called
  # Results from each method (i.e. a list for each contrast)
  results <- lapply(setNames(method, method), function(testfun){
    testfun <- paste0("test", toupper(testfun))
    do.call(what = testfun, args = c(list("data" = data,
                                          "design" = design,
                                          "formula" = formula,
                                          "contrasts" = contrasts,
                                          "logged" = logged,
                                          "p.adj.method" = p.adj.method, ...)))
  })


  ### Padjust for all contrasts ---
  results <- lapply(results, function(tmpres){ # for each method
    respadj <- sapply(tmpres, function(tmp) length(tmp$padj_all) ) # for each contrast
    if (all(respadj) == 0){
      tmpres <- lapply(tmpres, data.frame)
      pmat <- sapply(tmpres, function(tmp) tmp$pval )
      padj <- matrix(p.adjust(pmat, method = p.adj.method), ncol = length(tmpres), dimnames = dimnames(pmat))
      for (tmpname in names(tmpres)){ tmpres[[tmpname]]$padj_all <- padj[,tmpname] }
    }
    tmpres
  })




  # ### Add conf values ---
  # if (conf == TRUE & !is.null(TE@qcAssays$conf)){
  #
  #   results <- lapply(results, function(res){
  #     lapply(setNames(names(res), names(res)), function(tmpcontr){
  #
  #       tmpres <- res[[tmpcontr]]
  #       samples <- getContrastSamples(design, contrasts[[tmpcontr]])
  #       confdata <- data.matrix(TE@qcAssays$conf)
  #
  #       if (type == "iso"){
  #         confdata[is.na(data.matrix(data))] <- NA
  #
  #       } else {
  #         tmp <- .getAssays(TE, assay = assay_ref, type = "iso")
  #         confdata[is.na(data.matrix(tmp))] <- NA
  #         confdata <- data.frame(confdata)
  #       }
  #
  #       if (type == "met"){
  #         confdata$metabolite <- TE@isoData$metabolite
  #         confdata <- .sumAssay(confdata, FUN = min, na.rm = TRUE)
  #       }
  #
  #       groups <- unique(samples)
  #       g1 <- rowMeans(data.matrix(confdata[,names(samples)[samples == groups[1]]]))
  #       g2 <- rowMeans(data.matrix(confdata[,names(samples)[samples == groups[2]]]))
  #       conf <- matrixStats::rowMaxs(cbind(g1, g2), na.rm = TRUE)
  #       conf <- setNames(conf, rownames(confdata))
  #       tmpres$conf <- conf[rownames(tmpres)]
  #       tmpres
  #     })
  #   })
  #
  #
  # }



  ### Add to TE results ---
  for (re in names(results)){
    TE@results[[type]][[assay]][[re]][names(results[[re]])] <- results[[re]]
  }

  TE
}




### HELPER FUNCTIONS --------------------------------





#' Internal function to group contrasts that use the same data subset
#'
#' @param contrasts
#'
#' @return
#'
#' @examples
getSubsets <- function(contrasts){

  all_subsets <- lapply(contrasts, function(tmp){
    tmp <- tmp[sapply(tmp, length) == 1]
    unlist(tmp)
  })

  all_names <- lapply(all_subsets, function(tmp){
    if (is.null(tmp)){
      "all"
    } else {
      paste0(names(tmp), "==", tmp)
    }
  })
  all_names <- sapply(all_names, paste, collapse = "&")

  split(contrasts, all_names)
}




#' Internal function to get the data subset for a given contrast
#'
#' @param data
#' @param design
#' @param contrast
#'
#' @return
#'
#' @examples
getDataSubset <- function(data, design, contrast){

  factors <- unlist(lapply(setNames(contrast, NULL), function(tmp) tmp[sapply(tmp, length) == 1]), recursive = FALSE)
  ixl <- lapply(names(factors), function(fac) design[[fac]] == factors[[fac]])
  ix <- Reduce('&', ixl)

  if (is.null(ix)) return(data[rownames(design),,drop = FALSE])

  samples <- rownames(design)[ix]
  data[samples,,drop = FALSE]
}






#' Internal function to get the samples/sample names for a given contrast
#'
#' @param design
#' @param contrast
#' @param paired
#'
#' @return
#'
#' @examples
getContrastSamples <- function(design, contrast, paired = NULL){

  if (!is.null(paired)){
    if (paired == FALSE) paired <- NULL
  }

  to_subset <- sapply(contrast, length) == 1
  to_compare <- sapply(contrast, length) == 2

  # use only samples of given type (defined in design)
  if (any(to_subset)){
    cols <- names(contrast)[to_subset]
    tmp <- sapply(cols, function(col) design[[col]] %in% contrast[[col]] )
    ix_use <- apply(tmp, 1, mean) == 1
  } else {
    ix_use <- TRUE
  }

  # get samples for the comparison
  col <- names(contrast)[which(to_compare)]
  levels <- contrast[[col]]

  comp <- which(.naf(as.character(design[[col]]) == levels[1] & ix_use))
  ref <- which(.naf(as.character(design[[col]]) == levels[2] & ix_use))

  x <- factor(c(rep(levels[1], length(comp)), rep(levels[2], length(ref))))
  x <- relevel(x, ref = levels[2])
  x <- setNames(x, c(rownames(design)[comp], rownames(design)[ref]))

  if (!is.null(paired)){
    p <- setNames(design[[paired]], rownames(design))
    p1 <- p[names(x)[x == levels[1]]]
    p2 <- p[names(x)[x == levels[2]]]
    p2 <- p2[match(p1, p2)]
    x <- x[c(names(p1), names(p2))]
  }

  x
}





### METHODS FUNCTIONS --------------------------------


# Structure of testfun:
# testfun <- function(assaydata, sampledata, formula, contrasts, logged = FALSE, ...){
#
#   # assaydata: data matrix
#   # formula: model formula
#   # contrasts: list of contrasts
#   # logged: whether data are log-transformed or not
#
#
#   return(results)
# }




testTTEST <- function(data, design, formula = NULL, contrasts, logged = FALSE, var.equal = FALSE, p.adj.method = "fdr", paired = FALSE, FUN = NULL,  ...){

  if (paired == TRUE) message("Warning: Paired t-test requires matching non-NA samples .")
  if (paired == TRUE & !is.null(formula)) paired <- labels(terms(formula))

  if (logged == FALSE){
    data <- .logtrans(data, FUN = FUN)
  }

  results <- lapply(contrasts, function(contr){

    x <- getContrastSamples(design, contr, paired)

    data <- cbind(data[,names(x)[x == levels(x)[2]]], data[,names(x)[x == levels(x)[1]]]) ###

    ttest_res <- data.frame(t(apply(data, 1, function(conc){


      if (paired == TRUE){
        ix_na <- setNames(is.na(conc[names(x[x == levels(x)[1]])]) | is.na(conc[names(x[x == levels(x)[2]])]), NULL)
        conc[names(x[x == levels(x)[1]])][ix_na] <- NA
        conc[names(x[x == levels(x)[2]])][ix_na] <- NA
      }

      n <- sapply(unique(x), function(xx) sum(!is.na(conc[xx == x]))) # check if both groups have enough data

      if (all(n > 1)) {
        res <- t.test(conc ~ x, var.equal = var.equal, paired = paired, na.action = "na.omit")
        if (paired == FALSE){
          # ratio of means
          est <- setNames(res$estimate, gsub("mean in group ", "", names(res$estimate)))
          fc <- as.numeric(est[levels(x)[2]] / est[levels(x)[1]])
          diff <- as.numeric(est[levels(x)[2]] - est[levels(x)[1]])
        } else {
          # mean of ratios
          fc <- mean(conc[x == levels(x)[2]] / conc[x == levels(x)[1]], na.rm = TRUE)
          diff <- mean(conc[x == levels(x)[2]] - conc[x == levels(x)[1]], na.rm = TRUE)
        }


        c("pval" = res$p.value, "lfc" = fc, "diff" = diff)
      } else {
        c("pval" = NA, "lfc" = NA, "diff" = NA)
      }

    })))

    # if (logged == FALSE){ ttest_res$lfc <- log2(ttest_res$lfc) }

    ttest_res$padj <- p.adjust(ttest_res$pval, method = p.adj.method)

    ttest_res
  })


  results
}



# testLM(data = C13@metAssays$frac, design = C13@colData, formula = ~ Celltype + Donor, contrasts = contrasts)
# add missing log-transform in the beginning?
testLM <- function(data, design, formula, contrasts, logged = FALSE, p.adj.method = "fdr", FUN = NULL, ...){

  stopifnot(requireNamespace("contrast", quietly = TRUE))

  ### Input arguments ----

  data <- as.data.frame(t(data))
  data[is.na(data)] <- NA
  formula <- update(formula, conc ~  0 + .)

  if (logged == FALSE){
    data <- .logtrans(data, FUN = FUN)
  }

  ### Contrasts and formula ----

  subsets <- getSubsets(contrasts)


  ### Model fitting and testing ----

  results <- lapply(subsets, function(tmpsubset){

    subset_data <- getDataSubset(data, design, contrast = tmpsubset)
    subset_design <- design[rownames(subset_data), intersect(all.vars(formula), colnames(design)),drop = FALSE]
    # subset_design <- subset_design[,colnames(subset_design) %in% labels(terms(formula)), drop = FALSE]

    lapply(subset_data, function(conc){

      df <- data.frame(subset_design, conc)
      df <- droplevels(na.omit(df))

      if (all(sapply((df), function(x) length(.unique.na(x))) > 1)){

        lmfit <- lm(formula, data = df)


        sapply(tmpsubset, function(contr){

          if (all(contr[[1]] %in% lmfit$model[[names(contr)[1]]]) & !any(is.na(lmfit$coefficients))){

            all <- lapply(df[,colnames(subset_design), drop = FALSE], .unique.na)
            all <- all[!names(all) %in% names(contr)[1]]

            a <- c(setNames(contr[[1]][1], names(contr[1])), all)
            b <- c(setNames(contr[[1]][2], names(contr[1])), all)

            tmpres <- contrast::contrast(fit = lmfit, a, b, type = "average")

            c("pval" = as.numeric(tmpres$Pvalue), "lfc" = as.numeric(mean(tmpres$foldChange, na.rm = TRUE)))

          } else {
            c("pval" = NA, "lfc" = NA)
          }

        })

      } else {
        contr.res <- matrix(NA, nrow = 2, ncol = length(tmpsubset), dimnames = list(c("pval", "lfc"), names(tmpsubset)))
      }
    })
  })


  ### Results ----

  subset_result <- lapply(results, function(tmp){
    tmpdf <- lapply(names(tmp), function(tmptmp) data.frame(t(tmp[[tmptmp]]), contrast = colnames(tmp[[tmptmp]]), id = tmptmp))
    Reduce(rbind, tmpdf)
  })
  results <- Reduce(rbind, subset_result)
  # if (logged == FALSE){ results$lfc <- log2(results$lfc) }

  results$padj_all <- p.adjust(results$pval, method = p.adj.method)

  results <- split(results, results$contrast)
  results <- lapply(results, function(tmp){
    rownames(tmp) <- tmp$id
    tmp$padj <- p.adjust(tmp$pval, method = p.adj.method)
    tmp$id <- NULL
    tmp$contrast <- NULL
    tmp
  })

  results
}




.logtrans <- function(data, FUN = NULL){

  if (is.null(FUN)){
    FUN <- function(x){
      log2(x + min(min(x[.naf(x > 0)], na.rm = TRUE)/2, 0.5))
    }
  }

  FUN(data)
}



testLIMMA <- function(data, design, formula, contrasts, logged = FALSE, p.adj.method = "fdr", FUN = NULL, ...){

  stopifnot(requireNamespace("limma", quietly = TRUE))

  ### Input arguments ----

  fterms <- labels(terms(formula))
  block <- grep("|", fterms, fixed = TRUE, value = TRUE)
  block <- trimws(gsub("1|\\|", "", block))

  if (length(block) > 0) formula <- update(formula, as.formula(paste0("~ .- (1|", block, ")")))

  if (logged == FALSE){
    data <- .logtrans(data, FUN = FUN)
  }

  data <- as.data.frame(t(data))
  data[is.na(data)] <- NA
  formula <- update(formula,  ~  0 + .)


  ### Contrasts and formula ----

  subsets <- getSubsets(contrasts)


  ### Model fitting and testing ----

  results <- lapply(subsets, function(tmpsubset){

    subset_data <- getDataSubset(data, design, contrast = tmpsubset)
    subset_design <- design[rownames(subset_data), intersect(c(all.vars(formula), block), colnames(design)),drop = FALSE]
    subset_design <- droplevels(na.omit(subset_design))
    subset_data <- subset_data[rownames(subset_design),]

    mm <- model.matrix(formula , subset_design)
    mmcons <- lapply(tmpsubset, function(contr){ paste(paste0(names(contr)[1], contr[[1]]), collapse = " - ") })
    contrasts_limma <- do.call(limma::makeContrasts, c(mmcons, list(levels = mm)))

    if (length(block) > 0){
      blockfactor <- subset_design[[block]]
      dupcor <- limma::duplicateCorrelation(object = t(subset_data), block = blockfactor, design = mm)
      fit <- limma::lmFit(t(subset_data), mm, block = blockfactor, correlation = dupcor$consensus.correlation)
    } else {
      fit <- limma::lmFit(t(subset_data), mm)
    }


    cfit <- limma::contrasts.fit(fit, contrasts = contrasts_limma)
    efit <- limma::eBayes(cfit, ...) ### OPTIONS

    coefnames <- colnames(efit$coefficients)
    sub_res <- lapply(setNames(coefnames, coefnames), function(tmpcoef) limma::topTable(efit, coef = tmpcoef, number = Inf, adjust.method = p.adj.method) )
    sub_res <- lapply(sub_res, function(tmp) data.frame(row.names = rownames(tmp), "pval" = tmp$P.Value, "lfc" = tmp$logFC, "padj" = tmp$adj.P.Val) )

    sub_res
  })


  ### Results ----

  results <- unlist(setNames(results, NULL), recursive = FALSE)
  results

}




testLMM <- function(data, design, formula, contrasts, random = NULL, logged = FALSE, p.adj.method = "fdr", min_pval = .Machine$double.eps, FUN = NULL, ...){

  # variables must be factor

  stopifnot(requireNamespace("nlme", quietly = TRUE))
  stopifnot(requireNamespace("multcomp", quietly = TRUE))

  ### Input ----

  data <- as.data.frame(t(data))
  data[is.na(data)] <- NA

  if (logged == FALSE){
    data <- .logtrans(data, FUN = FUN)
  }

  form.fixed <- update(formula, conc ~ 0 + .)


  subsets <- getSubsets(contrasts)

  results <- lapply(subsets, function(tmpcontrasts){

    subset_data <- getDataSubset(data, design, contrast = tmpcontrasts)

    levels <- lapply(tmpcontrasts, function(tmp) setNames(paste0(tmp[[1]][1], " - ", tmp[[1]][2], " = 0"), names(tmp)[1]))

    class <- labels(terms(formula))
    class_rand <- setdiff(trimws(unlist(strsplit(labels(terms(random)), split = "\\|"))), "1")
    class <- gsub(".*\\(|\\).*", "", class)
    classes <- apply(sapply(class, function(tmp) as.character(design[rownames(subset_data),][[tmp]]) ), 1, paste, collapse = "")

    contr.names <- apply(sapply(tmpcontrasts, "[[", 1), 2, paste, collapse = " - ")

    nadf <- data.frame(row.names = names(tmpcontrasts))
    nadf$lfc <- NA_real_
    # nadf$diff <- NA_real_
    nadf$pval <- NA_real_
    nadf$lfc.mean <- NA_real_
    nadf <- as.matrix(nadf)

    res <- lapply(subset_data, function(conc){

      data.lme <- data.frame("conc" = as.numeric(conc), "classes" = classes, design[rownames(subset_data), union(class, class_rand), drop = FALSE])
      data.lme <- droplevels(na.omit(data.lme))

      if (nrow(data.lme) < 2) return(nadf)

      ix <- sapply(class, function(tmp){
        data.lme[,tmp] %in% names(table(data.lme[,tmp]))[table(data.lme[,tmp]) > 1]
      })
      ix <- apply(ix, 1, function(tmp) Reduce("&", tmp))
      data.lme <- droplevels(data.lme[ix,])


      tmpcontrasts2 <- lapply(tmpcontrasts, function(x) x[sapply(x, length) == 2])
      tmplevels <- levels[sapply(tmpcontrasts2, function(tmp) sapply(names(tmp), function(col) all(tmp[[col]] %in% data.lme[[col]]) ) )]
      if (length(tmplevels) == 0) return(nadf)
      lmecontrasts <- lapply(tmplevels, function(tmp) setNames(list(tmp), names(tmp)) )
      lmecontrasts <- unlist(unlist(lmecontrasts, recursive = FALSE, use.names = FALSE))
      if (length(lmecontrasts) == 1){
        lmecontrasts <- as.list(lmecontrasts)
      } else {
        lmecontrasts <- as.list(unstack(stack(lmecontrasts)))
      }


      if (all(data.lme$conc == 0) | any(table(data.lme[,class]) < 2 )) return(nadf)

      # fit

      if (length(lmecontrasts) == 0) return(nadf)

      lmetest <- FALSE

      tryCatch({
        suppressMessages({
        lmefit <- nlme::lme(fixed = form.fixed, random = random, data = data.lme,
                            control = nlme::lmeControl(maxIter = 1000,
                                                       msMaxEval = 1000,
                                                       msMaxIter = 1000,
                                                       niterEM = 2000,
                                                       returnObject = TRUE))
        lmetest <- multcomp::glht(lmefit, linfct = do.call(multcomp::mcp, lmecontrasts))
        })
      },
        error = function(x){
          message("Warning: Error in model!")
          return(lmetest)})

      if (is.logical(lmetest)) return(nadf)


      # tests
      conc.res <- summary(lmetest, test = multcomp::adjusted("none"))
      dfres <- data.frame("lfc" = conc.res$test$coefficients,
                          "pval" = conc.res$test$pvalues)

      cmat <- lapply(lmecontrasts[[1]], function(x) trimws(strsplit(gsub("=.*", "", x), split = " - ")[[1]]) )


      cmat <- sapply(seq_along(cmat), function(i) paste0(names(cmat)[i], cmat[[i]]) )

      means <- dplyr::group_by(.data = data.lme, classes) %>% dplyr::summarise(mean = mean(conc)) %>% dplyr::pull(mean, name = classes)
      cmat <- sapply(lmecontrasts[[1]], function(x) trimws(strsplit(gsub("=.*", "", x), split = " - ")[[1]]) )
      lfcs.means <- means[cmat[1,]] - means[cmat[2,]]
      dfres$lfc.mean <- as.numeric(lfcs.means)

      rownames(dfres) <- names(contr.names)[match(rownames(dfres), contr.names)]
      dfres <- dfres[names(contr.names),]
      rownames(dfres) <- names(contr.names)
      as.matrix(dfres)


    })



    res <- lapply(setNames(names(tmpcontrasts), names(tmpcontrasts)), function(tmp){
      t(sapply(res, function(x){ x[tmp,] }))
    })

    res


  })


  results <- unlist(setNames(results, NULL), recursive = FALSE)

  results <- lapply(results, function(tmp){
    tmp <- as.data.frame(tmp)
    tmp$pval[.naf(tmp$pval < min_pval)] <- min_pval
    tmp$padj <- p.adjust(tmp$pval, method = p.adj.method)
    tmp
  })

  results
  }










testDREAM <- function(data, design, formula, contrasts, logged = FALSE, p.adj.method = "fdr", FUN = FUN, ...){

  stopifnot(requireNamespace("variancePartition", quietly = TRUE))

  ### Input arguments ----

  data <- as.data.frame(t(data))
  data[is.na(data)] <- NA
  # design <- design[,colnames(design) %in% labels(terms(formula)), drop = FALSE]
  formula <- update(formula, ~  0 + .)


  if (logged == FALSE){
    data <- .logtrans(data, FUN = FUN)
  }

  ### Contrasts and formula ----

  subsets <- getSubsets(contrasts)


  ### Model fitting and testing ----

  results <- lapply(subsets, function(tmpsubset){

    subset_data <- getDataSubset(data, design, contrast = tmpsubset)
    subset_design <- design[rownames(subset_data), all.vars(formula),drop = FALSE]
    subset_design <- na.omit(subset_design)
    subset_data <- subset_data[rownames(subset_design),]

    subset_data <- subset_data[,colSums(is.na(subset_data)) == 0, drop = FALSE]

    vobjDream = variancePartition::voomWithDreamWeights(counts = t(subset_data), formula = formula, data = subset_design, normalize.method = "none")

    contrasts.dream <- sapply(tmpsubset, function(tmpcontr){
      variancePartition::getContrast(vobjDream, formula, subset_design, paste0(names(tmpcontr)[1], tmpcontr[[1]]) )
    })


    fitmm = variancePartition::dream(vobjDream, formula, subset_design, L = contrasts.dream, ...)

    subset_results <- lapply(setNames(colnames(contrasts.dream), colnames(contrasts.dream)), function(tmpname){
      tt <- variancePartition::topTable(fitmm, coef = tmpname, number = Inf, adjust.method = p.adj.method)
      colnames(tt) <- c("lfc", "ave.expr", "t", "pval", "padj", "z.std")
      tt
    })


  })


  ### Results ----

  results <- Reduce('c', results)
  results <- .addListNames(results, name = "contrast", format = "long")
  results$padj_all <- p.adjust(results$pval, method = p.adj.method)
  results <- split(results, results$contrast)
  results <- lapply(results, function(tmp) {
    rownames(tmp) <- tmp$row
    tmp$contrast <- NULL
    tmp$row <- NULL
    tmp} )
  results

}















#
# testBETAREG <- function(data, design, formula, contrasts, p.adj.method = "fdr", logged = FALSE, min_pval = .Machine$double.eps, ...){
#
#   stopifnot(requireNamespace("betareg", quietly = TRUE))
#   stopifnot(requireNamespace("multcomp", quietly = TRUE))
#
#   ### Input arguments ----
#
#   data <- as.data.frame(t(data))
#   data[is.na(data)] <- NA
#   stopifnot(min(data, na.rm = TRUE) >= 0 & max(data, na.rm = TRUE) <= 1)
#
#   formula <- update(formula, conc ~  0 + .)
#
#
#   ### Contrasts and formula ----
#
#   subsets <- getSubsets(contrasts)
#
#   ### Model fitting and testing ----
#
#   results <- lapply(subsets, function(tmpsubset){
#     subset_data <- getDataSubset(data, design, contrast = tmpsubset)
#     subset_design <- design[rownames(subset_data), intersect(all.vars(formula), colnames(design)),drop = FALSE]
#
#     nares <- setNames(rep(NA, length(tmpsubset)), names(tmpsubset))
#     nares <- data.frame(pval = nares, diff = nares)
#
#     res <- lapply(subset_data, function(conc){
#
#       if (length(.unique.na(conc)) <= 1) return(nares)
#       if (any(.naf(conc == 0)) | any(.naf(conc == 1))){
#         conc <- (conc * (length(conc)-1) + 0.5) / length(conc)
#       }
#
#       subset_design$conc <- conc
#       tmp <- subset(subset_design, !is.na(conc))
#
#       # car::linearHypothesis(fit, paste0(groups[1], " = ", groups[2]))
#
#       contr <- lapply(tmpsubset, function(tmpcontr){
#         ix <- sapply(tmpcontr, length) == 2
#         check_df <- tmp[tmp[[names(tmpcontr)[ix]]] %in% tmpcontr[[which(ix)]],, drop = FALSE]
#         if (length(unique(check_df[[which(ix)]])) == 2 & any(table(check_df[[which(ix)]]) > 1)){
#           paste0(names(tmpcontr)[ix], tmpcontr[[which(ix)]])
#         } else {
#           NULL
#         }
#       })
#
#       if (any(!sapply(contr, is.null))){
#         tryCatch({fit <- betareg::betareg(formula, data = subset_design)},
#                  error = function(x){return(nares)})
#       } else {
#         return(nares)
#       }
#
#       contr <- contr[sapply(contr, function(tmp) all(tmp %in% names(fit$coefficients$mean)))]
#
#       contr.glht <- sapply(contr, function(g){paste0(g[1], " - ", g[2], " = 0")})
#       contr.glht <- contr.glht[!sapply(contr, is.null)]
#
#       res <- summary(multcomp::glht(fit, linfct = contr.glht), test = multcomp::adjusted("none"))
#       pval <- setNames(as.numeric(res$test$pvalues), names(contr.glht))[names(tmpsubset)]
#
#       b <- fit$coefficients$mean[unique(unlist(contr))]
#       est <- exp(b)/(1+exp(b))
#       diff <- as.numeric(sapply(contr, function(tmp) setNames(est[tmp[1]] - est[tmp[2]], NULL) ))
#       #lfc <- log2(fc)
#
#
#       data.frame(pval = setNames(pval, names(tmpsubset)), diff)
#     })
#
#     contr.res <- rownames(res[[1]])
#     res <- lapply(setNames(contr.res, contr.res), function(i){
#       data.frame(t(sapply(res, function(tmp){ setNames(as.numeric(tmp[i,]), colnames(tmp)) })))
#     })
#
#     res
#   })
#
#
#   ### Results ----
#
#   results <- Reduce('c', results)
#   results
#
#
# }
#
#



testBETA <- function(data, design, formula, contrasts, zotrans = NULL, p.adj.method = "fdr", verbose = FALSE, logged = FALSE, min_pval = .Machine$double.eps, ...){

  stopifnot(requireNamespace("glmmTMB", quietly = TRUE))
  stopifnot(requireNamespace("multcomp", quietly = TRUE))

  ### Input arguments ----

  data <- as.data.frame(t(data))
  data[is.na(data)] <- NA
  stopifnot(min(data, na.rm = TRUE) >= 0 & max(data, na.rm = TRUE) <= 1)

  formula <- update(formula, conc ~  0 + .)

  if (is.null(zotrans)) zotrans <- is.null(list(...)[["ziformula"]])

  family <- glmmTMB::beta_family(link = "logit")


  ### Contrasts and formula ----

  subsets <- getSubsets(contrasts)

  ### Model fitting and testing ----

  results <- lapply(subsets, function(tmpsubset){
    subset_data <- getDataSubset(data, design[rownames(data),], contrast = tmpsubset)
    subset_design <- design[rownames(subset_data), intersect(all.vars(formula), colnames(design)),drop = FALSE]

    nares <- setNames(rep(NA, length(tmpsubset)), names(tmpsubset))
    nares <- data.frame(pval = nares, diff = nares, diff.mean = nares, lfc = nares)

    res <- lapply(subset_data, function(conc){


      conc_orig <- conc
      if (length(.unique.na(conc)) <= 1) return(nares)
      if (zotrans & (any(.naf(conc == 0)) | any(.naf(conc == 1)))){
        zotrans <- TRUE
        conc <- (conc * (length(conc)-1) + 0.5) / length(conc)
      }

      subset_design$conc <- conc
      tmp <- subset(subset_design, !is.na(conc))

      contr <- lapply(tmpsubset, function(tmpcontr){
        ix <- sapply(tmpcontr, length) == 2
        check_df <- tmp[tmp[[names(tmpcontr)[ix]]] %in% tmpcontr[[which(ix)]],, drop = FALSE]
        if (length(unique(check_df[[which(ix)]])) == 2 & any(table(check_df[[which(ix)]]) > 1)){
          paste0(names(tmpcontr)[ix], tmpcontr[[which(ix)]])
        } else {
          NULL
        }
      })

      if (any(!sapply(contr, is.null))){
        sres <- NULL
        tryCatch({

          # Model
          fit <- glmmTMB::glmmTMB(formula, subset_design, family = family, ...)
          sres <- summary(fit)

        },
        error = function(x){ NULL })
      } else {
        return(nares)
      }

      if (is.null(sres)) return(nares)

      contr <- contr[sapply(contr, function(tmp) all(tmp %in% rownames(sres$coefficients$cond)))]

      contr.glht <- sapply(contr, function(g){paste0(g[1], " - ", g[2], " = 0")})
      if (verbose) .colorcat(contr.glht)
      contr.glht <- contr.glht[!sapply(contr, is.null)]

      res <- summary(multcomp::glht(fit, linfct = contr.glht), test = multcomp::adjusted("none"))
      pval <- setNames(as.numeric(res$test$pvalues), names(contr.glht))[names(tmpsubset)]

      if (verbose) .colorcat(contr.glht, col = "red")
      if (verbose) .colorcat(pval, col = "blue")
      b <- res$coef[unique(unlist(contr))]

      if (zotrans){
        est <- exp(b)/(1+exp(b))
        est <- (est * sres$nobs - 0.5) / (sres$nobs - 1)
      } else {
        est <- b
      }

      diff <- as.numeric(sapply(contr, function(tmp) setNames(est[tmp[1]] - est[tmp[2]], NULL) ))

      # calculate also mean differences for comparison
      f <- lapply(tmpsubset, function(tmp) tmp[sapply(tmp, length) == 2] )
      subset_design$conc_orig <- conc_orig
      diff.mean <- sapply(f, function(tmp){
        tmpdf <- subset_design[subset_design[[names(tmp)]] %in% tmp[[1]],]
        means <- dplyr::group_by_at(tmpdf, names(tmp)) %>%
          dplyr::summarise(diff = mean(conc_orig, na.rm = TRUE))
        means <- setNames(means[["diff"]], means[[names(tmp)]])
        d <- means[tmp[[1]][1]] - means[tmp[[1]][2]]
        setNames(d, NULL)
      })
      if (verbose) .colorcat(diff.mean, col = "red")

      fc <- as.numeric(sapply(contr, function(tmp) setNames(est[tmp[1]] / est[tmp[2]], NULL) ))
      lfc <- log2(fc)

      data.frame(pval = setNames(pval, names(tmpsubset)), diff, diff.mean, lfc)
    })

    contr.res <- rownames(res[[1]])
    res <- lapply(setNames(contr.res, contr.res), function(i){
      df <- data.frame(t(sapply(res, function(tmp){ setNames(as.numeric(tmp[i,]), colnames(tmp)) })))
      df$pval[.naf(df$pval < min_pval)] <- min_pval
      df$padj <- p.adjust(df$pval, method = p.adj.method)
      df
    })

    res
  })


  ### Results ----

  results <- Reduce('c', results)
  results

}





.addListNames <- function(list, name = "sample", format = "list"){

  cols <- lapply(list, colnames)
  cols <- unlist(cols)[!duplicated(unlist(cols))]

  tmp <- lapply(names(list), function(tmpname) cbind(row = rownames(list[[tmpname]]), list[[tmpname]][,cols], tmpname) )
  res <- Reduce(rbind, tmp)
  colnames(res) <- c("row", cols, name)

  if (format %in% c("df", "long")) return(res)

  res <- split(res, res[[name]])

  res <- lapply(res, function(tmp){
    rownames(tmp) <- tmp$row
    tmp$row <- NULL
    tmp
  })

  res
}







