
# METHODS

#' @include AllGenerics.R


setMethod(f = "show", signature = "TracerExperiment", definition = function(object){
  .colorcat( paste0(class(object), ", ", format(object.size(object), units = "Mb")), col = "#090296")
  col <- "#CFCFCF"
  cat(.colorcat( paste0("isoAssays(", length(object@isoAssays), "):"), add = "", col = col), names(object@isoAssays), "\n")
  cat(.colorcat( paste0("metAssays(", length(object@metAssays), "):"), add = "", col = col), names(object@metAssays), "\n")
  cat(.colorcat( paste0("qcAssays(", length(object@qcAssays), "):"), add = "", col = col), names(object@qcAssays), "\n")
  cat(.colorcat( paste0("results(", length(object@results), "):"), add = "", col = col), names(object@results), "\n")
  col <- "#757575"
  # .colorcat("-use method 'astidy' to get data as a long dataframe", col = col)
  # .colorcat("-use method 'assay' to get data as a 2D matrix", col = col)
  # .colorcat("-use method 'array' to get data as a 3D array", col = col)
})


# Accessor methods

setMethod(f = "colnames", signature = "TracerExperiment", definition = function(x, do.NULL = TRUE, prefix = "col"){
  return(rownames(x@colData))
})

setMethod(f = "rownames", signature = "TracerExperiment", definition = function(x, do.NULL = TRUE, prefix = "row"){
  return(rownames(x@isoData))
})

setMethod(f = "metData", signature = "TracerExperiment", definition = function(object){
  return(object@metData)
})

#' colData
#'
#' @param TracerExperiment
#' @importMethodsFrom SummarizedExperiment colData
#' @return
#' @export
#'
#' @examples
setMethod(f = "colData", signature = "TracerExperiment", definition = function(x, ...){
  return(x@colData)
})

setMethod(f = "isoData", signature = "TracerExperiment", definition = function(object){
  return(object@isoData)
})


setMethod(f = "results", signature = "DESeqDataSet", definition = function(object, ...){
  return(DESeq2::results(object = object, ...))
})


setMethod(f = "results", signature = "TracerExperiment", definition = function(object, ...){

  ix <- list(...)
  res <- object@results

  subix <- 0
  for (i in ix){
    if (length(i) > 1){
      if (subix > 0){
        res <- napply(res, FUN = function(x, i){x[i]}, n = subix, i = i)
      } else {
        res <- res[i]
      }

      subix <- subix + 1
    } else {
      if (subix > 0){
        res <- napply(res, FUN = function(x, i){x[[i]]}, n = subix, i = i)
      } else {
        res <- res[[i]]
      }

    }
  }

  if (length(res) == 1) res <- res[[1]]
  res
})





setMethod(f = "assay", signature = "TracerExperiment", definition = function(x, i, type = "iso", met = NULL, ...){
  stopifnot(tolower(type) %in% c("iso", "met", "qc"))
  if (type == "iso") a <- x@isoAssays[i]
  if (type == "met") a <- x@metAssays[i]
  if (type == "qc") a <- x@qcAssays[[i]]
  if (length(a) == 1){ a <- a[[1]] }

  if (!is.null(met)) a <- subsetMet(a, met = met, isodata = x@isoData)
  a
})


# Replacement methods

setMethod(f = "metData<-", signature = "TracerExperiment", definition = function(object, value){
  object@metData <- value
  if (validObject(object)){
    return(object)
  }
})

setMethod(f = "colData<-", signature = "TracerExperiment", definition = function(x, ..., value){
  x@colData <- value
  if (validObject(x)){
    return(x)
  }
})

setMethod(f = "isoData<-", signature = "TracerExperiment", definition = function(object, value){
  object@isoData <- value
  if (validObject(object)){
    return(object)
  }
})


setMethod(f = "assay<-", signature = "TracerExperiment", definition = function(x, i, withDimnames, ..., value){

  args <- list(...)
  if ("type" %in% names(args)) type <- args[["type"]] else type <- "iso"
  stopifnot(tolower(type) %in% c("iso", "met", "qc"))

  if (type == "iso") x@isoAssays[[i]] <- value
  if (type == "met") x@metAssays[[i]] <- value
  if (type == "qc") x@qcAssays[[i]] <- value

  if (validObject(x)){
    return(x)
  }
})


setMethod(f = "astidy", signature = "TracerExperiment", definition = function(data, ...){
  .tidyTE(data)
})



.tidyTE <- function(TE, type = "iso"){

  if (type == "iso"){
    iso <- lapply(TE@isoAssays, .stretch)
    iso <- lapply(setNames(names(iso), names(iso)), function(tmp) dplyr::rename(.data = iso[[tmp]], !!tmp := Value) )
    isotidy <- Reduce(dplyr::full_join, iso)

    isotidy$JoinID <- isotidy$Sample
    TE@colData$JoinID <- rownames(TE@colData)

    res <- dplyr::full_join(x = isotidy, y = TE@colData, by = "JoinID", suffix = c(".isoData", ".colData"))
    res$JoinID <- NULL

    res <- dplyr::rename(res, isotopologue = Feature)
    TE@isoData$isotopologue <- rownames(TE@isoData)
    res <- dplyr::full_join(x = res, y = TE@isoData, by = "isotopologue")

    TE@metData$metabolite <- rownames(TE@metData)
    res <- dplyr::full_join(x = res, y = TE@metData, by = "metabolite")

  }



  res
}





# Merge and subset methods

# newC13 <- C13donor1 + C13donor2
#
# e1 <- C13donor1
# e2 <- C13donor2

setMethod(f = "+", signature = "TracerExperiment", definition = function(e1, e2){

  #inputnames <- c(rlang::as_name(rlang::enquo(e1)), rlang::as_name(rlang::enquo(e2)))
  #cat("Merging", paste(inputnames, collapse = " and "))
  inputnames <- c("x", "y")
  input <- list(e1, e2)

  # metData
  metdata <- dplyr::full_join(x = tibble::rownames_to_column(e1@metData),
                              y = tibble::rownames_to_column(e2@metData))
  metdata <- tibble::column_to_rownames(metdata, var = "rowname")

  # isoData
  isodata <- dplyr::full_join(x = tibble::rownames_to_column(e1@isoData),
                              y = tibble::rownames_to_column(e2@isoData),
                              by = c("rowname", "metabolite", "label"))
  isodata <- tibble::column_to_rownames(isodata, var = "rowname")
  colnames(isodata) <- gsub(".x$", paste0(".", inputnames[1]), colnames(isodata))
  colnames(isodata) <- gsub(".y$", paste0(".", inputnames[2]), colnames(isodata))

  # colData
  sampledata <- dplyr::full_join(x = tibble::rownames_to_column(e1@colData),
                                 y = tibble::rownames_to_column(e2@colData))
  sampledata <- tibble::column_to_rownames(sampledata, var = "rowname")


  # qcassays - keep only those that have matching colnames
  all_qcassays <- unique(unlist(lapply(input, function(tmp) names(tmp@qcAssays)) ))
  qcassaylist <- lapply(setNames(all_qcassays, all_qcassays), function(tmp) {
    merged <- dplyr::full_join(x = tibble::rownames_to_column(as.data.frame(e1@qcAssays[[tmp]])),
                               y = tibble::rownames_to_column(as.data.frame(e2@qcAssays[[tmp]])),
                               by = "rowname")
    rownames(merged) <- merged$rowname
    if (!all(rownames(sampledata) %in% colnames(merged))) return(NULL)
    merged[rownames(isodata), rownames(sampledata)]
  })
  qcassaylist <- qcassaylist[!sapply(qcassaylist, is.null)]


  # assays - need to have the same dimensions
  all_assays <- unique(unlist(lapply(input, function(tmp) names(tmp@isoAssays)) ))
  assaylist <- lapply(setNames(all_assays, all_assays), function(tmp) {
    merged <- dplyr::full_join(x = tibble::rownames_to_column(as.data.frame(e1@isoAssays[[tmp]])),
                               y = tibble::rownames_to_column(as.data.frame(e2@isoAssays[[tmp]])),
                               by = "rowname")
    rownames(merged) <- merged$rowname
    merged[rownames(isodata), rownames(sampledata)]
  })


  TE <- TracerExperiment(meta = setNames(list(e1@meta, e2@meta), inputnames),
                         metData = metdata,
                         colData = sampledata,
                         isoData = isodata,
                         isoAssays = assaylist,
                         qcAssays = qcassaylist,
                         metAssays = list(),
                         results = list())


  TE
})



setMethod(f = "array", signature = "ANY", definition = function(data, ...){
  A <- base::array(data, ...)
  return(A)
})



setMethod(f = "array", signature = "TracerExperiment", definition = function(data, ...){

  assaydata <- lapply(data@isoAssays, function(tmp){ cbind(data@isoData[,c("metabolite", "label")], tmp) })

  args <- list(...)
  if (length(args) == 0){ args <- list(1)}
  args <- lapply(args, function(tmp){
    if (all(is.numeric(tmp))){ names(assaydata)[tmp] } else { tmp }
  })

  i <- unlist(args)
  A <- lapply(setNames(i,i), function(j) .asArray(data = assaydata[[j]], formula = ~ metabolite + label))

  if (length(A) == 1){ A <- A[[1]] }
  A
})



t.array <- function(x){
  # transpose 3D array (along x-y direction, not altering z)
  if (length(dim(x)) == 2){
    xt <- t.default(x)
  } else {
    xt <- base::array(dim = dim(x)[c(2,1,3)], dimnames = dimnames(x)[c(2,1,3)])
    for (i in 1:nrow(x)){
      xt[,i,] <- x[i,,]
    }
  }
  xt
}



setMethod(f = "$", signature = "TracerExperiment", definition = function(x, name){
  colData(x)[[name]]
})

#' Subset TracerExperiment
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
subset.TracerExperiment <- function(x, factor, type = "col", ...){

  args <- rlang::enquo(factor)

  subsamples <- colnames(x)
  isodf <- isoData(x)
  metdf <- metData(x)

  if (type == "col"){
    # subset based on colData
    subres <- dplyr::filter(colData(x), !!args)
    subsamples <- rownames(subres)

  } else if (type == "met"){
    # subset based on metData
    metdf <- dplyr::filter(metdf, !!args)
    isodf <- isodf[isodf$metabolite %in% rownames(metdf),, drop = FALSE]

  } else if (type == "iso"){
    # subset based on isoData
    isodf <- dplyr::filter(isodf, !!args)
    metdf <- metdf[rownames(metdf) %in% isodf$metabolite,, drop = FALSE]

  }

  mets <- rownames(metdf)
  isos <- rownames(isodf)

  subiso <- lapply(x@isoAssays, function(tmp) tmp[isos, subsamples, drop=FALSE])
  submet <- lapply(x@metAssays, function(tmp) tmp[mets, subsamples, drop=FALSE])
  subqc <- lapply(x@qcAssays, function(tmp) tmp[isos, subsamples, drop=FALSE])

  subTE <- TracerExperiment(meta = x@meta,
                            metData = x@metData[mets,, drop = FALSE],
                            colData = x@colData[subsamples,,drop=FALSE],
                            isoData = x@isoData[isos,, drop = FALSE],
                            isoAssays = subiso,
                            qcAssays = subqc,
                            metAssays = submet)


  subTE
}






setMethod(f = "subset", signature = "TracerExperiment", definition = function(x, ...){

  subset.TracerExperiment(x, ...)

})


setMethod(f = "normalize", signature = "TracerExperiment", definition = function(object, method, assay = NULL, metassay = NULL, thres_LOQ = 1, min_groupfrac_per_iso = 0.9, split_by = ~ 1, ...){

  if (is.null(assay)){
    message("Using default assay 'corr'")
    assay <- "corr"
  }
  if (is.null(metassay)){
    message("Calculating metabolite sums using 'sumMets()'...")
    metsums <- sumMets(object, assay,
                       new_assay = NULL,
                       thres_LOQ = thres_LOQ,
                       max_nafrac_per_met = 0.9,
                       max_nafrac_per_group = 1,
                       min_rep_per_group = 1,
                       min_groupfrac_per_iso = min_groupfrac_per_iso,
                       split_by = split_by,
                       na_iso.rm = TRUE,
                       na.rm = TRUE)
  }

  object@isoAssays$norm <- normalizeIsoAssay(isoAssay = object@isoAssays[[assay]],
                                             metAssay = metsums,
                                             isodata = object@isoData,
                                             fractions = MID(object, assay = assay, new_assay = NULL),
                                             colData = object@colData,
                                             method = method, ...)
  object

})


setMethod(f = "impute", signature = "TracerExperiment", definition = function(object, ...){
  imp <- imputeTE(TE = object, ...)
  object@isoAssays[["imp"]] <- imp
  object
})





#' Split TE
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
split.TracerExperiment <- function(x, ...){

  by <- list(...)[["by"]]
  if (is.null(by)) return(list(x))

  design <- colData(x)
  cols <- labels(terms(by))
  if (length(cols) == 0) return(list(x))

  design <- design[,cols, drop = FALSE]
  rownames(design) <- NULL
  design <- unique(design)

  f <- data.frame(t(design)) %>% as.list()
  names(f) <- sapply(f, function(tmp) paste0(tmp, collapse = "_"))
  fstr <- sapply(f, function(tmp) paste(paste0(cols, " == ", "'", tmp, "'"), collapse = " & ") )

  splitTE <- lapply(fstr, function(tmp, x){
    eval(parse(text = paste0("subset.TracerExperiment(x, ", tmp, ")")), envir = environment())
  }, x = x)

  splitTE


}


setMethod(f = "split", signature = "TracerExperiment", definition = function(object, ...){
  split.TracerExperiment(object, ...)
})


setMethod(f = "split", signature = "ANY", definition = function(object, ...){
  base::split(x = object, ...)
})









