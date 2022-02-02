
# PREPROCESSING FUNCTIONS



#' Estimate LOD/LOQ values from QC samples
#'
#' @param TE
#' @param qc
#' @param exclude
#' @param FUN
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
estimateLOQs <- function(TE, qc = "blanks", exclude = "internal.standard", sdLOD = 0, sdLOQ = 2, FUN = NULL, ...){

  if (is.null(FUN)) FUN <- function(x, sds = 2, ...){  rowMeans(x, na.rm = TRUE) + sds*apply(x, 1, sd, na.rm = TRUE) }

  data <- TE@qcAssays[[qc]]
  if (is.null(data)) stop("Error: QC samples not found!")

  exclude <- intersect(exclude, rownames(data))

  blankmean <- rowMeans(data, na.rm = TRUE)
  LOD <- FUN(data, sds = sdLOD, ...)
  LOQ <- FUN(data, sds = sdLOQ, ...)

  # set LOD to 0 for 'excluded' metabolites
  LOD[exclude] <- 0
  LOQ[exclude] <- 0

  TE@isoData$blankmean <- blankmean[rownames(TE@isoData)]
  TE@isoData$LOD <- LOD[rownames(TE@isoData)]
  TE@isoData$LOQ <- LOQ[rownames(TE@isoData)]

  LODdf <- data.frame(matrix(rep(LOD, nrow(TE@colData)), ncol = nrow(TE@colData)))
  dimnames(LODdf) <- list(rownames(TE@isoData), rownames(TE@colData))
  TE@qcAssays$lod <- LODdf

  LOQdf <- data.frame(matrix(rep(LOQ, nrow(TE@colData)), ncol = nrow(TE@colData)))
  dimnames(LOQdf) <- list(rownames(TE@isoData), rownames(TE@colData))
  TE@qcAssays$loq <- LOQdf

  TE
}


# old version
preprocessLODs <- function(TE, assay = "raw", split_by = ~ 1, use_blanks = TRUE, blanks = "blanks", exclude = "internal.standard", ptres = 0.05, ...){


  data <- .getAssays(TE, assay = assay, type = "iso")
  design <- TE@colData
  if (is.null(data)) stop("Error: Data not found!")
  exclude <- intersect(exclude, rownames(data))
  if (length(exclude) == 0) exclude <- NULL

  # Split into groups
  if (length(labels(terms(split_by))) > 0){
    groups <- apply(design[,labels(terms(split_by)),drop = FALSE], 1, paste0,  collapse = "_")
    stopifnot(all(names(groups) == colnames(data)))
    data_grouped <- lapply(setNames(unique(groups), unique(groups)), function(g) data[,groups == g,drop = FALSE])
  } else {
    data_grouped <- list(data)
  }



  ### NaN: Test whether sample data are significantly above blanks ----

  blankdata <- TE@qcAssays[[blanks]]
  above_blank <- NULL

  if (use_blanks == TRUE & !is.null(blankdata)){

    pval_ab <- data.frame(sapply(data_grouped, function(tmp){
      fac <- relevel(factor( c(rep("blank", ncol(blankdata)), rep("sample", ncol(tmp)))), ref = "blank")
      testdata <- data.frame(blankdata[rownames(tmp),], tmp)
      res <- apply(testdata, 1, function(x){
        n <- sapply(unique(fac), function(ff) sum(!is.na(x[ff == fac])))
        if (all(n > 1)) t.test(x ~ fac, alternative = "less", na.action = "na.omit")$p.value else NA
      })
      res[matrixStats::rowAlls(testdata == 0)] <- 0
      if (!is.null(exclude)) res[rownames(testdata) %in% exclude] <- NA
      res
    }))


    # Data not above blanks:
    above_blank <- lapply(.usenames(data_grouped), function(g){
      tmp <- data_grouped[[g]]
      mat <- matrix( rep(nat(pval_ab[[g]] <= ptres), ncol(tmp)), ncol = ncol(tmp))
      dimnames(mat) <- dimnames(tmp)
      mat
    })

  }


  ### LODs: for each isotopologue ----
  LODs <- setNames(TE@isoData$LOD, rownames(TE@isoData))
  above_LOD <- lapply(data_grouped, function(tmp){
    tmp <- tmp > LODs[rownames(tmp)]
    if (!is.null(exclude)) tmp[rownames(tmp) %in% exclude] <- TRUE
    tmp
  })


  ### Combine

  if (is.null(above_blank)){ # if test with blank means is not done

    data_grouped <- lapply(.usenames(data_grouped), function(g){
      tmp <- data_grouped[[g]]
      tmp[!above_LOD[[g]][,colnames(tmp)]] <- NA # set all values below LOD to NA
      tmp
    })

  } else { # if test with blank means is done

    data_grouped <- lapply(.usenames(data_grouped), function(g){
      tmp <- data_grouped[[g]]
      tmp[!above_LOD[[g]][,colnames(tmp)] & !above_blank[[g]][,colnames(tmp)] & !is.na(tmp)] <- NaN # set all values below blank means and below LOD to NaN
      tmp[!above_LOD[[g]][,colnames(tmp)] & above_blank[[g]][,colnames(tmp)] & !is.na(tmp)] <- NA # set all values above blank means but below LOD to NA
      tmp
    })

  }

  data_grouped <- lapply(data_grouped, function(tmp) tmp[rownames(data),])
  data_lod <- Reduce(cbind, data_grouped)
  data_lod <- data_lod[rownames(data), colnames(data)]

  # keep zeros
  data_lod[naf(data == 0)] <- 0

  TE@isoAssays$lod <- data_lod

  # confidence values
  tmp <- TE@isoAssays[[assay]] - TE@isoData$blankmean
  tmp[tmp < 0] <- 0
  TE@qcAssays$conf <- tmp / TE@isoData$LOD
  TE@qcAssays$conf[tmp == 0] <- 0


  # add results in assay format
  n <- ncol(TE@isoAssays[[assay]])
  TE@qcAssays$blank <- data.frame(matrix(rep(TE@isoData$blankmean, n), ncol = n, dimnames = dimnames(TE@isoAssays[[1]])))
  TE@qcAssays$nan <- is.nan(data.matrix(TE@isoAssays$lod))
  TE@qcAssays$na <- is.na(data.matrix(TE@isoAssays$lod)) & !is.nan(data.matrix(TE@isoAssays$lod))

  TE
}






#' Preprocessing based on LOD/LOQ estimates
#'
#' @param TE
#' @param assay
#' @param LODassay
#' @param LOQassay
#' @param exclude
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
preprocessLOQs <- function(TE, assay = "raw", LODassay = "lod", LOQassay = "loq", exclude = "internal.standard", ...){

  data <- .getAssays(TE, assay = assay, type = "iso")
  if (is.null(data)) stop("Error: Data not found!")
  data_preprocessed <- data

  exclude <- intersect(exclude, rownames(data))
  if (length(exclude) == 0) exclude <- NULL

  LODs <- .getAssays(TE, assay = LODassay, type = "qc")
  LOQs <- .getAssays(TE, assay = LOQassay, type = "qc")

  above_LOQ <- NULL
  if (!is.null(LOQassay)){
    stopifnot(all.equal(dimnames(LOQs), dimnames(data)))
    if (!is.null(exclude)) LOQs[rownames(LOQs) %in% exclude,] <- 0
    above_LOQ <- .naf(data > LOQs)
    data_preprocessed[!above_LOQ] <- NA
  }

  above_LOD <- NULL
  if (!is.null(LODassay)){
    stopifnot(all.equal(dimnames(LODs), dimnames(data)))
    if (!is.null(exclude)) LODs[rownames(LODs) %in% exclude,] <- 0
    above_LOD <- .naf(data > LODs)
    data_preprocessed[!above_LOD] <- NaN
  }


  # keep zeros
  data_preprocessed[naf(data == 0)] <- 0
  TE@isoAssays$lod <- data_preprocessed

  # confidence values (0 = LOD, 1 = LOQ)
  conf <- (data - LODs)/(LOQs - LODs)
  conf[LOQs == 0] <- Inf
  TE@qcAssays$conf <- conf

  # add results in assay format
  TE@qcAssays$nan <- is.nan(data.matrix(TE@isoAssays$lod))
  TE@qcAssays$na <- is.na(data.matrix(TE@isoAssays$lod)) & !is.nan(data.matrix(TE@isoAssays$lod))

  TE
}














#' Cleaning of data
#'
#' @param TE
#' @param assay
#' @param split_by
#' @param exclude
#' @param maxNAfrac_per_group
#'
#' @return
#' @export
#'
#' @examples
clean <- function(TE, assay = "norm", split_by = ~ 1, exclude = "internal.standard", remove_imp = TRUE, max_na = 1/3, type = "iso", ...){

  if (is.null(assay)) assay <- .getAssays(TE, assay = assay, last = TRUE, names_only = TRUE, type = type)
  data <- TE@isoAssays[[assay]]
  design <- TE@colData
  if (is.null(data)) stop("Error: Data not found!")
  exclude <- intersect(exclude, rownames(data))
  if (length(exclude) == 0) exclude <- NULL



  # Remove imputed values
  if (remove_imp == TRUE & !is.null(TE@qcAssays$na)){
    .colorcat("Setting imputed values to NA...")
    imp <- as.matrix(TE@qcAssays$na)
    data[imp] <- NA
  }


  # Split into groups
  if (length(labels(terms(split_by))) > 0){
    groups <- apply(design[,labels(terms(split_by)),drop = FALSE], 1, paste0,  collapse = "_")
    stopifnot(all(names(groups) == colnames(data)))
    data_grouped <- lapply(setNames(unique(groups), unique(groups)), function(g) data[,groups == g,drop = FALSE])
  } else {
    data_grouped <- list(data)
  }

  ### NA: Missing values ----
  data_grouped <- lapply(data_grouped, function(tmp){
    tmp[rowMeans(is.na(tmp)) > max_na & !rownames(tmp) %in% exclude,] <- NA
    tmp
  })

  data_grouped <- lapply(data_grouped, function(tmp) tmp[rownames(data),])
  data_clean <- Reduce(cbind, data_grouped)
  data_clean <- data_clean[rownames(data), colnames(data)]
  data_clean[is.nan(data.matrix(data))] <- NaN


  TE@isoAssays[["clean"]] <- data_clean
  TE

}



