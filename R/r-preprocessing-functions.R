
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
estimateLODs <- function(TE, qc = "blanks", exclude = "internal.standard", FUN = NULL, ...){

  ### Adds LOD estimates to isoData

  if (is.null(FUN)) FUN <- function(x, sds = 2, ...){  rowMeans(x, na.rm = TRUE) + sds*apply(x, 1, sd, na.rm = TRUE) }

  data <- TE@qcAssays[[qc]]
  if (is.null(data)) stop("Error: QC samples not found!")

  exclude <- intersect(exclude, rownames(data))
  data[exclude] <- 0 # set LOD to 0 for 'excluded' metabolites

  blankmean <- rowMeans(data, na.rm = TRUE)
  LOD <- FUN(data, ...)

  TE@isoData$blankmean <- blankmean[rownames(TE@isoData)]
  TE@isoData$LOD <- LOD[rownames(TE@isoData)]

  TE
}



#' Preprocessing based on LOD/LOQ estimates
#'
#' @param TE
#' @param assay
#' @param split_by
#' @param use_blanks
#' @param blanks
#' @param exclude
#' @param ptres
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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
clean <- function(TE, assay = NULL, split_by = ~ 1, exclude = "internal.standard", maxNAfrac_per_group = 0.2, type = "iso", ...){

  if (is.null(assay)) assay <- .getAssays(TE, assay = assay, last = TRUE, names_only = TRUE, type = type)
  data <- TE@isoAssays[[assay]]
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

  ### NA: Missing values ----
  data_grouped <- lapply(data_grouped, function(tmp){
    tmp[rowMeans(is.na(tmp)) > maxNAfrac_per_group & !rownames(tmp) %in% exclude,] <- NA
    tmp
  })

  data_grouped <- lapply(data_grouped, function(tmp) tmp[rownames(data),])
  data_clean <- Reduce(cbind, data_grouped)
  data_clean <- data_clean[rownames(data), colnames(data)]
  data_clean[is.nan(data.matrix(data))] <- NaN


  TE@isoAssays[[assay]] <- data_clean
  TE

}





























