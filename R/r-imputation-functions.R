
# IMPUTATION FUNCTIONS



#' Imputation/missing value treatment of isotopologue abundance data
#'
#' @param TE
#' @param assay cleaned data with missing values
#' @param original raw values without missing values
#' @param nan replacement for NaN values (default = NA)
#' @param na replacement for NA values (number, method or function)
#' @param exclude
#' @param qc_LOD
#' @param qc_LOQ
#' @param type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' exampleTracerExperiment(add_qc = TRUE) |> estimateLOQs() |> preprocessLOQs() |> impute(nan = 0, na = "linsubLODs")
imputeTE <- function(TE, assay = "lod", original = "raw", nan = NA, na = c("subLODs", "linsubLODs", "LOD50", "missForest"), qc_LOD = "lod", qc_LOQ = "loq", exclude = "internal.standard", type = "iso", ...){


  data.org <-  .getAssays(TE, assay = assay, type = type)
  data <- data.org
  design <- TE@colData

  LOD <- TE@qcAssays[[qc_LOD]]
  LOQ <- TE@qcAssays[[qc_LOQ]]

  LOD <- LOD[rownames(data), colnames(data)]
  LOQ <- LOQ[rownames(data), colnames(data)]

  raw <- .getAssays(TE, assay = original, type = type)

  if (length(na) > 1) na <- na[1]


  # NaN values
  if (!is.null(nan)){
    message("Replacing nan")
    data[is.nan(data.matrix(data))] <- nan
  }

  # NA values
  if (is.numeric(na)){
    data[is.na(data.matrix(data))] <- na

  } else if (is.character(na)){
    # na should be the name of the function
    # any NaNs still present are treated as NAs
    imp <- do.call(what = paste0(".imp_", tolower(na)), args = list(data = data, LOD = LOD, LOQ = LOQ, raw = raw, ...))
    ix <- is.na(data.matrix(data))
    data[ix] <- imp[ix]

  } else if (is.function(na)){
    imp <- na(data, ...)
    ix <- is.na(data.matrix(data))
    data[ix] <- imp[ix]


  } else {
    data <- data
  }

  data[rownames(data) %in% exclude,] <- data.org[rownames(data.org) %in% exclude,]
  data.frame(data)

}





### FUNCTIONS ----


# Replace data by LOD-subtracted original values
.imp_sublods <- function(data, LOD, raw, ...){
  imp <- raw - LOD
  imp[.naf(imp < 0)] <- 0
  imp
}


# Weighted subtraction of LODs from raw
.imp_linsublods <- function(LOD, LOQ, raw, ...){
  imp <- raw - LOD / LOQ
  imp[.naf(imp < 0)] <- 0
  imp
}


# Replace data by LOD50 values
.imp_lod50 <- function(LOD, ...){
  LOD * 0.5
}


# MissForest imputation
.imp_missforest <- function(data, ...){
  res <- missForest::missForest(data)
  res$ximp
}






