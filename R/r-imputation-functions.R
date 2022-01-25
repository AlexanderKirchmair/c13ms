
# IMPUTATION FUNCTIONS



#' Imputation/missing value treatment of isotopologue abundance data
#'
#' @param TE
#' @param assay
#' @param original
#' @param nan
#' @param na
#' @param split_by
#' @param maxNAfrac_per_group
#' @param exclude
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
imputeTE <- function(TE, assay = "lod", original = "raw", nan = NA, na = "", split_by = ~ 1, maxNAfrac_per_group = 0.2, exclude = NULL, type = "iso", ...){

  data.org <-  .getAssays(TE, assay = assay, type = type)
  data <- data.org
  design <- TE@colData

  ### NaN replacement
  if (!is.null(nan)){
    message("Replacing nan")
    data[is.nan(data.matrix(data))] <- nan

  }

  ### NA imputation
  if (na == "missForest"){
    message("missForest imputation")
    res <- missForest::missForest(data)
    imp <- res$ximp
  } else if (na == "LOD50"){
    message("LOD50 imputation")
    lods <- TE@qcAssays$lod
    lods <- lods[rownames(data), colnames(data)]

    imp <- data
    imp[is.na(imp)] <- lods[is.na(imp)] * 0.5

  } else {
    message("No method selected, using original values")
    raw <- .getAssays(TE, assay = original, type = type)
    imp <- data
    imp[is.na(data.matrix(imp))] <- raw[is.na(data.matrix(imp))]
  }




  ### Post-process

  imp <- imp[rownames(data), colnames(data)]

  # Split into groups
  data_split <- split_by(data, design, formula = split_by)
  imp_split <- split_by(imp, design, formula = split_by)

  imp_split <- lapply(seq_along(imp_split), function(tmp){
    ix <- rowMeans(is.na(data_split[[tmp]])) > maxNAfrac_per_group
    imp_split[[tmp]][ix,] <- data_split[[tmp]][ix,]
    imp_split[[tmp]]
  })

  imp <- unsplit_by(imp_split)
  imp <- imp[,colnames(data.org)]
  imp[rownames(imp) %in% exclude,] <- data.org[rownames(imp) %in% exclude,]

  imp

}






















