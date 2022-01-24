
# CLASSES



#' TracerExperiment
#'
#' @slot meta List of metadata (optional)
#' @slot colData Data.frame with samples as rownames
#' @slot metData Data.frame
#' @slot isoData Data.frame
#' @slot isoAssays List
#' @slot qcAssays List
#' @slot metAssays List
#' @slot results List
#'
#' @return
#' @export
#'
#' @examples
TracerExperiment <- setClass(Class = "TracerExperiment",
                             slots = c("meta" = "list", # any metadata
                                       "colData" = "data.frame", # sample annotations
                                       "metData" = "data.frame", # metabolite annotations
                                       "isoData" = "data.frame", # isotopologue annotations
                                       "isoAssays" = "list", # isotopologue assays
                                       "qcAssays" = "list", # QC assays
                                       "metAssays" = "list", # metabolite assays
                                       "results" = "list") )




setValidity("TracerExperiment", method = function(object){

  valid <- TRUE
  msg <- NULL

  res <- sapply(object@isoAssays, function(tmp) all(colnames(tmp) == rownames(object@colData)) )
  if (length(res) > 0){
    if (any(res) != TRUE){
      valid <- FALSE
      msg <- "Error: Mismatching sample names."
    }
  }

  res <- sapply(object@isoAssays, function(tmp) all(rownames(tmp) == rownames(object@isoData)) )
  if (length(res) > 0){
    if (any(res) != TRUE){
      valid <- FALSE
      msg <- "Error: Mismatching isotopologue names."
    }
  }

  if (nrow(object@isoData) > 0){
    if (!all(object@isoData$metabolite %in% rownames(object@metData))){
      valid <- FALSE
      msg <- "Error: Missing metabolite annotation(s)."
    }
  }

  if (valid) TRUE else cat(crayon::red(msg))
})













