
# GENERICS




#' @export
setGeneric(name = "array", def = function(data, ...){ standardGeneric("array")} )

#' @export
setGeneric(name = "astidy", def = function(data, ...){ standardGeneric("astidy")} )

#' @export
setGeneric(name = "results", def = function(object, ...){ standardGeneric("results")} )

#' @export
setGeneric(name = "metData", def = function(object, ...){ standardGeneric("metData") })

#' @export
setGeneric(name = "metData<-", def = function(object, value){ standardGeneric("metData<-") })

#' @export
setGeneric(name = "isoData", def = function(object, ...){ standardGeneric("isoData") })

#' @export
setGeneric(name = "isoData<-", def = function(object, ...){ standardGeneric("isoData<-") })


#' @export
setGeneric(name = "impute", def = function(object, ...){ standardGeneric("impute")} )


#' colData
#' @importFrom SummarizedExperiment colData
#' @name colData
#' @export
colData <- colData

#' assay
#' @importFrom SummarizedExperiment assay
#' @name assay
#' @export
assay <- assay

#' @importFrom SummarizedExperiment assay<-
#' @name assay<-
#' @export
`assay<-` <- `assay<-`

#' rownames
#' @importFrom BiocGenerics rownames
#' @name rownames
#' @export
rownames <- rownames

#' colnames
#' @importFrom BiocGenerics colnames
#' @name colnames
#' @export
colnames <- colnames

#' @importFrom BiocGenerics normalize
#' @name normalize
#' @export
normalize <- normalize


