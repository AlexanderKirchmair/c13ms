
# GENERICS

#' @include AllClasses.R


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




