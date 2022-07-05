
# DATA HANDLING FUNCTIONS



#' Make new TracerExperiment
#'
#' @param data
#' @param metabolite
#' @param label
#' @param QC
#' @param metData
#' @param colData
#' @param meta
#' @param assay
#'
#' @return
#' @export
#'
#' @examples
makeTracerExperiment <- function(data, metabolite = metabolite, label = label, QC = NULL, metData = NULL, colData = NULL, meta = NULL, assay = "raw"){


  if (missing(data)){

    assayData <- list()
    isoData <- data.frame()
    colData <- data.frame()

  } else {

    lab <- rlang::enquo(label)
    met <- rlang::enquo(metabolite)

    data <- dplyr::rename(.data = data, label = !!lab, metabolite = !!met)
    data$label <- .orderLabels(data$label)
    isoData <- dplyr::select(data, c(metabolite, label))
    assayData <- list(dplyr::select(data, -c(metabolite, label)))
    names(assayData) <- assay

  }

  if (is.null(meta)) meta <- list()
  if (is.null(metData)) metData <- data.frame(row.names = unique(isoData$metabolite))
  if (is.null(colData)) colData <- data.frame(row.names = colnames(assayData[[1]]))
  if (is.null(QC)) QC <- list()

  TracerExperiment(meta = meta,
                   colData = colData,
                   metData = metData,
                   isoData = isoData,
                   isoAssays = assayData,
                   qcAssays = QC,
                   metAssays = list())

}












#' Order isotopologue labels numerically
#'
#' @param x
#' @param pattern
#'
#' @return
#'
.orderLabels <- function(x, pattern = "m"){
  xlevels <- sort(as.numeric(gsub(paste0(".*", pattern), "", unique(x))))
  xlevels <- paste0(pattern, xlevels)
  factor(x, ordered = TRUE, levels = xlevels)
}





#' Function to get default assay
#'
#' @param TE
#' @param type
#' @param last
#' @param assay
#' @param not
#' @param names_only
#' @param ...
#'
#' @return
#'
.getAssays <- function(TE, type = "iso", last = FALSE, assay = NULL, not = NULL, names_only = FALSE, ...){

  res <- list()
  if (tolower(type) == "iso") res <- TE@isoAssays
  if (tolower(type) == "met") res <- TE@metAssays
  if (tolower(type) == "qc") res <- TE@qcAssays


  if (is.character(assay)){
    if (!assay %in% names(res)){
      msg <- paste0("Error: ", paste0(paste0("'", assay, "'"), collapse = "/"), " not in available assays (",
                    paste0(names(res), collapse = "/"), ") of type ", type, "!")
      stop(msg)
    }
  }

  if (!is.null(not)){
    res <- res[!names(res) %in% not]
  }

  if (!is.null(assay)){
    res <- res[assay]
  } else if (last == TRUE){
    res <- res[length(res)]
  }

  if (names_only == TRUE){
    return(names(res))
  }

  if (length(res) == 1) res <- res[[1]]
  res
}



.usenames <- function(x){
  setNames(names(x), names(x))
}




.asArray <- function(data, formula = value ~ dim1 + dim2 + dim3){

  dims <- labels(terms(formula))
  var <- setdiff(all.vars(formula), dims)
  if (length(var) == 0) var <- setdiff(colnames(data), dims)

  df <- dplyr::select(data, c(!!dims, !!var))

  if (length(dims) > 2){
    df <- tidyr::pivot_wider(df, names_from = dims[2], values_from = var)
    i <- 3
  } else {
    i <- 2
  }

  dx <- unique(df[[dims[i]]])
  res <- lapply(setNames(dx, as.character(dx)), function(d) tidyr::pivot_wider(subset(df, df[[dims[i]]] == d), names_from = dims[i]) )
  res <- lapply(res, function(x) tibble::column_to_rownames(x, dims[1]))
  names(res) <- dx

  dall <- lapply(setNames(dims, dims), function(d) unique(df[[d]]) )
  dall$val <- var
  A <- array(NA, dim = sapply(dall, length), dimnames = dall)

  for (i in names(res)){
    A[,i,] <- as.matrix(res[[i]][dall[[1]], dall[[3]]])
  }

  A
}



#' Metabolite names
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
metnames <- function(x, ...){
  rownames(x@metData)
}









#' Simple pivot_longer wrapper for numeric matrices
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
.stretch <- function(mat){
  long <- as.data.frame(tidyr::pivot_longer(dplyr::mutate(mat, id = rownames(mat)), -id))
  rownames(long) <- paste0(long$name, "_", long$id)
  long <- dplyr::rename(.data = long, Sample = name, Feature = id, Value = value)
  long
}





#' Summarize isotopologues to metabolites
#'
#' @param TE
#' @param assay
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sumMets <- function(TE, assay = "norm", max_na = 0.1, max_imp = 0.25, ...){

  assaydata <- data.frame(TE@isoAssays[[assay]])
  mets <- TE@isoData[,c("metabolite"), drop = FALSE]

  sumdata <- .sumAssay(cbind(mets, assaydata), na.rm = TRUE)

  # remove metabolites with too many missing isotopologues
  nafraction <- .sumAssay(cbind(mets, is.na(assaydata)), FUN = mean)
  sumdata[nafraction > max_na] <- NA

  # remove metabolites with too many imputed isotopologues
  if (!is.null(TE@qcAssays$na) & !is.null(max_imp)){
    imp <- TE@qcAssays$na
    impfraction <- .sumAssay(cbind(mets, imp), FUN = mean)
    sumdata[impfraction > max_imp] <- NA
  }

  sumdata[unique(TE@isoData$metabolite),]
}






#' Mass isotopomer distribution (MID)
#'
#' @param TE
#' @param assay
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
MID <- function(TE, assay = "norm", max_na = 0, max_imp = 0.25, remove_imp = TRUE, ...){

  data <- cbind(TE@isoAssays[[assay]], TE@isoData[,c("metabolite"), drop = FALSE])
  sumdata <- sumMets(TE, assay = assay, max_na = max_na, max_imp = max_imp, ...)
  sumdata <- sumdata[data$metabolite,]

  data$metabolite <- NULL
  fractions <- data / sumdata

  fractions[.naf(data == 0)] <- 0
  fractions[is.na(data)] <- NA
  fractions[is.nan(data.matrix(data))] <- NaN

  if (remove_imp == TRUE & !is.null(TE@qcAssays$na)){
    imp <- as.matrix(TE@qcAssays$na)
    fractions[imp] <- NA
  }

  stopifnot(all.equal(dimnames(data), dimnames(fractions)))
  fractions
}








#' Fractional enrichment
#'
#' @param TE
#' @param assay
#' @param na.rm.iso
#' @param na.rm
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
isoEnrichment <- function(TE, assay = "norm", max_na = 0, max_imp = 0.25, na.rm.iso = TRUE, na.rm = FALSE, ...){

  fractions <- MID(TE, assay = assay, max_na = max_na, max_imp = max_imp, ...)
  fractions$id <- rownames(fractions)

  df <- data.frame(id = rownames(TE@isoData), metabolite = TE@isoData$metabolite)
  df <- dplyr::full_join(df, fractions, by = "id")
  rownames(df) <- df$id

  dfl <- split(df, df$metabolite)

  frac_enrich <- t(sapply(dfl, function(tmp){
    w <- seq(0, 1, length = nrow(tmp))
    if (na.rm.iso) tmp[rowMeans(is.na(tmp)) == 1,-(1:2)] <- 0
    colSums(tmp[,-(1:2)] * w, na.rm = na.rm)
  }))

  if (!is.null(TE@metData)) frac_enrich <- frac_enrich[rownames(TE@metData),]


  frac_enrich
}






#' Subset isotobologue data by metabolites
#'
#' @param data
#' @param met
#' @param isodata
#'
#' @return
#' @export
#'
#' @examples
subsetMet <- function(data, met = NULL, isodata = NULL){

  if (any(met %in% rownames(data))){
    tmp <- tmp[met,,drop = FALSE]
  } else {
    tmp <- data[rownames(isodata),,drop = FALSE]
    tmp <- tmp[isodata$metabolite %in% met,,drop = FALSE]
  }
  tmp
}





#
#
# split_by <- function(data, design, formula){
#
#   design <- design[colnames(data),]
#
#   if (length(labels(terms(formula))) > 0){
#     groups <- apply(design[,labels(terms(formula)), drop = FALSE], 1, paste0,  collapse = "_")
#     stopifnot(all(names(groups) == colnames(data)))
#     data_split <- lapply(setNames(unique(groups), unique(groups)), function(g) data[,groups == g, drop = FALSE])
#   } else {
#     data_split <- list(data)
#   }
#
#   data_split
# }
#
#
# unsplit_by <- function(list){
#
#   if (length(list) == 1) return(list[[1]])
#
#   stopifnot(all(Reduce(intersect, lapply(list, rownames)) %in% rownames(list[[1]])))
#   list <- lapply(list, function(tmp) tmp[rownames(list[[1]]),, drop = FALSE])
#   Reduce(cbind, list)
#
# }









.sumAssay <- function(data, var = metabolite, FUN = sum, ...){

  var <- rlang::enquo(var)

  sumdata <- data %>%
    dplyr::group_by(!!var) %>%
    dplyr::summarise(dplyr::across(.fns = FUN, ...)) %>%
    tibble::column_to_rownames(var = rlang::as_name(var)) %>%
    as.data.frame()

  sumdata
}






.naf <- function (data, ...){
  data[is.na(data)] <- FALSE
  data
}





.cutMID <- function(data){
  data[data < 0] <- 0

  if (is.null(dim(data))){
    if (any(!is.na(data))){
      ind <- 1:length(data)
      last <- max(ind[!is.na(data)])
      data[1:last]
    } else {
      NA
    }

  } else {
    if (any(!is.na(data))){
      ind <- 1:ncol(data)
      last <- max(ind[ as.logical(colSums(!is.na(data)))])
      data[,1:last, drop = FALSE]
    } else {
      NA
    }
  }

}



.colorcat <- function(str = "text", col = rgb(1,1,1), add = "\n"){
  cat( crayon::make_style(col)(paste0(str, add)) )
}



.unique.na <- function(x, ...){
  x <- unique(x,  ...)
  x[!is.na(x)]
}


