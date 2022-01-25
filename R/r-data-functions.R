
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








#' Make example TracerExperiment
#'
#' @param nsamples
#' @param nmets
#'
#' @return
#' @export
#'
#' @examples
exampleTracerExperiment <- function(nsamples = 10, nmets = 10){

  samples <- c(paste0(rep("A", floor(nsamples/2)), 1:floor(nsamples/2)), paste0(rep("B", ceiling(nsamples/2)), 1:ceiling(nsamples/2)))
  mets <- replicate(nmets, paste(sample(letters, sample(2:5)), collapse = ""))

  formulas <- sapply(1:length(mets), function(x){
    formula <- paste0("C", sample(1:9, 1), "H", sample(1:9, 1), "O", sample(0:3, 1), "N", sample(0:3, 1), "P", sample(0:3, 1), "S", sample(0:3, 1))
    formula <- gsub(".0", "", formula)
    formula
  })

  metData <- data.frame(row.names = mets, Molecule = gsub("1", "", formulas), MSion = NA)
  colData <- data.frame(row.names = samples, name = samples, group = sub("\\d", "", samples))

  nC <- as.numeric(sub("C", "", regmatches(formulas, regexpr("C[0-9]+", formulas))))
  nC <- setNames(nC, mets)
  isoData <- stack(lapply(nC, function(n) paste0("m", 0:n) ))[,c(2,1)]
  colnames(isoData) <- c("metabolite", "label")
  rownames(isoData) <- paste0(isoData$metabolite, "_", isoData$label)

  midl <- lapply(setNames(mets, mets), function(m){
    tmp <- rmid(ncol = nsamples, nrow = nrow(subset(isoData, metabolite == m)))
    rownames(tmp) <- paste(m, rownames(tmp), sep = "_")
    colnames(tmp) <- samples
    tmp
  })

  data <- data.frame(isoData, Reduce(f = rbind, midl))

  makeTracerExperiment(data, metData = metData, colData = colData)
}





#' Generate a random matrix of mass isotopomer distributions
#'
#' @param nrow
#' @param ncol
#' @param p
#' @param colnames
#' @param size
#'
#' @return
#' @export
#'
#' @examples
rmid <- function(nrow = NULL, ncol = 6, p = NULL, colnames = NULL, size = 10^3){

  if (is.null(nrow) & is.null(p)){
    nrow <- sample(1:20, 1)
  } else if (is.null(nrow) & !is.null(p)){
    nrow <- length(p)
  }

  if (is.null(colnames)){
    n <- round(ncol/2)
    colnames <- c(paste0("A", 1:n), paste0("B", (1):(ncol-n)))
  }

  if (is.null(p)){
    p <- runif(nrow)
    p <- p/sum(p)
  }

  stopifnot(length(p) == nrow)
  names(p) <- paste0("m", 1:nrow-1)

  mid <- rmultinom(n = ncol, size = size, prob = p)/size
  colnames(mid) <- colnames

  as.matrix(mid)
}







.orderLabels <- function(x, pattern = "m"){
  xlevels <- sort(as.numeric(gsub(paste0(".*", pattern), "", unique(x))))
  xlevels <- paste0(pattern, xlevels)
  factor(x, ordered = TRUE, levels = xlevels)
}





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






#
# checkNA <- function(TE, assay = "imp"){
#
#   data <- TE@isoAssays[[assay]]
#   lod <- TE@isoAssays[["lod"]]
#
#   data <- signif(data, 3)
#   data <-  matrix(as.character(unlist(data)), nrow = nrow(data), dimnames = dimnames(data))
#
#   data[is.na(lod)] <- paste0(data[is.na(lod)], "*")
#   data.frame(data)
#
# }




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





stretch <- function(mat){
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
sumMets <- function(TE, assay = "norm", ...){

  assaydata <- data.frame(TE@isoAssays[[assay]], TE@isoData[,c("metabolite"), drop = FALSE])
  assaydata <- assaydata[rowMeans(is.na(TE@isoAssays[[assay]])) != 1,, drop = FALSE]

  res <- .sumAssay(assaydata, ...)
  res[unique(TE@isoData$metabolite),]

}






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











#' Mass isotopomer distribution (MID)
#'
#' @param TE
#' @param assay
#' @param na.rm.sum
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
MID <- function(TE, assay = "norm", na.rm.sum = TRUE, ...){

  data <- cbind(TE@isoAssays[[assay]], TE@isoData[,c("metabolite"), drop = FALSE])
  sumdata <- sumMets(TE, assay = assay, na.rm = na.rm.sum)
  sumdata <- sumdata[data$metabolite,]

  data$metabolite <- NULL
  fractions <- data / sumdata

  fractions[.naf(data == 0)] <- 0
  fractions[is.na(data)] <- NA
  fractions[is.nan(data.matrix(data))] <- NaN

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
isoEnrichment <- function(TE, assay = "norm", na.rm.iso = TRUE, na.rm = FALSE, ...){

  fractions <- MID(TE, assay = assay)
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






subsetMet <- function(data, met = NULL, isodata = NULL){

  if (any(met %in% rownames(data))){
    tmp <- tmp[met,,drop = FALSE]
  } else {
    tmp <- data[rownames(isodata),,drop = FALSE]
    tmp <- tmp[isodata$metabolite %in% met,,drop = FALSE]
  }
  tmp
}







split_by <- function(data, design, formula){

  design <- design[colnames(data),]

  if (length(labels(terms(formula))) > 0){
    groups <- apply(design[,labels(terms(formula)), drop = FALSE], 1, paste0,  collapse = "_")
    stopifnot(all(names(groups) == colnames(data)))
    data_split <- lapply(setNames(unique(groups), unique(groups)), function(g) data[,groups == g, drop = FALSE])
  } else {
    data_split <- list(data)
  }

  data_split
}


unsplit_by <- function(list){

  if (length(list) == 1) return(list[[1]])

  stopifnot(all(Reduce(intersect, lapply(list, rownames)) %in% rownames(list[[1]])))
  list <- lapply(list, function(tmp) tmp[rownames(list[[1]]),, drop = FALSE])
  Reduce(cbind, list)

}


.colorcat <- function(str = "text", col = rgb(1,1,1), add = "\n"){
  cat( crayon::make_style(col)(paste0(str, add)) )
}




.unique.na <- function(x, ...){
  x <- unique(x,  ...)
  x[!is.na(x)]
}





