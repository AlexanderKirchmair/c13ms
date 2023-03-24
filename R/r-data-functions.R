
# DATA HANDLING FUNCTIONS



#' Make new TracerExperiment
#'
#' @param data raw isotopologue abundances
#' @param metabolite metabolite column
#' @param label label column
#' @param QC QC assays
#' @param metData metabolite annotation
#' @param colData sample annotation
#' @param meta metadata
#' @param assay initial assay name
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
#' @param TE TracerExperiment
#' @param assay norm
#' @param thres_LOQ LOQ threshold (values below are set to NA)
#' @param max_nafrac_per_met maximum allowed fraction of isotopologue NA values per metabolite
#' @param max_nafrac_per_group maximum allowed fraction of NA values per group
#' @param min_rep_per_group minimum required number of non-NA values per group
#' @param min_groupfrac_per_iso minimum required number of non-NA groups per isotopologue
#' @param split_by factor(s) from colData to define groups
#' @param exclude
#' @param na_iso.rm
#' @param na.rm
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' exampleTracerExperiment() |> sumMets(assay = "raw", qc_LOQ = NULL)
#' exampleTracerExperiment(add_qc = TRUE) |> estimateLOQs() |> preprocessLOQs() |> sumMets(assay = "lod", sum_qc = TRUE)
sumMets <- function(TE, assay = "norm", new_assay = "", thres_LOQ = 1.5, qc_LOQ = "loq", max_nafrac_per_met = 0.25, max_nafrac_per_group = 0.1, min_rep_per_group = 2, min_groupfrac_per_iso = 0.8, split_by = ~ 1, exclude = "internal.standard", sum_qc = NULL, na_iso.rm = TRUE, na.rm = TRUE, ...){

  if (!is.null(new_assay)){ if (new_assay == ""){ new_assay <- assay} }
  design <- TE@colData
  mets <- TE@isoData[,c("metabolite"), drop = FALSE]

  if (is.null(sum_qc)){sum_qc <- !is.null(qc_LOQ) }
  if (!qc_LOQ %in% names(TE@qcAssays)){
    qc_LOQ <- NULL
  .colorcat("Warning: No QC assay found.")
  }


  assaydata <- clean(TE, assay = assay,
                     new_assay = NULL,
                     qc_LOQ = qc_LOQ,
                     thres_LOQ = thres_LOQ,
                     max_nafrac_per_group = max_nafrac_per_group,
                     min_rep_per_group = min_rep_per_group,
                     split_by = split_by,
                     exclude = exclude)

  # exclude isotopologues with too many NA-only groups
  if (length(labels(terms(split_by))) > 0){
    groups <- apply(design[,labels(terms(split_by)),drop = FALSE], 1, paste0,  collapse = "_")
    stopifnot(all(names(groups) == colnames(assaydata)))
    data_grouped <- lapply(setNames(unique(groups), unique(groups)), function(g) assaydata[,groups == g,drop = FALSE])
  } else {
    data_grouped <- list(assaydata)
  }
  data_grouped <- lapply(data_grouped, function(tmp) tmp[rownames(assaydata),])
  assaydata_clean <- Reduce(cbind, data_grouped)
  assaydata_clean <- assaydata_clean[rownames(assaydata), colnames(assaydata)]
  if (!is.null(min_groupfrac_per_iso)){
    group_na <- sapply(data_grouped, function(tmp){ rowSums(!is.na(tmp)) })
    assaydata_clean[!rowMeans(group_na > min_rep_per_group) >= min_groupfrac_per_iso,] <- NA
  }

  # sum data
  tmp <- data.frame(mets, assaydata_clean)
  if (na_iso.rm == TRUE) tmp <- tmp[rowSums(!is.na(assaydata_clean)) > 0,]
  sumdata <- .sumAssay(tmp, na.rm = na.rm)

  # remove metabolites with too many missing isotopologues
  nafraction <- .sumAssay(data.frame(tmp[,1, drop = FALSE], is.na(tmp[,-1])), FUN = mean)
  sumdata[nafraction > max_nafrac_per_met] <- NA

  # set metabolite order
  mets_uni <- unique(mets[,1])
  if (all(mets_uni %in% rownames(TE@metData))) mets_uni <- rownames(TE@metData)
  sumdata <- sumdata[mets_uni,]
  rownames(sumdata) <- mets_uni

  if (is.null(new_assay)) return(sumdata)

  # combine LOQ values
  if (sum_qc == TRUE & !is.null(qc_LOQ)){
    LOQ <- .getAssays(TE, assay = qc_LOQ, type = "qc")
    LOQ <- LOQ[rownames(assaydata_clean), colnames(assaydata_clean)]
    LOQ[is.na(assaydata_clean)] <- NA
    tmp_loq <- data.frame(mets, LOQ)
    tmp_loq <- tmp_loq[rownames(tmp),]
    sumloq <- .sumAssay(tmp_loq, na.rm = TRUE, FUN = mean)
    sumloq <- sumloq[mets_uni,]
    rownames(sumloq) <- mets_uni

    TE@qcAssays[[paste0("met_loq_", new_assay)]] <- sumloq
  }


  TE@metAssays[[new_assay]] <- sumdata
  TE
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
MID <- function(TE, assay = "norm", metAssay = NULL, new_assay = "mid", ...){

  data <- cbind(TE@isoAssays[[assay]], TE@isoData[,c("metabolite"), drop = FALSE])

  if (is.null(metAssay)){
    sumdata <- sumMets(TE, assay = assay, new_assay = NULL, ...)
  } else {
    sumdata <- .getAssays(TE, type = "met", assay = metAssay)
  }

  sumdata <- sumdata[data$metabolite,]

  data$metabolite <- NULL
  fractions <- data / sumdata

  fractions[.naf(data == 0)] <- 0
  fractions[is.na(data)] <- NA
  fractions[is.nan(data.matrix(data))] <- NaN

  stopifnot(all.equal(dimnames(data), dimnames(fractions)))

  if (is.null(new_assay)) return(fractions)
  TE@isoAssays[[new_assay]] <- fractions
  TE
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
isoEnrichment <- function(TE, assay = "mid", new_assay = "frac", na.rm = TRUE, na.rm.iso = TRUE, ...){


  fractions <- .getAssays(TE, type = "iso", assay = assay)
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

  if (is.null(new_assay)) return(frac_enrich)
  TE@metAssays[[new_assay]] <- frac_enrich
  TE
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


