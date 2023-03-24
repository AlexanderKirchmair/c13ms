
# PREPROCESSING FUNCTIONS



#' Estimate LOD and sdLOD values from QC samples
#'
#' @param TE TracerExperiment
#' @param qc name of QC assay
#' @param exclude metabolites without LOD
#' @param FUN function to calculate LOD from QC assay
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' exampleTracerExperiment(add_qc = TRUE) |> estimateLOQs()
estimateLOQs <- function(TE, assay = "raw", qc = "blanks", qc_LOD = "lod", qc_LOQ = "loq", exclude = "internal.standard", sds = 1, FUN = NULL, ...){

  if (is.null(FUN)) FUN <- function(x, sds = 0, ...){  rowMeans(x, na.rm = TRUE) + sds*apply(x, 1, sd, na.rm = TRUE) }

  data <- .getAssays(TE, assay = assay, type = "iso")
  qcdata <- TE@qcAssays[[qc]]
  if (is.null(qcdata)) stop("Error: QC samples not found!")

  exclude <- intersect(exclude, rownames(qcdata))

  LOD <- FUN(qcdata, sds = 0, ...)
  sdLOD <- FUN(qcdata, sds = sds, ...)

  # set LOD to 0 for 'excluded' metabolites
  LOD[exclude] <- 0
  sdLOD[exclude] <- 0

  TE@isoData[[qc_LOD]] <- LOD[rownames(TE@isoData)]
  # TE@isoData[[qc_sdLOD]] <- sdLOD[rownames(TE@isoData)]

  LODdf <- data.frame(matrix(rep(LOD, nrow(TE@colData)), ncol = nrow(TE@colData)))
  dimnames(LODdf) <- list(rownames(TE@isoData), rownames(TE@colData))
  TE@qcAssays[[qc_LOD]] <- LODdf

  # calculate LOQ values (i.e., number of sds above LOD)
  LOQ <- (data - data.matrix(LODdf))/sdLOD
  LOQ[.naf(sdLOD == 0),] <- Inf
  TE@qcAssays[[qc_LOQ]] <- LOQ

  # sdLODdf <- data.frame(matrix(rep(sdLOD, nrow(TE@colData)), ncol = nrow(TE@colData)))
  # dimnames(sdLODdf) <- list(rownames(TE@isoData), rownames(TE@colData))
  # TE@qcAssays[[qc_sdLOD]] <- sdLODdf

  TE
}




#' Preprocessing based on LOD/LOQ estimates
#'
#' @param TE TracerExperiment
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
#' exampleTracerExperiment(add_qc = TRUE) |> estimateLOQs() |> preprocessLOQs()
preprocessLOQs <- function(TE, assay = "raw", new_assay = "lod", thres_LOQ = 1.5, qc_LOD = "lod", qc_LOQ = "loq", below_LOQ = NA, below_LOD = NaN, exclude = "internal.standard", ...){

  data <- .getAssays(TE, assay = assay, type = "iso")
  if (is.null(data)) stop("Error: Data not found!")
  data_preprocessed <- data

  exclude <- intersect(exclude, rownames(data))
  if (length(exclude) == 0) exclude <- NULL

  LOD <- .getAssays(TE, assay = qc_LOD, type = "qc")
  LOQ <- .getAssays(TE, assay = qc_LOQ, type = "qc")

  above_LOQ <- NULL
  if (!is.null(qc_LOQ)){
    stopifnot(all.equal(dimnames(LOQ), dimnames(data)))
    if (!is.null(exclude)) LOQ[rownames(LOQ) %in% exclude,] <- 0
    above_LOQ <- .naf(data > LOQ)
    data_preprocessed[!above_LOQ] <- below_LOQ
  }

  above_LOD <- NULL
  if (!is.null(qc_LOD)){
    stopifnot(all.equal(dimnames(LOD), dimnames(data)))
    if (!is.null(exclude)) LOD[rownames(LOD) %in% exclude,] <- 0
    above_LOD <- .naf(data > LOD)
    data_preprocessed[!above_LOD] <- below_LOD
  }

  # keep zeros
  data_preprocessed[.naf(data == 0)] <- 0
  TE@isoAssays[[new_assay]] <- data_preprocessed

  # add results in assay format
  TE@qcAssays$nan <- is.nan(data.matrix(TE@isoAssays[[new_assay]]))
  TE@qcAssays$na <- is.na(data.matrix(TE@isoAssays[[new_assay]])) & !is.nan(data.matrix(TE@isoAssays[[new_assay]]))

  TE
}












#' Cleaning of data
#'
#' @param TE TracerExperiment
#' @param assay norm
#' @param qc_LOQ LOQ values
#' @param new_assay clean
#' @param thres_LOQ LOQ threshold (values below are set to NA)
#' @param max_nafrac_per_group maximum allowed fraction of NA values per group
#' @param min_rep_per_group minimum required number of non-NA values per group
#' @param split_by factor(s) from colData to define groups
#' @param exclude internal standard
#' @param remove_imp optionally, remove all imputed values
#' @param ...
#'
#' @export
#'
#' @examples
#' exampleTracerExperiment(add_qc = TRUE) |> estimateLOQs() |> clean(assay = "raw")
#' exampleTracerExperiment(add_qc = TRUE) |> estimateLOQs() |> sumMets(assay = "raw", sum_qc = TRUE) |> clean(assay = "raw", type = "met", qc_LOQ = "met_loq_raw")
clean <- function(TE, assay = "norm", qc_LOQ = "loq", new_assay = "clean", thres_LOQ = 1.5, max_nafrac_per_group = 1/3, min_rep_per_group = 2, split_by = ~ 1, type = "iso", exclude = "internal.standard", soft = NULL, remove_imp = FALSE, ...){

  if (is.null(soft)){ soft <- !is.null(thres_LOQ) }

  if (is.null(assay)) assay <- .getAssays(TE, assay = assay, last = TRUE, names_only = TRUE, type = type)
  data <- assay(TE, assay, type = type)
  data.orig <- data

  data[is.na(data)] <- NA
  design <- TE@colData
  if (is.null(data)) stop("Error: Data not found!")

  if (!is.null(qc_LOQ)){
    LOQ <- .getAssays(TE, assay = qc_LOQ, type = "qc")
    stopifnot(all.equal(dim(data), dim(LOQ)))
  }

  exclude <- intersect(exclude, rownames(data))
  if (length(exclude) == 0) exclude <- NULL

  # Optionally set imputed values to NA
  if (remove_imp == TRUE & !is.null(TE@qcAssays$na)){
    .colorcat("Setting imputed values to NA...")
    imp <- as.matrix(TE@qcAssays$na)
    data[imp] <- NaN
  }

  # Set values below LOQ to NA
  # NA ... only temporary
  if (!is.null(qc_LOQ)){
    data[LOQ < thres_LOQ] <- NA
  }

  # Split into groups
  if (length(labels(terms(split_by))) > 0){
    groups <- apply(design[,labels(terms(split_by)),drop = FALSE], 1, paste0,  collapse = "_")
    stopifnot(all(names(groups) == colnames(data)))
    data_grouped <- lapply(setNames(unique(groups), unique(groups)), function(g) data[,groups == g,drop = FALSE])
  } else {
    data_grouped <- list(data)
  }

  # Missing values per group
  # NaN ... newly removed
  if (!is.null(max_nafrac_per_group)){
    data_grouped <- lapply(data_grouped, function(tmp){
      tmp[rowMeans(is.na(tmp)) > max_nafrac_per_group & !rownames(tmp) %in% exclude,] <- NaN
      tmp
    })
  }
  if (!is.null(min_rep_per_group)){
    data_grouped <- lapply(data_grouped, function(tmp){
      tmp[rowSums(!is.na(tmp)) < min_rep_per_group & !rownames(tmp) %in% exclude,] <- NaN
      tmp
    })
  }

  # Combine again
  data_grouped <- lapply(data_grouped, function(tmp) tmp[rownames(data),])
  data_clean <- Reduce(cbind, data_grouped)
  data_clean <- data_clean[rownames(data), colnames(data)]

  if (soft == TRUE){
    set_na <- is.nan(data.matrix(data_clean))
    data.orig[set_na] <- NA
    data_clean <- data.orig
  } else {
    data_clean[is.na(data_clean)] <- NA
  }

  data_clean[is.nan(data.matrix(data.orig))] <- NaN

  if (is.null(new_assay)) return(data_clean)

  assay(TE, new_assay, type = type) <- data_clean
  TE
}



