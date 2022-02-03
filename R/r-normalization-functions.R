
# NORMALIZATION FUNCTIONS



#' Normalization of isotopologue abundance data
#'
#' @param isoAssay
#' @param method
#' @param isodata
#' @param metAssay
#' @param fractions
#' @param colData
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
normalizeIsoAssay <- function(isoAssay, method = ~ IS + SUM, isodata = NULL, metAssay = NULL, fractions = NULL, colData = NULL, ...){

  ### Function for data normalization
  # provide custom functions in ellipsis
  # e.g. normalizeTE(data, method = ~ TMS + SUMFUN, SUMFUN = function(x){x/sum(x, na.rm = TRUE)})

  # Non-linear methods only work with summed metabolite abundances, e.g. log-based methods
  # These can be retrospectively weighted by the MID fractions to get isotopologue abundances

  ### INPUT ----

  data <- data.matrix(isoAssay)
  sumdata <- data.matrix(metAssay)

  args <- list(...)
  if (length(args) > 0){
    FUNs <- args[sapply(args, is.function)]
    FUNargs <- args[!sapply(args, is.function)]
  } else {
    FUNs <- NULL
  }

  methods <- strsplit(sub("~", "", deparse(method)), split = " + ", fixed = T)[[1]]
  methods.local <- methods[methods %in% c("COLSUM", "ISOSUM", "LOG", "TMS", "IS", "COLMEDIAN", "HKM", "ROWMEAN", "ROWMEDIAN", "TMM", "LEVEL", "AUTO", "SCALE", "VMN", "VSN")]
  methods.colData <- methods[!methods %in% methods.local & methods %in% colnames(colData)]
  methods.functions <- methods[!methods %in% methods.local & !methods %in% colnames(colData) & methods %in% names(FUNs)]

  # Apply methods in the order defined in the formula
  for (m in methods){

    if (m %in% methods.colData){

      # normalization factors from sample data
      cat("Normalization using", m, "\n")
      data <- normColData(data, colData[,m])

    } else if (m %in% methods.local){
      # apply locally defined functions
      m <- paste0("norm", m)
      cat("Normalization using", m, "\n")
      data <- do.call(m, c(list("data" = data.matrix(data), "sumdata" = sumdata, "fractions" = fractions, "isodata" = isodata, "coldata" = colData), list(...)))

    } else if (m %in% methods.functions){
      # apply user-defined functions
      cat("Normalization using", m, "\n")
      data <- FUNs[[m]](data, FUNargs)

    } else {
      stop(paste0("Method", m, " not found!"))
    }

  }


  stopifnot(dim(data) == dim(isoAssay))
  data
}











### FUNCTIONS ----

### Sample-wise linear transformations ----

# Simple normalization factors
normColData <- function(data, normfactor){
  stopifnot(length(normfactor) == ncol(data))
  t( t(data)/normfactor )
}


# Total metabolite sum
normCOLSUM <- function(data, ...){
  t( t(data)/colSums(data, na.rm = TRUE) )
}



normISOSUM <- function(data, isodata, FUN = colSums, ...){

  # Select consistently measured metabolites
  tmp <- data.frame(!is.na(data))
  tmp$metabolite <- isodata$metabolite
  namean <- .sumAssay(tmp, FUN = mean)
  mets <- setNames(matrixStats::rowMins(data.matrix(namean)), rownames(namean))
  mets <- sort(mets, decreasing = TRUE)

  data <- data.frame(data)
  data$metabolite <- isodata$metabolite
  imp_in <- data[data$metabolite %in% names(mets[mets>mean(mets)]),]
  imp_in$metabolite <- NULL
  imp <- missForest::missForest(imp_in)$ximp


  tmp <- data.frame(imp > 0)
  tmp$metabolite <- isodata[rownames(tmp),]$metabolite
  zeromean <- .sumAssay(tmp, FUN = mean)
  mets <- setNames(matrixStats::rowMins(data.matrix(zeromean)), rownames(zeromean))
  mets <- sort(mets, decreasing = TRUE)

  imp$metabolite <- isodata[rownames(imp),]$metabolite
  imp <- imp[imp$metabolite %in% names(mets[mets > mean(mets)]),]

  sumdata <- data.matrix(.sumAssay(imp))
  sizefactors <- FUN(sumdata)


  data$metabolite <- NULL
  normdata <- t( t(data)/ sizefactors )

  normdata
}




# Housekeeping metabolite
normHKM <- function(data, isodata, sumdata, ...){

  data <- data.frame(data)
  data$metabolite <- isodata$metabolite

  data$metabolite <- NULL

  tmp <- data[matrixStats::rowAlls(.naf(sumdata != 0)),]

  # Select housekeeping metabolites: low var, similar labelling
  vars <- matrixStats::rowVars(data.matrix(tmp)) / rowMeans(data.matrix(tmp), na.rm = TRUE)

  w <- colSums(data.matrix(tmp * 1/vars), na.rm = TRUE)
  w <- w / mean(w)

  normdata <- t( t(data.matrix(data))/ w )
  normdata[.naf(data == 0)] <- 0

  data.frame(normdata)
}






# Internal standard
normIS <- function(data, ...){
  ISmet <- list(...)[["ISmet"]]
  data <- data.matrix(data)
  t( t(data) / (mean(data[ISmet,])/data[ISmet,]) )
}

# Sample medians
normCOLMEDIAN <- function(sumdata, data, ...){
  sample_medians <- apply(sumdata, 2, median, na.rm = TRUE)
  t( t(data)/sample_medians )
}

# Trimmed mean of m-values
normTMM <- function(sumdata, data, ...){
  # ?????????????????????????????
  normfactors <- edgeR::calcNormFactors(sumdata[!is.na( rowSums(sumdata)),])
  t( t(data)/normfactors )
}





### Metabolite-wise transformation (i.e. weight metabolites by various scaling factors that are the same for all samples) ----


# Metabolite sum over samples
normROWMEAN <- function(data, ...){
  data / rowMeans(data, na.rm = TRUE)
}

normROWMEDIAN <- function(data, ...){
  data / matrixStats::rowMedians(data, na.rm = TRUE)
}



# Level-scaling
normLEVEL <- function(data, isodata, sumdata, ...){
  # Divide each metabolite by its mean over all samples
  data[rownames(isodata),] / rowMeans(sumdata, na.rm = TRUE)[isodata$metabolite]
}


# Auto-scaling
normSCALE <- function(sumdata, data, isodata, fractions, ...){

  scalefac <- matScale(sumdata, cols = TRUE, rows = FALSE, center = FALSE, scale = TRUE)
  fractions[rownames(isodata),] * scalefac[isodata$metabolite,]

}



### Non-linear transformations ----

# Log
normLOG <- function(sumdata, fractions, isodata, pseudocount = 1, ...){
  logsum <- log10(sumdata + pseudocount)
  fractions * logsum[isodata$metabolite,]
}



normVSN <- function(sumdata, fractions, isodata, ...){

  fitdata <- vsn::justvsn(x = data.matrix(sumdata))
  fractions[rownames(isodata),] * fitdata[isodata$metabolite,]

}







normVMN <- function(data, fractions, coldata, isodata, group, ...){

  norm <- vmn(data, groups = coldata[,group])
  fractions[rownames(isodata),] * norm[isodata$metabolite,]

}


vmn <- function(sumdata, groups, maxit = 10000, ...){

  data <- sumdata
  col_weights <- mean(colMeans(data, na.rm = TRUE)) / colMeans(data, na.rm = TRUE)
  row_weights <- rep(1, nrow(data))

  x <- setNames(c(col_weights, row_weights), NULL)
  f <- c(rep(TRUE, length(col_weights)), rep(FALSE, length(row_weights)))


  optfun <- function(par, data, return = FALSE){

    x <- exp(par)
    tmp <- data * x[f]
    tmp <- t( t(tmp) * x[!f] )
    if (return == TRUE) return(tmp)

    groups_uni <- unique(groups[!is.na(groups)])
    groups_uni <- setNames(groups_uni, as.character(groups_uni))
    means <- sapply(as.character(groups_uni), function(g){ rowMeans(tmp[,naf(groups == g)], na.rm = TRUE) })
    means[means == 0] <- NA

    # within-group
    var_min <- sapply(groups_uni, function(g) matrixStats::rowVars(tmp[,naf(groups == g)], na.rm = TRUE) / means[,g] )

    # between-group
    var_max <- matrixStats::rowVars(means, na.rm = TRUE) / rowMeans(means, na.rm = TRUE)

    val <- log(mean(var_min, na.rm = TRUE)) + 1 / log(mean(var_max, na.rm = TRUE))
    val
  }


  optres <- optim(par = log(x), fn = optfun, data = data, method = "SANN", control = list("maxit" = maxit))

  new <- optfun(optres$par, data = data, return = TRUE)
  new

}


















