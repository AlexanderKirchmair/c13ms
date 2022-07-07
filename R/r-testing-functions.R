

#' Make example TracerExperiment
#'
#' @param nsamples
#' @param nmets
#'
#' @return
#' @export
#'
#' @examples
exampleTracerExperiment <- function(nsamples = 10, nmets = 10, seed = 123){

  set.seed(seed)
  samples <- c(paste0(rep("A", floor(nsamples/2)), 1:floor(nsamples/2)), paste0(rep("B", ceiling(nsamples/2)), 1:ceiling(nsamples/2)))
  mets <- replicate(nmets, paste(sample(letters, sample(3:5, 1)), collapse = ""))
  mets[duplicated(mets)] <- paste0(mets[duplicated(mets)], 1:sum(duplicated(mets)))
  formulas <- replicate(nmets, rmolecule())
  stopifnot(sum(grepl("O0", formulas)) == 0)

  metData <- data.frame(row.names = mets, Molecule = formulas, MSion = NA)
  colData <- data.frame(row.names = samples, name = samples, group = sub("\\d", "", samples))

  nC <- as.numeric(sub("C", "", regmatches(formulas, regexpr("C[0-9]+", formulas))))
  nC <- setNames(nC, mets)
  isoData <- stack(lapply(nC, function(n) paste0("m", 0:n) ))[,c(2,1)]
  colnames(isoData) <- c("metabolite", "label")
  rownames(isoData) <- paste0(isoData$metabolite, "_", isoData$label)

  midl <- lapply(setNames(mets, mets), function(m){
    a <- sample(10:100, 1) * runif(min = 0.8, max = 1.2, n = nsamples)
    tmp <- rmid(ncol = nsamples, nrow = nrow(subset(isoData, metabolite == m)), abundance = a)
    rownames(tmp) <- paste(m, rownames(tmp), sep = "_")
    colnames(tmp) <- samples
    tmp
  })

  # add DE
  n <- nrow(subset(isoData, metabolite == mets[1]))
  nA <- sum(colData$group == "A")
  nB <- sum(colData$group == "B")
  tmpA <- rmid(p = (1:n)/n, ncol = nA, nrow = n, abundance = 100 * runif(min = 0.8, max = 1.2, n = nA))
  tmpB <- rmid(p = rev((1:n)/n), ncol = nB, nrow = n, abundance = 2 * runif(min = 0.8, max = 1.2, n = nB))
  tmp <- cbind(tmpA, tmpB)
  dimnames(tmp) <- dimnames(midl[[1]])
  midl[[1]] <- tmp

  data <- data.frame(isoData, Reduce(f = rbind, midl))

  makeTracerExperiment(data, metData = metData, colData = colData)
}





#' Generate random molecule formulas with n C atoms
#'
#' @param n
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
rmolecule <- function(n = NULL, ...){

  elements <- c("C", "H", "N", "O", "P", "S")
  p <- c(0.10, 0.62, 0.02, 0.24, 0.01, 0.01)

  if (is.null(n)) n <- sample(2:8, 1)

  nc <- 0
  while(nc != n){
    v <- table(sample(elements, n/p[1], prob = p, replace = TRUE))
    nc <- v["C"]
    if (is.null(nc) | is.na(nc)) nc <- 0
  }

  v <- setNames(as.character(v), names(v))
  v[v == "1"] <- ""
  v <- v[v != "0"]
  mol <- paste(paste0(names(v), v), collapse = "")

  mol
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
rmid <- function(nrow = NULL, ncol = 6, p = NULL, colnames = NULL, size = 10^3, abundance = NULL){

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

  if (!is.null(abundance)) mid <- t( t(mid) * abundance )
  as.matrix(mid)
}













