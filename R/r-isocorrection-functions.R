
# NATURAL ISOTOPE ABUNDANCE CORRECTION FUNCTIONS





#' Natural isotope abundance correction
#'
#' @param TE
#' @param assay
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' exampleTracerExperiment() %>% correctIso(assay = "raw")
correctIso <- function(TE, assay = "imp", ...){

  assay <- .getAssays(TE, assay = assay, type = "iso")
  metdata <- dplyr::select(TE@metData, c(Molecule, MSion))

  res <- isoCorr(assaydata = assay, molecules = metdata, isocomp = isocomp, ...)

  TE@isoAssays$corr <- res
  TE
}





#' Wrapper function for IsoCorrector
#'
#' @param assaydata
#' @param molecules
#' @param Molecule
#' @param MSion
#' @param tracer
#' @param tracer_purity
#' @param mode
#' @param tol_error
#' @param highres
#' @param tmpdir
#' @param allresults
#'
#' @return
#' @export
#'
#' @examples
isoCorr <- function(assaydata, molecules, Molecule = Molecule, MSion = MSion, tracer = "C", tracer_purity = NA, mode = "negative", tol_error = 0.05, highres = TRUE, tmpdir = NULL, isocomp = NULL, allresults = FALSE){

  stopifnot(requireNamespace("IsoCorrectoR"))

  .colorcat("Natural isotope abundance correction with IsoCorrectoR", col = "#77a4c7")
  .colorcat("https://doi.org/10.1038/s41598-018-36293-4", col = "#77a4c7")


  # Set tmp dir
  if (is.null(tmpdir)) tmpdir <- "./IsoCorrector"
  if (!dir.exists(tmpdir)){
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  }

  Molecule <- rlang::enquo(Molecule)
  MSion <- rlang::enquo(MSion)


  ### Prepare input files ----

  ### 1. Measurement file
  message("Preparing measurement file...")
  MeasurementFile <- file.path(tmpdir, "MeasurementFile.csv")
  measurement <- data.frame("Measurements/Samples" = gsub(pattern = "_m", replacement = paste0("_", tracer), x = rownames(assaydata)),
                            assaydata, check.names = FALSE)
  write.csv(measurement, file = MeasurementFile, row.names = FALSE, na = "")


  ### 2. Molecule information file
  # Three columns: "Molecule", "MS ion or MS/MS product ion" and "MS/MS neutral loss"
  message("Preparing molecule file...")
  MoleculeFile <- file.path(tmpdir, "MoleculeFile.csv")
  MolTable <- .makeMoleculeTable(molecules = molecules, Molecule = Molecule, MSion = MSion, mode = mode, tracer = tracer)
  write.csv(MolTable, file = MoleculeFile, row.names = FALSE, quote = FALSE)

  ### 3. Element information file
  # Four columns: "Element", "Isotope abundance_Mass shift", "Tracer isotope mass shift", and "Tracer purity"
  message("Preparing element file...")
  ElementFile <- file.path(tmpdir, "ElementFile.csv")
  if (is.null(isocomp)) isocomp <- getElementCompositions()
  ElemTable <- .makeElementsTable(df = MolTable, isocomp = isocomp, tracer = tracer, tracer_purity = tracer_purity)
  write.table(ElemTable, file = ElementFile, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)


  ### Run IsoCorrectoR ----

  icresults <- IsoCorrectoR::IsoCorrection(MeasurementFile = MeasurementFile,
                                           MoleculeFile = MoleculeFile,
                                           ElementFile = ElementFile,
                                           CorrectTracerImpurity = !is.na(tracer_purity),
                                           CorrectTracerElementCore = TRUE,
                                           CalculateMeanEnrichment = TRUE,
                                           UltraHighRes = highres,
                                           DirOut = tmpdir,
                                           FileOut = "Results",
                                           FileOutFormat = "csv",
                                           ReturnResultsObject = TRUE,
                                           CorrectAlsoMonoisotopic = FALSE,
                                           CalculationThreshold = 10^-8,
                                           CalculationThreshold_UHR = 8,
                                           verbose = FALSE,
                                           Testmode = FALSE)


  if (icresults$success == FALSE){
    message("Error in IsoCorrectoR::IsoCorrection!")
  } else if (allresults == TRUE){
    return(icresults)

  } else {


    # low values can be kept/replaced by zero

    resid <- icresults$results$Residuals
    rownames(resid) <- sub(paste0("_", tracer),"_m", rownames(resid))
    resid <- resid[rownames(assaydata),]

    res <- icresults$results$Corrected
    rownames(res) <- sub(paste0("_", tracer),"_m", rownames(res))
    res <- res[rownames(assaydata),]

    relerror <- icresults$results$RelativeResiduals
    rownames(relerror) <- sub(paste0("_", tracer),"_m", rownames(relerror))
    relerror <- relerror[rownames(assaydata),]

    # error relative to isotopologue
    relerror[is.na(relerror)] <- 0
    relerror[!is.finite(data.matrix(relerror))] <- 0
    relerror[.naf(assaydata == 0)] <- NA
    relerror[.naf(res == 0)] <- NA

    # error relative to metabolite
    resid[.naf(assaydata == 0)] <- NA
    resid[.naf(res == 0)] <- NA
    resid$met <- gsub("_m.*$", "", rownames(resid))
    sumdata <- resid %>% dplyr::group_by(met) %>% dplyr::summarise(dplyr::across(.fns = sum, na.rm = TRUE))
    sumdata <- as.data.frame(sumdata[match(resid$met, sumdata$met),])
    sumdata$met <- NULL
    rownames(sumdata) <- rownames(resid)
    resid$met <- NULL
    resid_rel_total <- resid / sumdata

    if (is.null(tol_error)) tol_error <- NA
    if (!is.na(tol_error)) res[.naf(abs(relerror) > tol_error) & .naf(abs(resid_rel_total) > tol_error)] <- NA

    res <- lapply(setNames(colnames(assaydata), colnames(assaydata)), function(tmp){
      if (!is.null(res[[tmp]])){ res[[tmp]] } else { rep(NA, nrow(res))}
    })
    res <- as.data.frame(res)
    dimnames(res) <- dimnames(assaydata)
    res[is.na(res) & assaydata == 0] <- 0
    res[is.na(assaydata)] <- NA
    res[is.nan(data.matrix(assaydata))] <- NaN
    stopifnot(all(dim(res) == dim(assaydata)))

    return(res)
  }

}




#' Internal function to generate molecule table for IsoCorrector
#'
#' @param molecules
#' @param Molecule
#' @param MSion
#' @param mode
#' @param tracer
#'
#' @return
#'
#' @examples
.makeMoleculeTable <- function(molecules, Molecule, MSion, mode, tracer){

  df <- data.frame("Molecule" = rownames(molecules),
                   "MS ion or MS/MS product ion" = dplyr::pull(.data = molecules, !!MSion),
                   "MS/MS neutral loss" = "", stringsAsFactors = FALSE, check.names = FALSE)

  # For molecules without ions, automatically derive ion
  chargeH <- ifelse(mode == "negative", -1, 1)
  msions <- sapply(molecules[[rlang::as_name(Molecule)]], FUN = function(X){
    if (!is.na(X)){
      X <- .add1(X)
      Hmol <- regmatches(X, regexpr("H[0-9]+", X))
      nH <- as.numeric(sub("H","", Hmol))
      Hion <- paste("H", nH + chargeH, sep = "")
      x <- sub(pattern = Hmol, replacement = Hion, x = X )
      gsub("H0", "", x)
    } else {
      NA
    }
  })
  df[[2]][is.na(df[[2]])] <- msions[is.na(df[[2]])]

  # Put formulas into correct format
  df[[2]] <- sapply(df[[2]], .add1)
  traceratoms <- regmatches(df[[2]], regexpr(paste0(tracer, "[0-9]+"), df[[2]]))
  df[[2]] <- paste(df[[2]], paste("Lab", traceratoms, sep = ""), sep = "")

  df
}



#' Internal function to generate elements table for IsoCorrector
#'
#' @param df
#' @param isocomp
#' @param tracer
#' @param tracer_purity
#'
#' @return
#'
#' @examples
.makeElementsTable <- function(df, isocomp, tracer, tracer_purity){

  all_elements <- unique(sub(unlist(strsplit(df[[2]], split = "[0-9]+")), pattern = "Lab", replacement = ""))
  elements <- data.frame(row.names = all_elements, matrix(nrow = length(all_elements), ncol = 4))
  colnames(elements) <- c("Element", "Isotope abundance_Mass shift", "Tracer isotope mass shift", "Tracer purity")
  elements$Element <- all_elements

  elements$`Isotope abundance_Mass shift` <- sapply(isocomp[elements$Element], function(tmp){
    tmp_use <- tmp[tmp$`Isotopic Composition` > 0,]
    ref <- which(tmp_use$`Isotopic Composition` == max(tmp_use$`Isotopic Composition`))
    tmp_use$shift <- tmp_use$`Mass Number` - tmp_use$`Mass Number`[ref]
    paste(paste(format(tmp_use$`Isotopic Composition`, scientific = FALSE), tmp_use$shift, sep = "_"), collapse = "/")
  })

  elements$`Tracer isotope mass shift` <- ifelse(elements$Element == tracer, "1", "") # only high-res
  elements$`Tracer purity` <- ifelse(elements$Element == tracer, tracer_purity, "") #

  elements.write <- c(paste(colnames(elements), collapse = ","),
                      gsub(",,", "", apply(elements, 1, function(tmp){paste(tmp, collapse = ",")})))

  elements.write
}





#' Retrieve element isotope compositions from NIST
#'
#' @return
#' @export
#'
#' @examples
getElementCompositions <- function(){

  file <- rvest::read_html(("https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=all"))
  rawtable <-  strsplit(x = rvest::html_text(file)[[1]], split = "\n")[[1]]

  start <- grep("Description of Quantities and Notes", rawtable)
  end <- grep("NIST | Physical Measurement Laboratory | Physical Reference Data | Atomic Weights and Isotopic Compositions Main Page", rawtable)
  cleantable <- rawtable[(start+2):(end-2)]

  isolist <- split(cleantable, f = cumsum(nchar(cleantable) == 0))
  isolist <- lapply(isolist, function(tmp) tmp[nchar(tmp) != 0])
  isolist <- lapply(isolist, function(tmpslot) setNames(sapply(strsplit(tmpslot, " = "), "[", 2), sapply(strsplit(tmpslot, " = "), "[", 1)) )

  names(isolist) <- sapply(isolist, function(tmp) sub("Atomic Number = ", "", tmp[[1]]) )


  isocomp <- lapply(unique(names(isolist)), function(n){

    df <- as.data.frame(t(as.data.frame(isolist[names(isolist) == n])))
    df[["Atomic Number"]] <- as.numeric(df[["Atomic Number"]])
    df[["Mass Number"]] <- as.numeric(df[["Mass Number"]])
    df[["Relative Atomic Mass"]] <- as.numeric(gsub("\\(.*\\)", "", df[["Relative Atomic Mass"]]))
    df[["Isotopic Composition"]] <- as.numeric(gsub("\\(.*\\)", "", df[["Isotopic Composition"]]))
    df[["Isotopic Composition"]][is.na(df[["Isotopic Composition"]])] <- 0

    df[["Standard Atomic Weight"]] <- lapply(strsplit(gsub("\\[|\\]", "", df[["Standard Atomic Weight"]]), split = ","), "c")
    suppressWarnings({ df[["Standard Atomic Weight"]] <- lapply(df[["Standard Atomic Weight"]], as.numeric) })
    df[["Notes"]] <- NULL

    rownames(df) <- paste0(df$`Atomic Symbol`[1], df$`Mass Number`)
    return(df)
  })

  names(isocomp) <- sapply(isocomp, function(df) df$`Atomic Symbol`[1] )
  isocomp
}




#' Helper function for formatting molecule formulas for IsoCorrector
#'
#' @param tmp
#'
#' @return
#'
#' @examples
.add1 <- function(tmp){
  # split string and get elements
  tmp.el <- unlist(strsplit(tmp, "(?<=.)(?=[[:upper:]])", perl = T)) # split at each uppercase letter
  tmp.nel <- gsub("[^[:digit:]]", "", tmp.el)
  tmp.nel[nchar(tmp.nel) == 0] <- 1
  # put back together
  tmp.ellett <- gsub("[[:digit:]]", "", tmp.el)
  tmp.res <- paste0(tmp.ellett, tmp.nel)
  res <- paste(tmp.res, collapse = "")
  res
}






