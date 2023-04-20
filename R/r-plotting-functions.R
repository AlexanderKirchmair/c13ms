
# PLOTTING FUNCTIONS



#' variancePartition
#'
#' @param data
#' @param design
#' @param formula
#' @param rm.na
#'
#' @return
#' @export
#'
#' @examples
variancePartition <- function(data, design, formula = NULL, rm.na = FALSE){

  stopifnot(requireNamespace("variancePartition", quietly = TRUE))

  design <- design[, sapply(design, function(tmpcol) length(unique(tmpcol))>1 )]
  design <- design[, sapply(design, function(tmpcol) any(table(tmpcol) > 1) )]
  design <- design[!matrixStats::rowAnys(is.na(design)),]

  cols <- intersect(colnames(data), rownames(design))
  data <- data[,cols]
  design <- design[cols,]


  if (is.null(formula)){
    fac <- colnames(design)[sapply(design, function(tmp) is.factor(tmp) | is.character(tmp) )]
    formula <- paste(c(setdiff(colnames(design), fac), paste0("(1|", fac, ")")), collapse = " + ")
    formula <- as.formula(paste0("~ ", formula))
  }

  # Either remove NA features or set them to zero
  if (any(is.na(data))){
    if (rm.na == TRUE){
      message("Removing NA features...")
      data <- data[rowSums(is.na(data)) == 0,]
    } else {
      message("Setting NA values to zero...")
      data <- data[rowSums(!is.na(data)) >= 2,]
      data[is.na(data)] <- 0
    }
  }

  if (any(matrixStats::rowVars(data.matrix(data), na.rm = TRUE) == 0)){
    message("Removing zero-variance features...")
    data <- data[matrixStats::rowVars(data.matrix(data)) != 0,]
  }


  varpfit <- variancePartition::fitExtractVarPartModel(exprObj = data, formula = formula, data = design)
  varp <- variancePartition::sortCols(varpfit)
  gg1 <- variancePartition::plotPercentBars(varp)
  gg2 <- variancePartition::plotVarPart(varp)

  list("results" = varp, "plot_metabolites" = gg1, "plot_factors" = gg2)
}






#' RLE plots
#'
#' @param data
#' @param topN
#' @param title
#' @param cluster
#'
#' @return
#' @export
#'
#' @examples
ggrle <- function(data, topN = 10, title = NULL, cluster = TRUE){

  # Relative log expression

  rledata <- data.matrix(data)
  rledata <- rledata[rowSums(!is.na(rledata)) != 0,]

  min_val <- mean(sort(rledata[rledata>0])[1:10])/2
  min_val <- round(min_val, abs(round(log10(min_val))))
  rledata[rledata == 0] <- min_val

  logratios <- log10( rledata / matrixStats::rowMedians(rledata, na.rm = TRUE) )
  logratios[is.nan(logratios)] <- NA

  tmp <- logratios
  tmp[matrixStats::rowMedians(rledata, na.rm = TRUE) == 0,] <- sd(logratios[is.finite(logratios)], na.rm = TRUE)
  tmp[matrixStats::rowMedians(rledata, na.rm = TRUE) == 0,] <- tmp[matrixStats::rowMedians(rledata, na.rm = TRUE) == 0,] * sign(rledata[matrixStats::rowMedians(rledata, na.rm = TRUE) == 0,])

  tmp[tmp == -Inf] <- min(logratios[is.finite(logratios)], na.rm = T) - sd(logratios[is.finite(logratios)], na.rm = TRUE)
  tmp[tmp == -Inf] <- max(logratios[is.finite(logratios)], na.rm = T) + sd(logratios[is.finite(logratios)], na.rm = TRUE)


  ggdf <- reshape2::melt(tmp)
  colnames(ggdf) <- c("metabolite", "sample", "logratios")

  lab_mets <- names(sort(abs(setNames(ggdf$logratios, ggdf$metabolite)), decreasing = TRUE)[1:topN])

  if (cluster == TRUE) {
    tmp[is.na(tmp)] <- 0
    clust <- dendsort::dendsort(stats::hclust(stats::dist(t(tmp))))
    ddata <- ggdendro::dendro_data(as.dendrogram(clust), type = "rectangle")
    dendro <- ggplot2::ggplot() +
      cowplot::theme_nothing() +
      ggplot2::geom_text(data = ddata$labels, color = NA,  hjust = 0.5, vjust = 0, mapping = ggplot2::aes(x = label, y = 0, label = label)) + # this adds the correct y scale
      ggplot2::geom_segment(data = ddata$segments, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
      ggplot2::scale_y_continuous(expand = c(0,0))

    ggdf$sample <- factor(ggdf$sample, ordered = TRUE, levels = clust$labels[clust$order]) # reorder data according to the clustering

  }


  gg <- ggplot(ggdf, aes(x = sample, y = logratios)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_boxplot(fill = rgb(0.8, 0.8, 0.8), outlier.shape = NA) +
    geom_jitter(aes(color = metabolite), show.legend = TRUE) +
    scale_color_manual(breaks = lab_mets, values = circlize::rand_color(n = length(unique(ggdf$metabolite)), luminosity = "dark")) +
    xlab("") +
    ggtitle(title)


  if (cluster == TRUE) gg <- patchwork::wrap_plots(dendro, gg,  ncol = 1, heights = c(0.3,1))

  gg
}






#' Isotopologue plotting
#'
#' @param TE
#' @param assay
#' @param mets
#' @param summarise_by
#' @param dir
#' @param filename
#' @param title
#' @param title_size
#' @param order_by
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' exampleTracerExperiment(nmets = 2) |> MID("raw") |> isoplot(title = "Molecule", colorscale = c("white", "darkblue"), nacolor = "gray50", linewidth = 1, col_cex = 0.5)
isoplot <- function(TE, assay = "mid", mets = NULL, summarise_by = NULL, dir = NULL, filename = NULL, title = NULL, title_size = NULL, height = 10, order_by = NULL, dev = "png", ...){

  isodata <- isoData(TE)
  if (is.null(mets)) mets <- rownames(metData(TE))
  mids <- lapply(setNames(mets, mets), function(m) assay(TE, assay, met = m) )

  by <- rlang::enquo(summarise_by)
  if (!rlang::quo_is_null(by)){
    by <- rlang::as_name(by)
    coldata <- colData(TE)[,by, drop = FALSE]
    mids <- lapply(mids, function(mid){
      df <- data.frame(t(mid))
      df$rep <- coldata[rownames(df),, drop = TRUE]
      mid <- .sumAssay(df, var = rep, FUN = mean, na.rm = TRUE)
      data.frame(t(mid))
    })

  }

  order_by <- rlang::enquo(order_by)
  if (!rlang::quo_is_null(order_by)){
    order <- colData(TE) %>% dplyr::arrange(!!order_by) %>% rownames()
    mids <- lapply(mids, function(mid){ mid[,order, drop = FALSE] })
  }


  gg <- lapply(mids, function(mid){ ggiso(mid = mid, isodata = isodata, title_size = title_size, ...) })

  if (!is.null(title)){
    gg <- lapply(setNames(names(gg), names(gg)), function(m){
      if (!is.null(gg[[m]])){
        gg[[m]] + ggplot2::ggtitle(metData(TE)[m,title])
      } else {
        NULL
      }

    })
  }


  if (is.null(dir) & length(gg) == 1) return(gg[[1]])
  if (is.null(dir) & length(gg) > 1) return(gg)

  if (!dir.exists(dir)) dir.create(dir)


  p <- invisible(lapply(setNames(names(gg), names(gg)), function(m){
    if (is.null(gg[[m]])) return(NULL)
    # nc <- ncol(mids[[m]])
    # nr <- nrow(mids[[m]])
    #
    # tmp <- gg[[m]] + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none")
    # gtmp <- ggplotGrob(tmp)
    # wtmp <- sum(grid::convertWidth(gtmp$widths, "in", TRUE))
    # htmp <- sum(grid::convertHeight(gtmp$heights, "in", TRUE))
    #
    # g <- ggplotGrob(gg[[m]])
    # w <- sum(grid::convertWidth(g$widths, "in", TRUE))
    # h <- sum(grid::convertHeight(g$heights, "in", TRUE))
    #
    # w0 <- w - wtmp
    # h0 <- h - htmp

    #wtmp/nc
    #htmp/nr
    d <- DeLuciatoR::get_dims(gg[[m]], maxheight = height)
    w <- d$width / d$height

    if (dev == "png"){
      ggplot2::ggsave(plot = gg[[m]],
                      filename = file.path(dir, paste0(m, ".", dev)),
                      device = dev, type = "cairo", dpi = 300,
                      width = height * w,
                      height = height)
    } else {
      ggplot2::ggsave(plot = gg[[m]],
                      filename = file.path(dir, paste0(m, ".", dev)),
                      device = dev, dpi = 300,
                      width = height * w,
                      height = height)
    }


  }))

  invisible(NULL)
}








#' Isotopologue plotting
#'
#' @param mid
#' @param isodata
#' @param label
#' @param cumulative
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggiso <- function(mid, isodata, label = TRUE, cumulative = FALSE, title_size = NULL, ...){

  if (all(is.na(mid))) return(NULL)

  mid.org <- mid
  if (cumulative == TRUE){
    tmp <- mid
    tmp[is.na(tmp)] <- 0
    tmp <- apply(tmp, 2, function(x) rev(cumsum(rev(x))))
    dim(tmp) <- dim(mid)
    dimnames(tmp) <- dimnames(mid)
    tmp[is.na(mid)] <- NA
    mid <- tmp
  }

  df <- .stretch(data.frame(mid))
  df.org <- .stretch(data.frame(mid.org))
  df$metabolite <- isodata[df$Feature,]$metabolite
  df$iso <- isodata[df$Feature,]$label
  df$Value <- df$Value * 100
  df$label <- paste0(round(df.org$Value * 100), "%")
  df$Sample <- factor(df$Sample, ordered = TRUE, levels = colnames(mid))

  if (cumulative == TRUE){
    df <- droplevels(subset(df, iso != "m0"))
    if (nrow(df) == 0) return(NULL)
  }

  if (label == TRUE){
    mapping <- ggplot2::aes(x = iso, y = Sample, fill = Value, label = label)
  } else {
    mapping <- ggplot2::aes(x = iso, y = Sample, fill = Value)
  }


  gg <- ggcircles(data = droplevels(df), mapping = mapping, title_size = title_size, ...) + ggplot2::labs(fill = "label %")

  gg
}




#' Plot circles with ggplot
#'
#' @param data
#' @param mapping
#' @param sym
#' @param r
#' @param fontsize
#' @param title_size
#' @param labelsize
#' @param labelface
#' @param labelcolor
#' @param col_cex
#' @param colorscale
#' @param nacolor
#' @param legend
#' @param linewidth
#' @param plot_margin
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' data.frame(x = rep(1:3,3), y = rep(c("A","B","C"), each=3), value = runif(9)*100) |> ggcircles()
ggcircles <- function(data, mapping = aes(x = x, y = y, fill = value), sym = TRUE,
                      r = 0.5, fontsize = 10, title_size = 10, labelsize = 10, labelface = 1, labelcolor = "black", col_cex = 1,
                      colorscale = c("white", "mediumblue"), nacolor = "gray90", legend = TRUE, linewidth = 0.5, plot_margin = 1, ...){

  ### Input ----

  xcol <- rlang::as_name(mapping[["x"]])
  ycol <- rlang::as_name(mapping[["y"]])

  if (!is.factor(data[[xcol]])) data[[xcol]] <- factor(data[[xcol]], ordered = TRUE)
  if (!is.factor(data[[ycol]])) data[[ycol]] <- factor(data[[ycol]], ordered = TRUE)

  data[[ycol]] <- revOrder(data[[ycol]])

  data$.x <- as.numeric(data[[xcol]])
  data$.y <- as.numeric(data[[ycol]]) * 1.1
  data$.radius <- r

  xbreaks <- setNames(data$.x, data[[xcol]])
  xbreaks <- xbreaks[!base::duplicated(xbreaks)]

  ybreaks <- setNames(data$.y, data[[ycol]])
  ybreaks <- ybreaks[!base::duplicated(ybreaks)]


  ### Mappings ----

  # global
  aes1 <- ggplot2::aes(x = .x, y = .y)

  # geom_circle
  aes2 <- aes1
  aes2["r"] <- ggplot2::aes(r = .radius)
  tmp <- mapping[names(mapping) %in% c("fill", "r")]
  aes2[names(tmp)] <- tmp
  names(aes2)[names(aes2) == "x"] <- "x0"
  names(aes2)[names(aes2) == "y"] <- "y0"

  # geom_text
  aes3 <- aes1
  aes3["label"] <- mapping["label"]


  ### GGPLOT ----

  gg <- ggplot2::ggplot(data = data, mapping = aes1) +
    theme_circles(base_size = fontsize, title_size = title_size, col_cex = col_cex, plot_margin = plot_margin) +
    ggforce::geom_circle(aes2, size = linewidth, inherit.aes = FALSE, ...)

  gg %<>% + ggplot2::scale_x_continuous(breaks = xbreaks,
                                        expand = ggplot2::expansion(mult = c(0,0)),
                                        limits = range(data$.x) + c(-0.5, 0.5),
                                        labels = names(xbreaks),
                                        position = "top")

  if (sym == TRUE){
    gg %<>% + ggplot2::scale_y_continuous(breaks = ybreaks,
                                          expand = ggplot2::expansion(mult = c(0,0)),
                                          limits = range(data$.y) + c(-0.5, 0.5),
                                          labels = names(ybreaks),
                                          sec.axis = ggplot2::sec_axis(~., breaks = ybreaks, labels = names(ybreaks)))
  } else {
    gg %<>% + ggplot2::scale_y_continuous(breaks = ybreaks,
                                          expand = ggplot2::expansion(mult = c(0,0)),
                                          limits = range(data$.y) + c(-0.5, 0.5),
                                          labels = names(ybreaks))
  }


  gg %<>% + ggplot2::scale_fill_gradientn(colours = colorscale,
                                          limits = c(0,100),
                                          na.value = nacolor,
                                          values = (seq_along(colorscale)-1)/(length(colorscale)-1),
                                          guide = ifelse(legend, "colourbar", "none"))

  gg <- gg + ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::coord_fixed(clip = "off")


  if ("label" %in% names(mapping)) gg <- gg + ggplot2::geom_text(aes3, size = labelsize/ggplot2:::.pt, fontface = labelface, color = labelcolor)

  gg
}





revOrder <- function(x){
  if (!is.factor(x)){
    x <- rev(x)
  } else {
    x <- factor(x, ordered = TRUE, levels = rev(levels(x)))
  }
  x
}




theme_circles <- function(base_size = 10, title_size = 10, col_cex = 1, plot_margin = 1, base_family = "", base_line_size = base_size/80, base_rect_size = base_size/80){

  th1 <- ggplot2::theme_minimal(base_size = base_size,
                                base_family = base_family,
                                base_line_size = base_line_size)

  th2 <- ggplot2::theme(panel.spacing = grid::unit(c(0,0,0,0), units = "mm"),
                        plot.margin = grid::unit(c(1,1,1,1)*plot_margin, "mm"),
                        plot.title = ggplot2::element_text(size = title_size,
                                                  color = rgb(0,0,0),
                                                  face = "bold",
                                                  angle = 0,
                                                  hjust = 0.5,
                                                  vjust = 1,
                                                  margin = ggplot2::margin(c(0,0,0,0)) ),
                        axis.title = ggplot2::element_blank(),
                        axis.text = ggplot2::element_text(margin = ggplot2::margin(c(0,0,0,0)),
                                                 color = rgb(0,0,0),
                                                 size = ggplot2::rel(1)),
                        axis.text.x = ggplot2::element_text(size = base_size * col_cex, hjust = 0.5, vjust = 1),
                        axis.text.y.left = ggplot2::element_text(size = base_size, hjust = 1, vjust = 0.5),
                        axis.text.y.right = ggplot2::element_text(colour = NA, size = base_size, hjust = 1, vjust = 0.5), # symmetric
                        panel.grid.major = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        legend.margin = ggplot2::margin(c(0,0,0,0), unit = "mm"),
                        legend.justification = "center",
                        legend.text = ggplot2::element_text(size = base_size,
                                                   margin = ggplot2::margin(c(0,0,0,0)),
                                                   color = rgb(0,0,0)),
                        legend.title = ggplot2::element_text(size = base_size,
                                                     margin = ggplot2::margin(c(0,0,0,0)),
                                                     color = rgb(0,0,0)),
                        complete = TRUE)

  ggplot2::`%+replace%`(th1, th2)

}





