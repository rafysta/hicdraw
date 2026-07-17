#' @keywords internal
getRotate <- function(df) {
  df <- dplyr::mutate(df,
    x1 = (start1 + start2) / 2, x2 = (start1 + end2) / 2,
    x3 = (end1 + end2) / 2,
    y1 = (start2 - start1), y2 = end2 - start1, y4 = start2 - end1)
  df <- dplyr::mutate(df, x4 = x2, y3 = y1, g = dplyr::row_number())
  df <- dplyr::select(df, x1, x2, x3, x4, y1, y2, y3, y4, g, score)
  df.x <- tidyr::gather(dplyr::select(df, x1, x2, x3, x4, g, score),
                        key = "ord", value = "x", -g, -score)
  df.y <- tidyr::gather(dplyr::select(df, y1, y2, y3, y4, g, score),
                        key = "ord", value = "y", -g, -score)
  df.x <- dplyr::mutate(df.x,
    ord = as.integer(factor(ord, levels = c("x1", "x2", "x3", "x4"), ordered = TRUE)))
  df.y <- dplyr::mutate(df.y,
    ord = as.integer(factor(ord, levels = c("y1", "y2", "y3", "y4"), ordered = TRUE)))
  df <- dplyr::left_join(df.x, df.y, by = c("g", "score", "ord"))
  dplyr::arrange(df, g, ord)
}

#' Draw a Hi-C contact map as a rotated (45 degree) triangle
#'
#' Works with both `.hic` and `hic200` maps attached via [addHiCdata()].
#'
#' @param hic A Hi-C object.
#' @param dataName Name of the attached map.
#' @param col Fill colour, or `"matlab"` for the jet-like palette.
#' @param ymax Maximum distance from the diagonal to draw, in bp.
#' @param zmax Quantile used for the top of the colour scale.
#' @param xaxis Draw the x axis?
#' @param unit x-axis unit (e.g. 1e3 for kb).
#' @param ntick Approximate number of y ticks.
#' @param title,ylabel Plot title and strip label.
#' @param scale Numeric multiplier applied to contact counts (depth
#'   normalisation). Default 1.
#' @param rotate If `TRUE` (default) draw the 45-degree rotated triangle
#'   (x = position, y = distance). If `FALSE` draw the standard **square**
#'   matrix (x and y are both genomic position).
#' @param fixed Only used when `rotate = FALSE`. If `TRUE` (default) apply
#'   `coord_fixed(ratio = 1)` so the square map has a 1:1 x/y aspect. When such a
#'   square map is stacked with tracks, combine them with [drawStack()] /
#'   [draw_from_config()], which use `egg::ggarrange()` to keep the x-axis widths
#'   aligned (plain `cowplot`/`patchwork` cannot align fixed-aspect panels).
#'
#' @return A ggplot object.
#' @export
drawHiC <- function(hic, dataName = "hic1", col = "red", ymax = NULL,
                    zmax = 0.95, xaxis = TRUE, unit = 1e6, ntick = 3,
                    title = NULL, ylabel = "", scale = 1,
                    rotate = TRUE, fixed = TRUE) {
  rrn <- hic[["hic"]][[dataName]][["resolution"]]
  pchr <- hic[["drawing"]][["chr"]]
  pstart <- hic[["drawing"]][["start"]]
  pend <- hic[["drawing"]][["end"]]

  # -------- standard SQUARE matrix (x = y = position); coord_fixed by default --------
  if (!rotate) {
    unit_lab <- if (unit == 1e6) "Mb" else if (unit == 1e3) "kb" else "bp"
    raw <- .get_region_counts(hic, dataName, pchr, max(pstart, 1), pend, scale = scale)
    full <- rbind(
      data.frame(x = raw$x, y = raw$y, score = raw$counts),
      data.frame(x = raw$y, y = raw$x, score = raw$counts))
    full <- full[!duplicated(full[c("x", "y")]), , drop = FALSE]
    Tv <- as.numeric(stats::quantile(full$score, probs = zmax, na.rm = TRUE))
    df <- dplyr::mutate(full, start1 = x, end1 = x + rrn,
                        start2 = y, end2 = y + rrn, cate = ylabel)
    co <- if (fixed) {
      ggplot2::coord_fixed(ratio = 1, xlim = c(pstart / unit, pend / unit),
                           ylim = c(pend / unit, pstart / unit))
    } else {
      ggplot2::coord_cartesian(xlim = c(pstart / unit, pend / unit),
                               ylim = c(pend / unit, pstart / unit))
    }
    grob1 <- grid::grobTree(grid::textGrob("■", x = 0.06, y = 0.95, hjust = 1,
      vjust = 0.5, gp = grid::gpar(col = ifelse(col == "matlab", "red", col),
                                   fontsize = 18, fontface = "bold")))
    grob2 <- grid::grobTree(grid::textGrob(paste0(" ", round(Tv, 0)), x = 0.06,
      y = 0.95, hjust = 0, vjust = 0.5, gp = grid::gpar(col = "black", fontsize = 12)))
    p <- ggplot2::ggplot(df, ggplot2::aes(fill = score)) +
      ggplot2::geom_rect(ggplot2::aes(xmin = start1 / unit, xmax = end1 / unit,
                                      ymin = start2 / unit, ymax = end2 / unit)) +
      co +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::annotation_custom(grob1) + ggplot2::annotation_custom(grob2) +
      ggplot2::labs(x = paste0("Position in chr ", pchr, " (", unit_lab, ")"),
                    y = paste0("Position (", unit_lab, ")"), title = title) +
      ggplot2::facet_wrap(~cate, ncol = 1, strip.position = "right") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
        strip.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "black"))
    if (col == "matlab" && requireNamespace("colorRamps", quietly = TRUE)) {
      p <- p + ggplot2::scale_fill_gradientn(colours = colorRamps::matlab.like(10),
        labels = scales::scientific, limits = c(0, Tv), oob = scales::squish,
        name = "score", na.value = "white")
    } else {
      p <- p + ggplot2::scale_fill_gradient(name = "score", low = "white",
        high = ifelse(col == "matlab", "red", col), space = "Lab",
        limits = c(0, Tv), oob = scales::squish)
    }
    if (!xaxis) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                              axis.title.x = ggplot2::element_blank())
    }
    return(p)
  }
  # -------- rotated 45-degree triangle (default; unchanged) --------

  if (!is.null(ymax)) {
    ymax <- ymax / 1e3
  } else {
    ymax <- (pend - pstart) / 1e3
  }

  load_start <- max(pstart - ymax * 1e3, 1)
  load_end <- pend + ymax * 1e3
  raw <- .get_region_counts(hic, dataName, pchr, load_start, load_end, scale = scale)

  df <- getRotate(dplyr::mutate(raw,
    start1 = x, end1 = x + rrn, start2 = y, end2 = y + rrn, score = counts))
  Tv <- as.numeric(stats::quantile(df$score, probs = zmax, na.rm = TRUE))

  grob1 <- grid::grobTree(grid::textGrob("■", x = 0.06, y = 0.95,
    hjust = 1, vjust = 0.5,
    gp = grid::gpar(col = ifelse(col == "matlab", "red", col),
                    fontsize = 18, fontface = "bold")))
  grob2 <- grid::grobTree(grid::textGrob(paste0(" ", round(Tv, 0)),
    x = 0.06, y = 0.95, hjust = 0, vjust = 0.5,
    gp = grid::gpar(col = "black", fontsize = 12)))

  p <- ggplot2::ggplot(dplyr::mutate(df, cate = ylabel),
      ggplot2::aes(x = x / unit, y = y / 1e3, fill = score, group = g)) +
    ggplot2::geom_polygon() +
    ggplot2::coord_cartesian(xlim = c(pstart / unit, pend / unit),
                             ylim = c(0, ymax)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
      breaks = scales::pretty_breaks(n = ntick)) +
    ggplot2::annotation_custom(grob1) + ggplot2::annotation_custom(grob2) +
    ggplot2::labs(y = "Distance (kb)",
      x = paste0("Position in chr ", pchr, " (kb)"), title = title) +
    ggplot2::facet_wrap(~cate, ncol = 1, strip.position = "right") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black"))

  if (col == "matlab" && requireNamespace("colorRamps", quietly = TRUE)) {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = colorRamps::matlab.like(10), labels = scales::scientific,
      limits = c(0, Tv), oob = scales::squish, name = "score", na.value = "white")
  } else {
    p <- p + ggplot2::scale_fill_gradient(name = "score", low = "white",
      high = ifelse(col == "matlab", "red", col), space = "Lab",
      limits = c(0, Tv), oob = scales::squish)
  }

  if (!xaxis) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.title.x = ggplot2::element_blank())
  }
  p
}
