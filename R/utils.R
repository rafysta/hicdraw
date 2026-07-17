#' Remove inter-panel spacing so stacked tracks align cleanly
#'
#' Strips axis titles/margins from every plot except the last in a list, so a
#' vertical stack shares one x axis.
#'
#' @param plist A list of ggplot objects (top to bottom).
#' @return The modified list.
#' @export
removeSpace <- function(plist) {
  n <- length(plist)
  for (i in seq_len(n)) {
    plist[[i]] <- plist[[i]] +
      ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0.2, 0, 0.2), "cm"))
    if (i != n) {
      plist[[i]] <- plist[[i]] +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }
  }
  plist
}

#' Stack a Hi-C map with a boundary track (and optional gene track)
#'
#' Convenience wrapper that draws the map, its boundary-strength track and an
#' optional gene model, aligned vertically. Requires `cowplot`.
#'
#' @param hic A Hi-C object with a map and a boundary track attached.
#' @param hic_name Name of the attached map.
#' @param boundary_name Name of the attached boundary track (or `NULL`).
#' @param gene If `TRUE`, add the gene track (requires [addGenedata()]).
#' @param col,zmax,ymax,unit,title Passed to [drawHiC()].
#' @param heights Relative panel heights.
#' @param scale Depth-normalisation multiplier for the map.
#' @param rotate,fixed Passed to [drawHiC()]. Use `rotate = FALSE` for a square
#'   map; its x-axis is then kept aligned with the tracks via `egg::ggarrange()`.
#' @return A combined ggplot/cowplot object.
#' @export
drawStack <- function(hic, hic_name = "hic1", boundary_name = "boundary",
                      gene = FALSE, col = "matlab", zmax = 0.99, ymax = NULL,
                      unit = 1e3, title = NULL, heights = NULL, scale = 1,
                      rotate = TRUE, fixed = TRUE) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("package 'cowplot' is required for drawStack()")
  }
  plist <- list(drawHiC(hic, dataName = hic_name, col = col, zmax = zmax,
                        ymax = ymax, unit = unit, title = title,
                        xaxis = FALSE, scale = scale,
                        rotate = rotate, fixed = fixed))
  labs <- "hic"
  if (!is.null(boundary_name)) {
    plist <- c(plist, list(drawBW(hic, dataName = boundary_name,
      type = "BorderStrength", fill = "red", fill2 = "blue",
      xaxis = FALSE, xtitle = FALSE, ylabel = "BS", unit = unit)))
    labs <- c(labs, "bs")
  }
  if (gene) {
    plist <- c(plist, list(drawGene(hic, unit = unit, xaxis = TRUE)))
    labs <- c(labs, "gene")
  }
  plist <- removeSpace(plist)
  if (is.null(heights)) {
    heights <- c(4, rep(1.7, length(plist) - 2), if (gene) 2 else NULL)
    heights <- heights[seq_along(plist)]
  }
  # egg::ggarrange gives every panel the same plot-area width even when the map
  # uses coord_fixed (square). Fall back to cowplot if egg is not installed.
  if (requireNamespace("egg", quietly = TRUE)) {
    cowplot::plot_grid(egg::ggarrange(plots = plist, ncol = 1,
                                      heights = heights, draw = FALSE))
  } else {
    cowplot::plot_grid(plotlist = plist, ncol = 1, align = "v",
                       rel_heights = heights)
  }
}
