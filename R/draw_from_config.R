#' Generate many figures from a configuration table
#'
#' Reads a configuration table (one row per output figure) and a sample sheet,
#' then draws and saves each requested figure. Modelled on the parameter-file
#' driven workflow used for the .hic maps, but aware of the hic200 format and
#' of boundary tracks.
#'
#' Configuration columns (missing optional columns are tolerated):
#' \describe{
#'   \item{id}{Row identifier.}
#'   \item{type}{`"original"` (default) or `"subtraction"` (map1 - map2).}
#'   \item{data_format}{`"hic200"` (default) or `"hic"`.}
#'   \item{input_id1, input_id2}{`sample_id` values referencing the sample sheet.}
#'   \item{normalization}{Only used for `.hic` maps.}
#'   \item{resolution}{e.g. `"200bp"`.}
#'   \item{chr, start, end}{Region to draw.}
#'   \item{max_distance}{Distance from the diagonal to draw (bp); default region width.}
#'   \item{color_scale_max_value}{Explicit colour-scale maximum.}
#'   \item{color_scale_max_percentile}{Quantile for the colour-scale maximum (default 0.99).}
#'   \item{with_boundary}{`TRUE`/`FALSE` — stack a boundary-strength panel.}
#'   \item{scale}{`"none"` (default) or `"depth"` (normalise by total_read).}
#'   \item{output_name}{Output file stem.}
#' }
#'
#' @param config_file Path to a `.xlsx`, `.tsv` or `.csv` configuration table.
#' @param sample_sheet Path to the sample sheet (`.tsv`) or a data.frame with
#'   columns `sample_id`, `matrix_file`, `boundary_file`, `total_read`.
#' @param matrix_dir Directory holding the matrix files.
#' @param boundary_dir Directory holding the boundary files.
#' @param out_dir Output directory (created if needed).
#' @param gene_file Optional gene model file to add a gene panel.
#' @param chrom_sizes Chromosome sizes (default [pombe_chrom_sizes]).
#' @param depth_target Reference total read count for `scale = "depth"`
#'   (default 1e8).
#' @param base_width,base_height Saved figure size in cm.
#' @param dpi Saved figure resolution.
#' @param rows Optional integer vector to restrict which config rows are drawn.
#'
#' @return Invisibly, a character vector of written file paths.
#' @export
draw_from_config <- function(config_file, sample_sheet, matrix_dir,
                             boundary_dir, out_dir, gene_file = NULL,
                             chrom_sizes = pombe_chrom_sizes,
                             depth_target = 1e8, base_width = 22,
                             base_height = 16, dpi = 200, rows = NULL) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("package 'cowplot' is required for draw_from_config()")
  }
  cfg <- .read_table(config_file)
  ss <- if (is.data.frame(sample_sheet)) sample_sheet else .read_table(sample_sheet)
  ss$sample_id <- as.character(ss$sample_id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (is.null(rows)) rows <- seq_len(nrow(cfg))

  written <- character(0)
  for (i in rows) {
    row <- cfg[i, , drop = FALSE]
    type <- .col(row, "type", "original")
    fmt <- .col(row, "data_format", "hic200")
    res <- as.character(.col(row, "resolution", "200bp"))
    chr <- as.character(row[["chr"]])
    pstart <- as.numeric(row[["start"]])
    pend <- as.numeric(row[["end"]])
    out_name <- as.character(.col(row, "output_name", paste0("fig_", i)))
    with_bnd <- isTRUE(as.logical(.col(row, "with_boundary", TRUE)))
    scale_mode <- as.character(.col(row, "scale", "none"))
    max_dist <- suppressWarnings(as.numeric(.col(row, "max_distance", pend - pstart)))
    if (is.na(max_dist)) max_dist <- pend - pstart

    s1 <- ss[ss$sample_id == as.character(row[["input_id1"]]), , drop = FALSE]
    if (nrow(s1) == 0) { message("row ", i, ": input_id1 not in sample sheet"); next }
    f1 <- file.path(matrix_dir, s1$matrix_file[1])
    tr1 <- suppressWarnings(as.numeric(s1$total_read[1]))
    sc1 <- if (identical(scale_mode, "depth") && !is.na(tr1)) depth_target / tr1 else 1

    hic <- createHiCmap()
    if (!is.null(gene_file)) hic <- addGenedata(hic, gene_file)
    hic <- addHiCdata(hic, f1, resolution = res, dataName = "m1",
                      format = fmt, chrom_sizes = chrom_sizes, total_read = tr1)

    if (with_bnd && !is.na(s1$boundary_file[1])) {
      hic <- addBoundary(hic, file.path(boundary_dir, s1$boundary_file[1]),
                         dataName = "b1")
    }
    hic <- setRegion(hic, chr, pstart, pend)

    zmax <- .col(row, "color_scale_max_percentile", 0.99)
    zmax <- suppressWarnings(as.numeric(zmax)); if (is.na(zmax)) zmax <- 0.99

    p_hic <- drawHiC(hic, dataName = "m1", col = "matlab", zmax = zmax,
                     ymax = max_dist, unit = 1e3, title = out_name,
                     xaxis = FALSE, scale = sc1)
    plist <- list(p_hic)
    if (with_bnd && !is.null(hic[["bw"]][["b1"]])) {
      plist <- c(plist, list(drawBW(hic, dataName = "b1",
        type = "BorderStrength", fill = "red", fill2 = "blue",
        xaxis = !(!is.null(gene_file)), xtitle = FALSE, ylabel = "BS", unit = 1e3)))
    }
    if (!is.null(gene_file)) {
      plist <- c(plist, list(drawGene(hic, unit = 1e3, xaxis = TRUE)))
    }
    plist <- removeSpace(plist)
    heights <- c(4, rep(1.7, length(plist) - 1))
    # egg::ggarrange keeps panel widths aligned even for a coord_fixed (square)
    # map; fall back to cowplot if egg is unavailable.
    if (requireNamespace("egg", quietly = TRUE)) {
      pmix <- cowplot::plot_grid(egg::ggarrange(plots = plist, ncol = 1,
                                                heights = heights, draw = FALSE))
    } else {
      pmix <- cowplot::plot_grid(plotlist = plist, ncol = 1, align = "v",
                                 rel_heights = heights)
    }
    fout <- file.path(out_dir, paste0(out_name, ".png"))
    cowplot::save_plot(fout, pmix, base_width = base_width / 2.54,
                       base_height = base_height / 2.54, dpi = dpi)
    written <- c(written, fout)
    message(sprintf("[%d/%d] %s", i, length(rows), fout))
  }
  invisible(written)
}

#' @keywords internal
.read_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "xlsx") {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("package 'readxl' is required to read .xlsx config")
    }
    as.data.frame(readxl::read_excel(path))
  } else {
    as.data.frame(data.table::fread(path))
  }
}

#' @keywords internal
.col <- function(row, name, default) {
  if (name %in% colnames(row) && !is.na(row[[name]])) row[[name]] else default
}
