#' Attach a Hi-C contact map to the object
#'
#' Two input formats are supported:
#' * `format = "hic"` — a Juicer `.hic` file, read on demand through
#'   [strawr::straw()] (the classic behaviour).
#' * `format = "hic200"` — a sparse `bin1 bin2 score` matrix produced by
#'   hic200-cpp. The file is read once with [read_hic200()] and cached inside
#'   the object.
#'
#' @param hic A Hi-C object.
#' @param file Path to the map file.
#' @param resolution Resolution string, e.g. `"200bp"`, `"5kb"`.
#' @param dataName Name used to refer to this map later.
#' @param normalization Normalisation for `.hic` (`"ICE"` is mapped to
#'   `"NONE"`); ignored for `hic200`.
#' @param format `"hic"` or `"hic200"`.
#' @param chrom_sizes Chromosome sizes for `hic200` (default [pombe_chrom_sizes]).
#' @param total_read Optional total read count for depth normalisation; stored
#'   for later use by [drawHiC()] / [draw_from_config()].
#'
#' @return The updated Hi-C object.
#' @export
addHiCdata <- function(hic, file, resolution, dataName = "hic1",
                       normalization = "ICE", format = c("hic", "hic200"),
                       chrom_sizes = pombe_chrom_sizes, total_read = NA_real_) {
  format <- match.arg(format)
  if (!file.exists(file)) {
    message(file, " is not found")
    return(hic)
  }

  res_num <- .parse_resolution(resolution)
  hic[["hic"]][[dataName]][["file"]] <- file
  hic[["hic"]][[dataName]][["resolution_string"]] <- resolution
  hic[["hic"]][[dataName]][["resolution"]] <- res_num
  hic[["hic"]][[dataName]][["format"]] <- format
  hic[["hic"]][[dataName]][["total_read"]] <- total_read

  if (format == "hic") {
    if (normalization == "ICE") normalization <- "NONE"
    hic[["hic"]][[dataName]][["normalization"]] <- normalization
  } else {
    hic[["hic"]][[dataName]][["normalization"]] <- "NONE"
    hic[["hic"]][[dataName]][["obj"]] <-
      read_hic200(file, chrom_sizes = chrom_sizes, res = res_num)
  }
  hic
}

#' @keywords internal
.parse_resolution <- function(resolution) {
  if (is.numeric(resolution)) return(as.integer(resolution))
  r <- sub("kb", "000", resolution)
  as.integer(sub("bp", "", r))
}

#' @keywords internal
#' Fetch a region as a data.frame with columns x, y, counts (bp coordinates).
.get_region_counts <- function(hic, dataName, chr, pstart, pend, scale = 1) {
  d <- hic[["hic"]][[dataName]]
  res <- d[["resolution"]]
  if (identical(d[["format"]], "hic200")) {
    slice_hic200(d[["obj"]], chr, pstart, pend, scale = scale)
  } else {
    region <- paste(chr, max(pstart, 1), pend, sep = ":")
    df <- strawr::straw(d[["normalization"]], d[["file"]], region, region,
                        "BP", res)
    data.frame(x = df$x, y = df$y, counts = df$counts * scale)
  }
}
