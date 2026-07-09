#' Create an empty Hi-C drawing object
#'
#' The object is a plain list that accumulates data tracks (Hi-C maps,
#' bigwigs, boundary scores, gene models) and the current drawing region.
#'
#' @return An empty list used by the `add*` and `draw*` functions.
#' @export
createHiCmap <- function() {
  list()
}

#' @keywords internal
set_target <- function(hic, n = 1000) {
  pchr <- hic[["drawing"]][["chr"]]
  pstart <- hic[["drawing"]][["start"]]
  pend <- hic[["drawing"]][["end"]]

  binSize <- as.integer((pend - pstart) / n / 50) * 50
  if (binSize < 1) binSize <- 1L
  D_target <- data.frame(chr = pchr, start = seq(pstart, pend, by = binSize))
  D_target <- dplyr::mutate(D_target, end = start + binSize)
  D_target <- dplyr::filter(D_target, end <= pend)
  hic[["drawing"]][["df"]] <- D_target
  if (requireNamespace("GenomicRanges", quietly = TRUE)) {
    hic[["drawing"]][["grange"]] <- GenomicRanges::GRanges(D_target)
  }
  hic
}

#' Set the genomic region to draw
#'
#' @param hic A Hi-C object from [createHiCmap()].
#' @param chr Chromosome name.
#' @param start,end Region in bp.
#' @return The updated Hi-C object.
#' @export
setRegion <- function(hic, chr, start, end) {
  hic[["drawing"]][["chr"]] <- chr
  hic[["drawing"]][["start"]] <- start
  hic[["drawing"]][["end"]] <- end
  set_target(hic)
}

#' Print the current drawing region
#' @param hic A Hi-C object.
#' @export
showRegion <- function(hic) {
  cat(hic[["drawing"]][["chr"]], ":", hic[["drawing"]][["start"]], "-",
      hic[["drawing"]][["end"]], "\n")
  invisible(hic)
}
