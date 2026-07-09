#' Build a genome-wide bin index for hic200-cpp matrices
#'
#' hic200-cpp stores contacts as `bin1 bin2 score`, where the bin is a
#' genome-wide, 0-based index at a fixed resolution. Chromosomes are
#' concatenated in the order given by `chrom_sizes`, so the first bin of a
#' chromosome starts at the cumulative bin count of all preceding chromosomes.
#'
#' @param chrom_sizes data.frame with columns `chr` and `length` (bp), in the
#'   same order used when the matrix was generated.
#' @param res Bin size in bp (default 200).
#'
#' @return data.frame with columns `chr`, `length`, `nbin`, `bin_offset`.
#' @export
#' @examples
#' hic200_index(pombe_chrom_sizes, res = 200)
hic200_index <- function(chrom_sizes, res = 200L) {
  stopifnot(all(c("chr", "length") %in% colnames(chrom_sizes)))
  nbin <- as.integer(ceiling(chrom_sizes$length / res))
  offset <- cumsum(c(0L, nbin[-length(nbin)]))
  data.frame(
    chr = as.character(chrom_sizes$chr),
    length = as.numeric(chrom_sizes$length),
    nbin = nbin,
    bin_offset = as.integer(offset),
    stringsAsFactors = FALSE
  )
}

#' Convert a genomic region to genome-wide bin indices
#'
#' @param index Output of [hic200_index()].
#' @param chr Chromosome name.
#' @param start,end Region start/end in bp.
#' @param res Bin size in bp.
#' @return Integer vector `c(bin0, bin1)` (inclusive).
#' @export
bins_for_region <- function(index, chr, start, end, res = 200L) {
  row <- index[index$chr == as.character(chr), , drop = FALSE]
  if (nrow(row) == 0) stop("chromosome '", chr, "' not found in index")
  b0 <- row$bin_offset + as.integer(floor(start / res))
  b1 <- row$bin_offset + as.integer(floor(end / res))
  c(b0, b1)
}

#' Read a hic200-cpp sparse matrix file
#'
#' Reads the whole `bin1 bin2 score` table into memory and keeps the bin index
#' so that regions can be sliced quickly and repeatedly.
#'
#' @param file Path to a (optionally gzipped) hic200-cpp matrix file.
#' @param chrom_sizes data.frame with `chr` and `length`; defaults to
#'   [pombe_chrom_sizes].
#' @param res Bin size in bp (default 200).
#'
#' @return An object of class `hic200` (a list with `dt`, `index`, `res`,
#'   `file`).
#' @export
read_hic200 <- function(file, chrom_sizes = pombe_chrom_sizes, res = 200L) {
  if (!file.exists(file)) stop(file, " is not found")
  dt <- data.table::fread(file, header = TRUE,
                          col.names = c("bin1", "bin2", "score"))
  index <- hic200_index(chrom_sizes, res = res)
  structure(list(dt = dt, index = index, res = as.integer(res), file = file),
            class = "hic200")
}

#' Slice a genomic region out of a hic200 object
#'
#' Returns a symmetric (both triangles) long table of contacts within the
#' region, with genomic bp coordinates.
#'
#' @param obj An `hic200` object from [read_hic200()].
#' @param chr Chromosome name.
#' @param start,end Region in bp.
#' @param scale Numeric multiplier applied to `counts` (e.g. for depth
#'   normalisation). Default 1.
#'
#' @return data.frame with columns `x`, `y` (bp) and `counts`.
#' @export
slice_hic200 <- function(obj, chr, start, end, scale = 1) {
  stopifnot(inherits(obj, "hic200"))
  res <- obj$res
  row <- obj$index[obj$index$chr == as.character(chr), , drop = FALSE]
  if (nrow(row) == 0) stop("chromosome '", chr, "' not found in index")
  bb <- bins_for_region(obj$index, chr, start, end, res = res)
  b0 <- bb[1]; b1 <- bb[2]
  dt <- obj$dt
  sub <- dt[dt$bin1 >= b0 & dt$bin1 <= b1 & dt$bin2 >= b0 & dt$bin2 <= b1, ]
  off <- row$bin_offset
  x <- (sub$bin1 - off) * res
  y <- (sub$bin2 - off) * res
  counts <- sub$score * scale
  # symmetrise: add lower triangle for off-diagonal pairs
  offdiag <- x != y
  data.frame(
    x = c(x, y[offdiag]),
    y = c(y, x[offdiag]),
    counts = c(counts, counts[offdiag])
  )
}
