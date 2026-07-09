#' Chromosome sizes for Schizosaccharomyces pombe
#'
#' Lengths (bp) of the three nuclear chromosomes, in the order used by the
#' hic200-cpp matrices in this project.
#'
#' @format data.frame with columns `chr` and `length`.
#' @export
pombe_chrom_sizes <- data.frame(
  chr = c("I", "II", "III"),
  length = c(5579133, 4539804, 2452883),
  stringsAsFactors = FALSE
)
