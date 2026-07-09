# Source data for pombe_chrom_sizes.
# The object is defined directly in R/data.R (as an exported data.frame) so no
# .rda build step is required. This file documents the provenance.
#
# Lengths from ASM294v2 (PomBase) nuclear chromosomes:
#   I   = 5,579,133 bp
#   II  = 4,539,804 bp
#   III = 2,452,883 bp
pombe_chrom_sizes <- data.frame(
  chr = c("I", "II", "III"),
  length = c(5579133, 4539804, 2452883),
  stringsAsFactors = FALSE
)
