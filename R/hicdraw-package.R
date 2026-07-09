#' @keywords internal
"_PACKAGE"

#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @importFrom utils head
NULL

# Silence R CMD check notes about non-standard-evaluation column names used
# throughout the dplyr / ggplot2 pipelines.
utils::globalVariables(c(
  "bin1", "bin2", "score", "start", "end", "start1", "end1", "start2", "end2",
  "x", "y", "x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4", "g", "counts",
  "chr", "strand", "symbol", "mix_id", "scoreCate", "boundary", "cate", "ord",
  "BS.norm"
))

#' @export
magrittr::`%>%`
