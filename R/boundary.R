#' Attach a hic200-cpp boundary-score track
#'
#' Reads a `*_BS.txt` file (columns `chr start end BS BS.norm boundary TADid TAD`)
#' and stores it as a signal track usable with `drawBW(type = "BorderStrength")`.
#'
#' @param hic A Hi-C object.
#' @param file Path to a `*_BS.txt` file.
#' @param dataName Track name.
#' @param score_column Which column to use as the drawn score
#'   (default `"BS.norm"`).
#' @return The updated Hi-C object.
#' @export
addBoundary <- function(hic, file, dataName = "boundary",
                        score_column = "BS.norm") {
  if (!file.exists(file)) {
    message(file, " is not found")
    return(hic)
  }
  df <- data.table::fread(file)
  df <- as.data.frame(df)
  df$score <- suppressWarnings(as.numeric(df[[score_column]]))
  keep <- intersect(c("chr", "start", "end", "score", "boundary", "TAD"),
                    colnames(df))
  addBW(hic, df[, keep, drop = FALSE], dataName = dataName, type = "df")
}

#' Attach a hic200-cpp domain-call track
#'
#' Reads a `*_BS_domains.txt` file (columns
#' `chr start end domain_id n_bins is_TAD`).
#'
#' @param hic A Hi-C object.
#' @param file Path to a `*_BS_domains.txt` file.
#' @param dataName Track name.
#' @return The updated Hi-C object.
#' @export
addDomains <- function(hic, file, dataName = "domain") {
  if (!file.exists(file)) {
    message(file, " is not found")
    return(hic)
  }
  df <- as.data.frame(data.table::fread(file))
  addBW(hic, df, dataName = dataName, type = "df")
}
