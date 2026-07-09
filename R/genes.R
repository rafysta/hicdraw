#' Attach a gene model track
#'
#' @param hic A Hi-C object.
#' @param file Path to a BED-like gene file.
#' @param columnName Column names to assign (default matches a 6-column BED
#'   with `mix_id` = `gene_symbol`).
#' @return The updated Hi-C object.
#' @export
addGenedata <- function(hic, file,
                        columnName = c("chr", "start", "end", "mix_id",
                                       "dummy", "strand")) {
  df <- data.table::fread(file, header = FALSE, col.names = columnName,
                          select = seq_along(columnName))
  if ("mix_id" %in% colnames(df)) {
    df <- tidyr::separate(df, "mix_id", c("gene", "symbol"), "_",
                          convert = TRUE, extra = "merge")
  }
  hic[["gene"]] <- df
  hic
}

#' Draw the gene model track
#'
#' @param hic A Hi-C object.
#' @param label Draw gene symbols?
#' @param unit x-axis unit.
#' @param xaxis Show x axis text?
#' @return A ggplot object.
#' @export
drawGene <- function(hic, label = FALSE, unit = 1e6, xaxis = FALSE) {
  if (!requireNamespace("gggenes", quietly = TRUE)) {
    stop("package 'gggenes' is required for drawGene()")
  }
  pstart <- hic[["drawing"]][["start"]]
  pend <- hic[["drawing"]][["end"]]
  df <- dplyr::filter(hic[["gene"]], chr == hic[["drawing"]][["chr"]],
                      end > pstart, start < pend)

  p <- ggplot2::ggplot(df, ggplot2::aes(xmin = start / unit, xmax = end / unit,
      y = strand, fill = strand, x = (start + end) / 2 / unit,
      forward = ifelse(strand == "+", 1, -1))) +
    ggplot2::coord_cartesian(xlim = c(pstart / unit, pend / unit)) +
    ggplot2::scale_fill_manual(values = c("#fe9929", "#fee391")) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    gggenes::geom_gene_arrow(
      arrowhead_width = grid::unit(3, "mm"),
      arrowhead_height = grid::unit(5, "mm"), col = "black") +
    ggplot2::facet_wrap(~"", ncol = 1, strip.position = "right") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black")) +
    ggplot2::labs(y = "")

  if (!xaxis) {
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                            axis.text.x = ggplot2::element_blank())
  }
  if (label && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = symbol),
      seed = 42, box.padding = 0.5, size = 3)
  }
  p
}
