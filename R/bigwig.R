#' Attach a signal track (bigWig or data.frame)
#'
#' @param hic A Hi-C object.
#' @param data Either a path to a bigWig file (`type = "bw"`) or a data.frame
#'   with columns `chr, start, end, score` (and optionally `boundary`) when
#'   `type = "df"`.
#' @param dataName Name of the track.
#' @param type `"bw"` or `"df"`.
#' @return The updated Hi-C object.
#' @export
addBW <- function(hic, data, dataName = "data1", type = "bw") {
  if (type == "bw") {
    if (!file.exists(data)) {
      message(data, " is not found")
      return(hic)
    }
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      stop("package 'rtracklayer' is required for type = 'bw'")
    }
    hic[["bw"]][[dataName]][["data"]] <- rtracklayer::BigWigFile(data)
  } else if (type == "df") {
    hic[["bw"]][[dataName]][["data"]] <- data
  }
  hic[["bw"]][[dataName]][["type"]] <- type
  hic
}

#' Draw a signal track
#'
#' @param hic A Hi-C object.
#' @param dataName Name of an attached track.
#' @param type One of `"area"`, `"BorderStrength"`, `"tad"`, `"line"`.
#' @param fill,fill2 Fill colours for positive/negative values.
#' @param col Line colour (for `type = "line"`).
#' @param xaxis,xtitle Show x axis text / title.
#' @param ymin,ymax Y-axis limits.
#' @param ylabel Strip label (defaults to `dataName`).
#' @param unit x-axis unit.
#' @param ntick Approximate number of y ticks.
#' @return A ggplot object.
#' @export
drawBW <- function(hic, dataName = "data1", type = "area", fill = "#354863",
                   fill2 = NULL, col = "#354863", xaxis = TRUE, xtitle = TRUE,
                   ymin = NULL, ymax = NULL, ylabel = NULL, unit = 1e6,
                   ntick = 3) {
  pchr <- hic[["drawing"]][["chr"]]
  pstart <- hic[["drawing"]][["start"]]
  pend <- hic[["drawing"]][["end"]]

  tt <- hic[["bw"]][[dataName]][["type"]]
  if (tt == "bw") {
    bav <- rtracklayer::summary(hic[["bw"]][[dataName]][["data"]],
                                hic[["drawing"]][["grange"]], type = "max")
    df <- data.frame(dplyr::select(hic[["drawing"]][["df"]], start, end),
                     score = bav@unlistData$score)
  } else if (tt == "df") {
    df <- dplyr::filter(hic[["bw"]][[dataName]][["data"]],
                        chr == pchr, start >= pstart, end <= pend)
  } else {
    stop("data type not found")
  }

  if (is.null(fill2)) fill2 <- fill
  if (is.null(ylabel)) ylabel <- dataName

  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8))

  if (type %in% c("area", "BorderStrength")) {
    df <- dplyr::mutate(df, scoreCate = ifelse(score > 0, "p", "n"))
    p <- ggplot2::ggplot(df,
        ggplot2::aes(xmin = start / unit, xmax = end / unit, ymin = 0, ymax = score)) +
      ggplot2::geom_rect(alpha = 0.7, ggplot2::aes(fill = scoreCate)) +
      ggplot2::coord_cartesian(xlim = c(pstart / unit, pend / unit)) +
      ggplot2::scale_fill_manual(values = c(fill, fill2), breaks = c("p", "n")) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::facet_wrap(~"", strip.position = "right") +
      ggplot2::labs(x = paste0("Position in ", pchr), y = ylabel) + base_theme
    if (type == "BorderStrength" && "boundary" %in% colnames(df)) {
      bx <- dplyr::pull(dplyr::mutate(dplyr::filter(df, boundary == 1),
                                      b = (start + end) / 2 / unit), b)
      p <- p + ggplot2::geom_vline(xintercept = bx, col = "#354863", lty = 2)
    }
  } else if (type == "line") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = start / unit, y = score)) +
      ggplot2::geom_step(alpha = 0.7, col = col) +
      ggplot2::coord_cartesian(xlim = c(pstart / unit, pend / unit)) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::facet_wrap(~"", strip.position = "right") +
      ggplot2::labs(x = paste0("Position in ", pchr), y = ylabel) + base_theme
  } else if (type == "tad") {
    df <- dplyr::mutate(df, scoreCate = ifelse(dplyr::row_number() %% 2 == 0, "p", "n"))
    p <- ggplot2::ggplot(df,
        ggplot2::aes(xmin = start / unit, xmax = end / unit, ymin = -1, ymax = 1)) +
      ggplot2::geom_rect(alpha = 0.7, ggplot2::aes(fill = scoreCate)) +
      ggplot2::coord_cartesian(xlim = c(pstart / unit, pend / unit)) +
      ggplot2::scale_fill_manual(values = c(fill, fill2), breaks = c("p", "n")) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::facet_wrap(~"", strip.position = "right", drop = FALSE) +
      ggplot2::labs(x = paste0("Position in ", pchr), y = ylabel) + base_theme +
      ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank())
    ymin <- -1; ymax <- 1
  } else {
    stop("unknown draw type: ", type)
  }

  if (type != "tad") {
    if (is.null(ymax)) ymax <- max(df$score, na.rm = TRUE) * 1.05
    if (is.null(ymin)) ymin <- min(df$score, na.rm = TRUE)
  }
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0),
      limits = c(ymin, ymax), breaks = scales::pretty_breaks(n = ntick))

  if (!xaxis) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  if (!xtitle) p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  p
}
