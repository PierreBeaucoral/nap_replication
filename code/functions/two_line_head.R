# ==============================================================================
# code/functions/two_line_head.R
# Two-line LaTeX column heads for the wide results tables (03, 04, 05), so the
# tables fit the text width without being shrunk. Needs \usepackage{makecell}.
# ==============================================================================

#' Split a column label over two lines with \makecell
#'
#' Breaks at the word boundary closest to the middle of the label; a one-word
#' label is returned unchanged.
#'
#' @param x a single (already LaTeX-escaped) column label
#' @return the label wrapped as \makecell{first line \\ second line}
two_line_head <- function(x) {
  words <- strsplit(x, " ", fixed = TRUE)[[1L]]
  if (length(words) < 2L) return(x)
  ends <- cumsum(nchar(words) + 1L)[-length(words)]  # break after word k
  k    <- which.min(abs(ends - nchar(x) / 2))
  paste0("\\makecell{", paste(words[seq_len(k)], collapse = " "), " \\\\ ",
         paste(words[-seq_len(k)], collapse = " "), "}")
}
