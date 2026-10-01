# ==============================================================================
# code/functions/crs_positive.R
# The single "is this summed CRS amount positive?" test used where a stage
# classifies recipient-years as zero or positive (stages 03, 04, 07 and 11).
# ==============================================================================

#' Positive CRS sum, beyond floating-point residue
#'
#' CRS activity records can carry negative amounts (see log1p_crs() in
#' 01_prepare_data.R), so a recipient-year sum of offsetting records can be a
#' residue such as 1e-15 instead of an exact 0. A sum counts as positive only
#' above `tol` (USD millions) and as zero otherwise; a sum below -tol stops the
#' stage, as 01 does for the summed columns it logs.
#'
#' The check below stops the stage if the tolerance classifies any sum
#' differently from the exact tests (x > 0, x == 0), so that a CRS vintage for
#' which it matters is noticed rather than silently shifting the counts.
#'
#' @param x numeric vector of summed CRS amounts (USD millions); NA allowed
#' @param tol amounts in [-tol, tol] count as zero
#' @return logical vector, same length as x (NA where x is NA)
crs_positive <- function(x, tol = 1e-9) {
  if (any(x < -tol, na.rm = TRUE)) stop("Negative CRS sum: ", deparse(substitute(x)))
  pos <- x > tol
  if (!identical(pos, x > 0) || !identical(!pos, x == 0))
    stop("The zero tolerance changes the zero/positive split of ", deparse(substitute(x)),
         " (values in [-", tol, ", ", tol, "] other than 0); review the counts that depend on it.")
  pos
}
