# ==============================================================================
# code/functions/make_country_id.R
# The single recipient-name -> integer unit ID used as att_gt()'s idname by
# stages 03, 04, 05, 07, 08, 09, 10, 11, 13 and 14.
# ==============================================================================

#' Locale-independent integer unit ID from recipient names
#'
#' Why not `as.integer(factor(name))`: factor() sorts its levels with the
#' session's collation, and names with non-ASCII letters sort differently under
#' the C locale ("Croatia" < "Cuba" < "Côte d'Ivoire", "Turkmenistan" <
#' "Tuvalu" < "Türkiye") and under language locales such as en_US or fr_FR
#' ("Côte d'Ivoire" < "Croatia", "Türkiye" < "Turkmenistan"). The IDs of those
#' recipients would then depend on the machine, and with them the unit labels
#' of every fit and the randomization-inference design hash.
#'
#' The published results were produced with the language-locale order: R on
#' macOS takes LC_COLLATE from the system locale (here fr_FR) when LANG and
#' LC_ALL are unset, even where the shell's `locale` reports "C". The ICU
#' collation for "en_US" (stringi, which carries its own collation rules)
#' gives that order for these names whatever the session locale.
#'
#' The ranks are computed on the distinct names and mapped back by exact
#' matching, so the result does not depend on how the session re-encodes
#' non-ASCII strings (stringi reads them as UTF-8, which the panels are).
#'
#' @param name character vector of recipient names (one entry per panel row)
#' @return integer vector, same length as name: the rank of each name in the
#'   collation-ordered set of distinct names
make_country_id <- function(name) {
  stopifnot(is.character(name), !anyNA(name))
  u    <- unique(name)
  rank <- stringi::stri_rank(u, locale = "en_US")
  if (anyDuplicated(rank)) stop("make_country_id(): distinct names collate as equal")
  rank[match(name, u)]
}
