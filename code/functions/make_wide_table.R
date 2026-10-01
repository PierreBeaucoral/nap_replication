# ==============================================================================
# code/functions/make_wide_table.R
# The single builder of the wide "simple ATT across outcomes" tables, sourced by
# stages 03 (Table 2), 04 (retained cohorts, not-yet-treated, India, balanced
# panel, 2010 panel) and 05 (donor type).
#
# Uses objects every calling stage defines at its top: BITERS,
# compute_pretrend_test() (functions/pretrend_test.R), two_line_head()
# (functions/two_line_head.R), esc_header(), wpval_reconciliation(),
# write_tex_float(), and the did, dplyr and xtable packages.
# ==============================================================================

#' Estimate one CS (2021) specification for several outcomes and write the wide
#' table of their simple ATTs
#'
#' Main specification (retain_thin = FALSE): doubly robust, multiplier-bootstrap
#' SEs (BITERS replications, set.seed(1242) immediately before each bootstrap
#' fit); the pre-trend test uses an analytical twin (bstrap = FALSE), whose
#' dynamic aggregation reports the influence-function SEs that the test's
#' covariance reproduces. Retained-cohorts specification (retain_thin = TRUE):
#' outcome regression with analytical SEs, one fit for the ATT and the test.
#' Any estimation failure stops the stage: no column is written as "---".
#'
#' @param did_panel_in estimation panel (thin cohorts already filtered as needed)
#' @param retain_thin FALSE for the main specification, TRUE for the
#'   retained-cohorts specification (see above)
#' @param outcomes list of list(var, label, ...), one per column
#' @param dir_tabs output folder of the table
#' @param tex_label LaTeX label of the table
#' @param caption_spec specification text appended to the caption
#' @param control_grp "nevertreated" or "notyettreated"
#' @param out_filename file name of the table in dir_tabs
#' @param after_column optional function(oc, col) called once per outcome after
#'   every reported quantity of that column has been computed (stage 03 stores
#'   the headline fits with it). `col` holds the fits and aggregations of the
#'   column: gt_boot, gt_analytical, agg_simple, agg_simple_analytic,
#'   agg_dyn_analytic, att, se, pretrend, n_obs, n_country, est_method, panel.
#' @return invisibly, the per-outcome statistics (formatted cells, unrounded
#'   ATT/SE and the analytical influence function used by stage 05's contrasts)
make_wide_table <- function(did_panel_in, retain_thin, outcomes,
                            dir_tabs, tex_label, caption_spec,
                            control_grp  = "nevertreated",
                            out_filename = "att_combined_wide.tex",
                            after_column = NULL) {

  stopifnot(control_grp %in% c("nevertreated", "notyettreated"),
            is.null(after_column) || is.function(after_column))

  # Log outcome -> raw column in levels (USD millions; pp for the share), for the
  # pre-treatment mean and the implied effect.
  raw_var_map <- c(
    log_commits           = "commitments",
    share_adapt           = "share_adapt",
    lcommitments_all      = "commitments_all",
    lcommitments_nonadapt = "commitments_nonadapt",
    ldisbursements        = "disbursements",
    log_commits_dac       = "commits_dac",
    log_commits_multi     = "commits_multi",
    log_commits_other     = "commits_other"
  )

  use_dr <- if (retain_thin) "reg" else "dr"
  control_grp_note <- if (control_grp == "notyettreated")
    "not-yet-treated control group" else "never-treated control group"

  fit_outcome <- function(yname, bstrap) {
    tryCatch(
      att_gt(
        yname         = yname,
        tname         = "year",
        idname        = "country_id",
        gname         = "cohort_year",
        xformla       = ~ ge_est + log_population,
        data          = did_panel_in,
        est_method    = use_dr,
        bstrap        = bstrap,
        biters        = BITERS,
        cband         = FALSE,
        control_group = control_grp,
        anticipation  = 0,
        base_period   = "universal",
        panel         = TRUE,
        allow_unbalanced_panel = TRUE
      ),
      error = function(e) stop("att_gt() failed for ", yname, " [", caption_spec, "]: ",
                               conditionMessage(e), call. = FALSE)
    )
  }

  table_stats <- setNames(vector("list", length(outcomes)),
                          vapply(outcomes, `[[`, character(1L), "var"))
  # Seed rule (reproduces the published SEs): seed before the outcomes loop and
  # immediately before each bootstrap fit.
  set.seed(1242)

  for (oc in outcomes) {
    message(sprintf("  Table stats [%s]: %s", caption_spec, oc$label))

    # Analytical fit first (consumes no random numbers).
    gt_tab_analytical <- fit_outcome(oc$var, bstrap = FALSE)
    if (!retain_thin) {
      set.seed(1242)
      gt_tab <- fit_outcome(oc$var, bstrap = TRUE)
    } else {
      gt_tab <- gt_tab_analytical  # same fit for the ATT and the pre-trend test
    }

    # Aggregations of the analytical fit consume no random numbers; aggte() on
    # the bootstrap fit draws its own multiplier bootstrap.
    agg_dyn_analytical <- aggte(gt_tab_analytical, type = "dynamic", na.rm = TRUE,
                                min_e = -5, max_e = Inf)
    agg_s            <- aggte(gt_tab, type = "simple", na.rm = TRUE)
    agg_s_analytical <- aggte(gt_tab_analytical, type = "simple", na.rm = TRUE)

    att <- agg_s$overall.att
    se  <- agg_s$overall.se
    t_v <- if (!is.na(att) && !is.na(se) && se > 0) att / se else NA_real_
    stars <- if (is.na(t_v)) "" else
      if (abs(t_v) > 2.576) "***" else
      if (abs(t_v) > 1.960) "**"  else
      if (abs(t_v) > 1.645) "*"   else ""

    pt <- compute_pretrend_test(agg_dyn_analytical, gt_tab_analytical)

    rows_oc   <- did_panel_in[!is.na(did_panel_in[[oc$var]]), ]
    n_obs     <- nrow(rows_oc)
    n_country <- length(unique(rows_oc$country_id))

    raw_var  <- raw_var_map[oc$var]
    is_share <- (oc$var == "share_adapt")

    if (!is.na(raw_var) && raw_var %in% names(did_panel_in)) {
      pre_rows <- did_panel_in %>%
        filter(cohort_year > 0, year < cohort_year,
               !is.na(.data[[raw_var]]))
      mean_pre <- mean(pre_rows[[raw_var]], na.rm = TRUE)
    } else {
      mean_pre <- NA_real_
    }

    # The naive back-transform exp(ATT)-1 x pre-treatment mean is not
    # additive across outcomes (it violates the paper's own additivity check),
    # so it is reported only when the underlying ATT is significant at 5%
    # (|t| > 1.960); otherwise the cell is suppressed ("---").
    is_sig5 <- !is.na(t_v) && abs(t_v) > 1.960
    if (!is_share && !is.na(att) && !is.na(mean_pre) && is_sig5) {
      implied_usd <- (exp(att) - 1) * mean_pre
    } else {
      implied_usd <- NA_real_
    }

    mean_pre_fmt <- if (is_share) {
      sprintf("%.2f pp", mean_pre)
    } else if (!is.na(mean_pre)) {
      sprintf("%.1f", mean_pre)
    } else {
      "---"
    }

    implied_fmt <- if (is_share) {
      "---"
    } else if (!is.na(implied_usd)) {
      sprintf("%.1f", implied_usd)
    } else {
      "---"
    }

    table_stats[[oc$var]] <- list(
      att_fmt      = paste0(sprintf("%.4f", att), stars),
      se_fmt       = sprintf("(%.4f)", se),
      t_fmt        = sprintf("%.3f", t_v),
      mean_pre_fmt = mean_pre_fmt,
      implied_fmt  = implied_fmt,
      n_obs        = format(n_obs, big.mark = ","),
      n_country    = as.character(n_country),
      pt_stat      = sprintf("%.3f", pt$stat),
      pt_pval      = sprintf("%.3f", pt$pval),
      pt_df        = pt$df,
      pt_df_did    = pt$df_did,
      pt_leads     = pt$leads,
      pt_ginv      = pt$ginv_used,
      pt_reason    = pt$wpval_reason,
      pt_pval_num  = pt$pval,
      pt_wpval_num = pt$Wpval_did,
      pt_label     = oc$label,
      pt_wpval_did = if (is.na(pt$Wpval_did)) "---" else sprintf("%.3f", pt$Wpval_did),
      # Unrounded values and the unit-level influence function of the
      # analytical twin (stage 05 contrasts donor-type ATTs estimated on the
      # same recipient-years with them); `ids` is the unit order of its rows.
      att_num      = att,
      se_num       = se,
      inf_func     = agg_s_analytical$inf.function$simple.att,
      se_analytic  = agg_s_analytical$overall.se,
      att_analytic = agg_s_analytical$overall.att,
      ids          = sort(unique(did_panel_in$country_id))
    )
    message(sprintf("    ATT = %s  |  Pre-trend chi2(%d) = %.3f  p = %.3f  |  did Wpval = %s",
                    table_stats[[oc$var]]$att_fmt, pt$df, pt$stat, pt$pval,
                    table_stats[[oc$var]]$pt_wpval_did))
    message(sprintf("    Mean pre-treat (raw) = %s  |  Implied effect (USD M) = %s",
                    mean_pre_fmt, implied_fmt))

    # Runs after every reported quantity of this column, so it cannot perturb
    # the draws behind them; the next column re-seeds before its bootstrap fit.
    if (!is.null(after_column)) {
      after_column(oc, list(
        gt_boot = gt_tab, gt_analytical = gt_tab_analytical,
        agg_simple = agg_s, agg_simple_analytic = agg_s_analytical,
        agg_dyn_analytic = agg_dyn_analytical, att = att, se = se, pretrend = pt,
        n_obs = n_obs, n_country = n_country, est_method = use_dr,
        panel = did_panel_in))
    }
  }

  row_labels <- c(
    "ATT",
    "SE",
    "$t$-statistic",
    "Pre-treatment mean",
    "Implied effect (USD M)",
    "Observations",
    "Countries",
    "Pre-trend $\\chi^2$",
    "Pre-trend $p$",
    "\\texttt{did} pre-test $p$"
  )

  col_short <- vapply(outcomes, `[[`, character(1L), "label")
  tab_wide  <- data.frame(` ` = row_labels, check.names = FALSE, stringsAsFactors = FALSE)

  for (i in seq_along(outcomes)) {
    oc  <- outcomes[[i]]
    s   <- table_stats[[oc$var]]
    stopifnot(!is.null(s))
    col <- c(
      s$att_fmt, s$se_fmt, s$t_fmt,
      s$mean_pre_fmt, s$implied_fmt,
      s$n_obs, s$n_country,
      s$pt_stat, s$pt_pval, s$pt_wpval_did
    )
    tab_wide[[two_line_head(esc_header(col_short[i]))]] <- col
  }

  # Estimator, control group and inference are stated in the note, and a
  # sample-size row that is identical in every column is stated there once,
  # so the table fits the text width without being shrunk.
  n_row_words <- c(Observations = "observations (recipient-years)", Countries = "countries")
  n_row_const <- Filter(function(r)
    length(unique(unlist(tab_wide[tab_wide[[1L]] == r, -1L]))) == 1L, names(n_row_words))
  const_note  <- if (length(n_row_const) == 0L) "" else
    paste0("Every column: ", paste(vapply(n_row_const, function(r)
      paste(tab_wide[tab_wide[[1L]] == r, 2L], n_row_words[[r]]), character(1L)),
      collapse = ", "), ". ")
  tab_wide <- tab_wide[!tab_wide[[1L]] %in% n_row_const, , drop = FALSE]

  xtab <- xtable(tab_wide, label = tex_label)
  align(xtab) <- paste0("ll", paste(rep("c", length(col_short)), collapse = ""))

  raw_lines <- capture.output(
    print(xtab, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small",
          floating = FALSE)
  )

  # Degrees of freedom behind the two pre-trend columns, taken from the first
  # outcome (they are identical across outcomes here, as all are estimated on
  # the same panel), and per-column inputs for the pre-trend note: the lead
  # window used, both p-values for every outcome, whether a generalized inverse
  # was needed, and did's own reason when it returned no statistic.
  first_pt_df_did <- table_stats[[1L]]$pt_df_did
  first_n_country <- as.integer(table_stats[[1L]]$n_country)
  first_pt_leads  <- table_stats[[1L]]$pt_leads
  col_wpval <- vapply(table_stats, function(s)
    if (is.null(s$pt_wpval_num)) NA_real_ else s$pt_wpval_num, numeric(1L))
  col_pwald <- vapply(table_stats, function(s)
    if (is.null(s$pt_pval_num)) NA_real_ else s$pt_pval_num, numeric(1L))
  col_labels_pt <- vapply(table_stats, function(s) s$pt_label, character(1L))
  col_ginv  <- vapply(table_stats, function(s) isTRUE(s$pt_ginv), logical(1L))
  col_df    <- vapply(table_stats, function(s) as.integer(s$pt_df), integer(1L))
  col_reason <- vapply(table_stats, function(s)
    if (is.null(s$pt_reason)) NA_character_ else s$pt_reason, character(1L))

  est_method_label <- if (use_dr == "dr") "DR" else "OR"
  se_label <- if (!retain_thin)
    paste0("multiplier-bootstrap SE (", BITERS,
           " reps, seed 1242; \\texttt{did} 2.5.0; CRS Apr.\\ 2026)") else
    "analytical (IF) SE; \\texttt{did} 2.5.0; CRS Apr.\\ 2026"

  notes_txt <- paste0(
    "CS\\,(2021) ", est_method_label,
    "; WGI gov.\\ effectiveness + log population; ",
    control_grp_note, "; ", se_label, ". ",
    wpval_reconciliation(pre_egt = first_pt_leads, df_did = first_pt_df_did,
                         n_clusters = first_n_country, wpval_did = col_wpval,
                         pval_wald = col_pwald, labels = col_labels_pt,
                         wpval_reason = col_reason, ginv_used = col_ginv,
                         df = col_df),
    "Pre-treatment mean: mean of the outcome in levels over treated recipient-years before ",
    "adoption, USD M", if ("share_adapt" %in% vapply(outcomes, `[[`, character(1L), "var"))
      " (pp for the share)" else "", ". Implied effect: $(\\exp(\\mathrm{ATT})-1) \\times$ ",
    "pre-treatment mean, a naive back-transform, not additive across outcomes; suppressed ",
    "if the ATT is insignificant at 5\\%. ", const_note,
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  )
  source_txt <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

  cap_title <- paste0(
    "Effect of NAP adoption on climate finance: simple ATT across outcomes (",
    caption_spec, ")"
  )

  write_tex_float(file.path(dir_tabs, out_filename), cap_title, tex_label,
                  raw_lines, notes_txt, source_txt)

  invisible(table_stats)
}
