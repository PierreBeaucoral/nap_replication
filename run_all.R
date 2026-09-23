# =============================================================================
# run_all.R  —  Master replication script
# Climate Finance and National Adaptation Plans (NAP)
# =============================================================================
# Reproduces every table and figure of the paper. See README.md for data
# sources, computational requirements and the exhibit-by-exhibit map.
#
# USAGE (from the project root)
#   Rscript run_all.R                     default: start from the shipped,
#                                         CRS-derived panels in data/processed/
#   NAP_FROM_RAW=1 Rscript run_all.R      full rebuild from the raw OECD CRS
#                                         files placed in data/raw/CRS/
#
# MODES
#   Default     Stage 01 (raw CRS -> data/processed/) is skipped: the analysis
#               starts at stage 02 from the processed panels shipped in
#               data/processed/. No raw CRS files are needed.
#   From raw    NAP_FROM_RAW=1
#               (a) verifies the SHA-256 checksum of every raw CRS file that
#                   stages 01 and 09 read (2007-2024) against
#                   data/raw/CRS/crs_checksums.csv and warns loudly when a file
#                   differs from the April 2026 vintage used in the paper;
#               (b) moves the shipped panels to data/processed_shipped/ and
#                   rebuilds data/processed/ from scratch with stage 01 (and,
#                   for the activity-level extract, stage 09);
#               (c) compares every rebuilt processed file with the shipped one
#                   and prints PASS/FAIL per file; the analysis then continues
#                   on the rebuilt data.
#
# PIPELINE (each stage runs in a fresh R process, in this order)
#   01_prepare_data.R              raw CRS/NAP/WDI/WGI -> data/processed/,
#                                  plus the three scope tables (from-raw mode only)
#   02_descriptive_stats.R         descriptive tables and figures
#   03_main_results.R              main Callaway-Sant'Anna results; stores the
#                                  headline fits in output/fits/ for later stages
#   04_robustness.R                robustness, HonestDiD, placebo, mitigation,
#                                  dCDH, Goodman-Bacon
#   05_heterogeneity.R             heterogeneity (donor type, LDC, governance,
#                                  income) and difference tests
#   13_cohort_anticipation.R       cohort battery, anticipation, event date,
#                                  placebo ladder, conditioning sets
#   14_hazard_napa.R               EM-DAT hazard controls and NAP-timing checks,
#                                  prior-NAPA split, NAPA falsification check
#   07_principal_and_share.R       principal-marker-only outcome, within-country
#                                  share, global-share reconciliation
#   08_randomization_inference.R   Fisher randomization inference
#   09_remarking_decomposition.R   activity-level re-marking decomposition
#   10_model_tests.R               testable implications of the model
#   11_base_year_sensitivity.R     2021-cohort base-year disclosure, full-window
#                                  pre-trend tests
#   12_group_figures.R             the two descriptive group figures
#
# ENVIRONMENT VARIABLES
#   NAP_FROM_RAW=1                 full rebuild from raw CRS (see MODES).
#   NAP_START_AT=<script file>     resume from that stage (earlier stages are
#                                  assumed complete), e.g.
#                                  NAP_START_AT=07_principal_and_share.R
#   NAP_SKIP_RI=1                  skip stage 08 when its cached permutation
#                                  draws (output/tables/randomization/ri_draws.csv)
#                                  exist. Without this flag stage 08 still uses
#                                  the cached draws and only rebuilds its table
#                                  and figure; delete the CSV to force the
#                                  full ~29-minute recompute.
#
# OUTPUT
#   data/processed/          analysis panels
#   output/tables/           LaTeX tables (.tex)
#   output/figures/          figures (.png, .pdf)
#   paper/Tables, paper/Figures   assembled from output/ only when a paper/
#                            folder exists next to this script (not part of
#                            the replication package).
# =============================================================================

library(here)

t0 <- Sys.time()

scripts <- c("01_prepare_data.R",
             "02_descriptive_stats.R",
             "03_main_results.R",
             "04_robustness.R",
             "05_heterogeneity.R",
             "13_cohort_anticipation.R",
             "14_hazard_napa.R",
             "07_principal_and_share.R",
             "08_randomization_inference.R",
             "09_remarking_decomposition.R",
             "10_model_tests.R",
             "11_base_year_sensitivity.R",
             "12_group_figures.R")

from_raw       <- identical(Sys.getenv("NAP_FROM_RAW"), "1")
skip_ri        <- identical(Sys.getenv("NAP_SKIP_RI"), "1")
ri_draws_cache <- here("output", "tables", "randomization", "ri_draws.csv")
processed_dir  <- here("data", "processed")
shipped_dir    <- here("data", "processed_shipped")

start_at <- Sys.getenv("NAP_START_AT", "")
if (nzchar(start_at) && !start_at %in% scripts) {
  stop("NAP_START_AT is not a pipeline stage: ", start_at)
}
if (!nzchar(start_at) && !from_raw) start_at <- "02_descriptive_stats.R"

# ---- Helpers for the from-raw mode ------------------------------------------

#' Verify the raw CRS files against the shipped SHA-256 table
#'
#' Stops if a required file is missing; warns (does not stop) if a file's
#' checksum differs, because a later OECD download can legitimately differ
#' from the April 2026 vintage behind the published results.
verify_crs_checksums <- function() {
  ref <- utils::read.csv(here("data", "raw", "CRS", "crs_checksums.csv"),
                         colClasses = c("integer", "character", "numeric", "character"))
  status <- character(nrow(ref))
  for (i in seq_len(nrow(ref))) {
    y  <- ref$year[i]
    fp <- here("data", "raw", "CRS", c(paste0("CRS ", y, " Data.txt"), paste0("CRS ", y, " data.txt")))
    fp <- fp[file.exists(fp)]
    if (length(fp) == 0L) {
      stop("Raw CRS file for ", y, " not found in data/raw/CRS/ (expected 'CRS ", y,
           " data.txt'). See data/raw/CRS/README.md.")
    }
    sha <- digest::digest(file = fp[1L], algo = "sha256")
    status[i] <- if (identical(sha, ref$sha256[i])) "MATCH" else "DIFFERENT VINTAGE"
    message(sprintf("  CRS %d  %-17s %s", y, status[i], sha))
  }
  n_bad <- sum(status != "MATCH")
  if (n_bad > 0L) {
    warning(sprintf(paste0(
      "%d of %d raw CRS files do NOT match the checksums of the April 2026 OECD vintage ",
      "used in the paper. The OECD revises past years, so results rebuilt from these files ",
      "can differ from the published ones; the panel comparison below shows where."),
      n_bad, nrow(ref)), call. = FALSE, immediate. = TRUE)
  } else {
    message("  All ", nrow(ref), " raw CRS files match the published vintage.")
  }
  invisible(status)
}

#' Compare one rebuilt data file with its shipped counterpart
#'
#' PASS if the files are byte-identical, or if their contents agree after
#' sorting rows on all columns (numeric tolerance 1e-10); FAIL otherwise.
compare_data_file <- function(rebuilt, shipped, reader) {
  if (unname(tools::md5sum(rebuilt)) == unname(tools::md5sum(shipped))) {
    return("PASS (byte-identical)")
  }
  a <- reader(rebuilt)
  b <- reader(shipped)
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) {
    return(sprintf("FAIL (columns or row count differ: %d vs %d rows)", nrow(a), nrow(b)))
  }
  data.table::setorderv(a, names(a), na.last = TRUE)
  data.table::setorderv(b, names(b), na.last = TRUE)
  eq <- all.equal(as.data.frame(a), as.data.frame(b), tolerance = 1e-10,
                  check.attributes = FALSE)
  if (isTRUE(eq)) "PASS (equal within 1e-10 after sorting)" else paste("FAIL:", eq[1L])
}

read_csv_dt <- function(f) data.table::fread(f, na.strings = c("NA", ""), encoding = "UTF-8")
read_gz_dt  <- function(f) {
  con <- gzfile(f, encoding = "UTF-8")
  on.exit(close(con))
  data.table::fread(text = readLines(con, warn = FALSE), na.strings = "NA",
                    colClasses = "character", encoding = "UTF-8")
}

#' Compare every shipped data file that has been rebuilt and not yet compared
#'
#' @param done character vector of files already compared
#' @return named character vector: result per newly compared file
compare_rebuilt <- function(done) {
  shipped <- list.files(shipped_dir, recursive = TRUE, pattern = "\\.csv(\\.gz)?$")
  todo    <- setdiff(shipped, done)
  todo    <- todo[file.exists(file.path(processed_dir, todo))]
  res     <- setNames(character(length(todo)), todo)
  for (f in todo) {
    reader <- if (grepl("\\.gz$", f)) read_gz_dt else read_csv_dt
    res[f] <- compare_data_file(file.path(processed_dir, f), file.path(shipped_dir, f), reader)
    message(sprintf("  %-62s %s", f, res[f]))
  }
  res
}
comparison <- character(0)  # file -> PASS/FAIL result, filled in from-raw mode

# ---- From-raw mode: checksums, then set the shipped panels aside -------------
if (from_raw) {
  message("\n=== NAP_FROM_RAW=1: full rebuild from raw OECD CRS files ===")
  message("Verifying SHA-256 checksums of the raw CRS files (2007-2024) ...")
  verify_crs_checksums()

  if (dir.exists(shipped_dir)) {
    # A previous from-raw run already set the shipped panels aside; what is in
    # data/processed/ now is that run's rebuild, so it is rebuilt again.
    message("Reference panels: data/processed_shipped/ (kept from a previous run).")
    unlink(processed_dir, recursive = TRUE)
  } else if (dir.exists(processed_dir)) {
    if (!file.rename(processed_dir, shipped_dir)) stop("Could not move data/processed/ aside.")
    message("Shipped panels moved to data/processed_shipped/ for comparison.")
  }
  dir.create(processed_dir, showWarnings = FALSE)
}

# ---- Run the pipeline stages in order, each in a fresh R process -------------
# Rscript (not source()) keeps stages independent: no objects leak between
# stages, each script's own seeds apply, and a failing stage stops the run.
n_steps <- length(scripts) + 1L
for (i in seq_along(scripts)) {
  script <- scripts[i]
  if (nzchar(start_at) && i < match(start_at, scripts)) {
    message(sprintf("\n[%d/%d] skipping %s (starting at %s; %s)", i, n_steps, script, start_at,
                    "earlier outputs, including the shipped data/processed/ panels, are used as is"))
    next
  }

  if (identical(script, "08_randomization_inference.R") &&
      skip_ri && file.exists(ri_draws_cache)) {
    message(sprintf(
      "\n[%d/%d] skipping %s (NAP_SKIP_RI=1 and cached draws found at %s)",
      i, n_steps, script, ri_draws_cache))
    next
  }

  message(sprintf("\n[%d/%d] running %s ...", i, n_steps, script))
  t_stage <- Sys.time()
  status <- system2("Rscript", shQuote(file.path(here("code"), script)))
  if (status != 0L) {
    stop(sprintf("Stage failed (exit status %d): %s\n", status, script),
         "Fix the error above, then re-run run_all.R.")
  }
  message(sprintf("[%d/%d] %s finished in %.1f min", i, n_steps, script,
                  as.numeric(difftime(Sys.time(), t_stage, units = "mins"))))

  if (from_raw && dir.exists(shipped_dir)) {
    comparison <- c(comparison, compare_rebuilt(names(comparison)))
  }
}

if (from_raw && dir.exists(shipped_dir)) {
  shipped_all <- list.files(shipped_dir, recursive = TRUE, pattern = "\\.csv(\\.gz)?$")
  not_rebuilt <- setdiff(shipped_all, names(comparison))
  failures    <- names(comparison)[startsWith(comparison, "FAIL")]
  message("\n=== Rebuilt vs. shipped data/processed/ ===")
  message(sprintf("  %d file(s) compared: %d PASS, %d FAIL; %d not rebuilt by this run%s",
                  length(comparison), length(comparison) - length(failures), length(failures),
                  length(not_rebuilt),
                  if (length(not_rebuilt)) paste0(": ", paste(not_rebuilt, collapse = ", ")) else ""))
  if (length(failures) > 0L) {
    warning("Rebuilt data differ from the shipped panels: ",
            paste(failures, collapse = ", "), call. = FALSE)
  }
}

# ---- Final step: normalise table floats; assemble output/ into paper/ --------
# Make generated result tables fit the page and stay where referenced:
# wrap each tabular in \adjustbox{max width=\textwidth} (shrinks only if too
# wide) and set the float to [H]. Skips longtables. Idempotent.
fix_result_tables <- function(dir) {
  files <- list.files(dir, pattern = "[.]tex$", recursive = TRUE, full.names = TRUE)
  for (f in files) {
    txt <- readLines(f, warn = FALSE)
    if (any(grepl("begin{longtable}", txt, fixed = TRUE))) next
    txt <- sub("\\\\begin\\{table\\}\\[[A-Za-z!]+\\]", "\\\\begin{table}[H]", txt)
    if (!any(grepl("adjustbox", txt, fixed = TRUE))) {
      out <- character(0)
      for (ln in txt) {
        if (grepl("^\\\\begin\\{tabular\\}", ln)) {
          out <- c(out, "\\adjustbox{max width=\\textwidth}{%", ln)
        } else if (grepl("^\\\\end\\{tabular\\}", ln)) {
          out <- c(out, ln, "}")
        } else out <- c(out, ln)
      }
      txt <- out
    }
    writeLines(txt, f)
  }
}
fix_result_tables(here("output", "tables"))

copy_tree <- function(src, dst) {
  # Only exhibit file types are assembled into paper/: caches (.csv) written
  # next to the tables stay in output/.
  files <- list.files(src, recursive = TRUE, full.names = FALSE,
                      pattern = "\\.(tex|pdf|png)$")  # excludes dotfiles
  for (f in files) {
    to <- file.path(dst, f)
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(src, f), to, overwrite = TRUE)
  }
  length(files)
}

if (dir.exists(here("paper"))) {
  message(sprintf("\n[%d/%d] assembling output/ -> paper/ ...", n_steps, n_steps))
  n_fig <- copy_tree(here("output", "figures"), here("paper", "Figures"))
  n_tab <- copy_tree(here("output", "tables"),  here("paper", "Tables"))
  message(sprintf("  assembled %d figures and %d tables into paper/.", n_fig, n_tab))
} else {
  message(sprintf("\n[%d/%d] no paper/ folder: exhibits are in output/tables and output/figures.",
                  n_steps, n_steps))
}

message(sprintf("\nDone in %.1f min.", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
