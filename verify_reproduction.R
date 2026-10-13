## =============================================================================
##  verify_reproduction.R
##
##  Compares the output of a new run of Screening_Simulation.R with the
##  reference output distributed in reference_output/, cell by cell.
##  Every CSV file present in both folders is compared; numeric cells must
##  agree within TOL, other cells exactly. A summary is printed and any
##  differing cells are written to Comparison_Report.csv.
##
##  Usage (from the repository root, after running Screening_Simulation.R):
##    Rscript verify_reproduction.R
##  Base R only.
## =============================================================================
MINE      <- "Screening_Output/tables"     # output of the new run
REFERENCE <- "reference_output/tables"     # reference output
TOL       <- 1e-6                          # results are seeded and should agree exactly

files <- intersect(list.files(MINE, "\\.csv$"), list.files(REFERENCE, "\\.csv$"))
only_mine <- setdiff(list.files(MINE, "\\.csv$"), files)
only_ref  <- setdiff(list.files(REFERENCE, "\\.csv$"), files)
cat("Files compared:", length(files), "\n")
if (length(only_ref))  cat("Missing from the new run:", paste(only_ref, collapse = ", "), "\n")
if (length(only_mine)) cat("Only in the new run:", paste(only_mine, collapse = ", "), "\n")

report <- list(); summary_rows <- list()
for (f in files) {
  a <- read.csv(file.path(MINE, f), check.names = FALSE, stringsAsFactors = FALSE)
  b <- read.csv(file.path(REFERENCE, f), check.names = FALSE, stringsAsFactors = FALSE)
  status <- "identical"; n_diff <- 0; max_diff <- 0
  if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) {
    status <- sprintf("STRUCTURE DIFFERS (new %dx%d, reference %dx%d)", nrow(a), ncol(a), nrow(b), ncol(b))
  } else {
    for (j in seq_along(a)) {
      x <- a[[j]]; y <- b[[j]]
      if (is.numeric(x) && is.numeric(y)) {
        d <- abs(x - y); bad <- which((is.na(x) != is.na(y)) | (!is.na(d) & d > TOL))
        if (length(bad)) {
          max_diff <- max(max_diff, d[bad], na.rm = TRUE); n_diff <- n_diff + length(bad)
          report[[length(report) + 1]] <- data.frame(file = f, row = bad, column = names(a)[j],
                                                     new = x[bad], reference = y[bad], abs_diff = d[bad])
        }
      } else {
        bad <- which(as.character(x) != as.character(y))
        if (length(bad)) {
          n_diff <- n_diff + length(bad)
          report[[length(report) + 1]] <- data.frame(file = f, row = bad, column = names(a)[j],
                                                     new = as.character(x[bad]), reference = as.character(y[bad]), abs_diff = NA)
        }
      }
    }
    if (n_diff) status <- "DIFFERENCES"
  }
  summary_rows[[f]] <- data.frame(file = f, status = status, cells_differing = n_diff, max_abs_diff = max_diff)
}
SUMMARY <- do.call(rbind, summary_rows); rownames(SUMMARY) <- NULL
options(width = 200); print(SUMMARY, right = FALSE)
if (length(report)) {
  REP <- do.call(rbind, report)
  write.csv(REP, "Comparison_Report.csv", row.names = FALSE)
  cat("\nDetails of", nrow(REP), "differing cells written to", normalizePath("Comparison_Report.csv"), "\n")
} else cat("\nAll compared files agree with the reference output.\n")

## Unit checks of the new run
chk <- file.path(MINE, "table_S1_unit_checks.csv")
if (file.exists(chk)) {
  u <- read.csv(chk); pc <- grep("pass", names(u), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(pc)) cat("\nUnit checks passed in the new run:", sum(as.logical(u[[pc]])), "of", nrow(u), "\n")
}
