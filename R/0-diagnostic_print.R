# One display contract for the diagnostics: every diagnostic prints the header
# convergence_diagnostics() uses -- an overall status, then a one-line verdict,
# then a compact table -- so Mohn's rho, OSA residuals and a jitter all answer
# "is this acceptable?" in the same four-level vocabulary.
#
# Return values are unchanged: each object is still the data frame or list it
# was. Only class() gains an entry and print() gets an opinion.

# Worst of a set of severities, using the convergence vocabulary so a reader
# meets the same four words everywhere. .CONV_SEVERITY / .conv_overall live in
# R/0-convergence.R.
.rce_worst <- function(sev) {
  sev <- sev[!is.na(sev)]
  if (!length(sev)) return("OK")
  .CONV_SEVERITY[max(match(sev, .CONV_SEVERITY))]
}

# The shared first line: `<Rceattle osa>  status: WARN`, then a one-line verdict.
.rce_diag_header <- function(what, status, verdict = NULL) {
  cat("<Rceattle ", what, ">  status: ", status, "\n", sep = "")
  if (!is.null(verdict) && nzchar(verdict)) cat("  ", verdict, "\n", sep = "")
  invisible(NULL)
}

# Severity tag in the same [WARN] column form print.Rceattle_convergence uses.
.rce_sev_tag <- function(sev) sprintf("[%-4s]", sev)

# Print a compact data frame without the row names or the wrapping. `cols` is a
# named character vector: names are the headers shown, values the columns taken.
.rce_diag_table <- function(df, cols, digits = 3, indent = "  ") {
  if (!nrow(df)) return(invisible(NULL))
  out <- lapply(names(cols), function(h) {
    v <- df[[cols[[h]]]]
    if (is.numeric(v)) formatC(v, format = "g", digits = digits) else as.character(v)
  })
  names(out) <- names(cols)
  w <- vapply(names(out), function(h) max(nchar(c(h, out[[h]])), na.rm = TRUE), 1L)
  cat(indent, paste(mapply(formatC, names(out), width = w, flag = "-"),
                    collapse = "  "), "\n", sep = "")
  for (i in seq_len(nrow(df))) {
    cat(indent, paste(mapply(function(h, k) formatC(out[[h]][i], width = w[[k]],
                                                    flag = "-"),
                             names(out), seq_along(out)),
                      collapse = "  "), "\n", sep = "")
  }
  invisible(NULL)
}

# "3 of 11" / "none of 11", for verdict lines.
.rce_n_of <- function(n, total, none = "none") {
  if (n == 0L) paste0(none, " of ", total) else paste0(n, " of ", total)
}
