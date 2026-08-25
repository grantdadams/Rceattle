# One display contract for the diagnostics.
#
# convergence_diagnostics() already reads well: a classed object, a print method,
# a four-level severity, and an overall status on the first line. Nothing else in
# the family followed it. osa_diagnostics() returned sixteen unclassed columns at
# seven significant figures, which wraps across three screen-widths in an
# 80-column terminal -- so answering "is fleet 2's index residual acceptable?"
# meant reassembling one row from three blocks -- and said nothing about how many
# sources had failed. retrospective() returned Mohn's rho as a bare number with
# no reference band, while osa_diagnostics() shipped its null intervals: two
# diagnostics, opposite conventions on the same question. jitter() returned the
# objective values without the fraction-at-the-optimum that is the whole point of
# running it.
#
# These helpers give all of them the header convergence_diagnostics() uses. The
# RETURN VALUES are unchanged -- each object still is the data frame or list it
# was, and $mohns, $nll, $Rceattle_list, $species and every column still work.
# Only class() gains an entry and print() gets an opinion.

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
