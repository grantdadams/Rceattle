# -----------------------------------------------------------------------------
# A testthat reporter that appends one line per test file, as each file ends.
#
# Why this exists: tools/ci/coverage-shard.R runs the suite inside
# covr::package_coverage(code = ), which evaluates that code in a SUBPROCESS and
# surfaces its output only once it returns. The progress reporter already prints
# a per-file time, but on the run that matters -- the one killed by
# `timeout-minutes` -- none of it is ever flushed, so the job log jumps straight
# from "* DONE (Rceattle)" to "The operation was canceled" with no indication of
# which files ran or how long they took. That is how the static cost proxy in
# coverage-shard.R went unchecked against reality: the shard that overruns is
# exactly the shard that reports nothing.
#
# Writing to a file instead, reopened in append mode per line, puts each timing
# on disk the instant that file finishes, so it survives the job being
# cancelled. The workflow uploads it with `if: always()`.
#
# One tab-separated `file<TAB>seconds` row plus a header, so a later pass can
# read the artifacts straight into a data frame and replace the proxy weights
# with measured medians.
# -----------------------------------------------------------------------------

CoverageTimingReporter <- R6::R6Class(
  "CoverageTimingReporter",
  inherit = testthat::Reporter,
  public = list(
    path = NULL,
    file_name = NULL,
    started = NULL,

    initialize = function(path) {
      super$initialize()
      self$path <- path
      # Truncate and write the header once, so a re-run does not append to a
      # previous shard's table.
      cat("file\tseconds\n", file = self$path)
      invisible(self)
    },

    start_file = function(filename) {
      self$file_name <- filename
      self$started   <- Sys.time()
      invisible(self)
    },

    end_file = function() {
      if (is.null(self$started)) return(invisible(self))
      secs <- as.numeric(difftime(Sys.time(), self$started, units = "secs"))
      # append = TRUE reopens and closes the connection for each row, so it is
      # durable before the next file starts. ~150 rows a run, so the cost of
      # that is irrelevant next to what it buys on a cancelled shard.
      cat(sprintf("%s\t%.1f\n", self$file_name, secs),
          file = self$path, append = TRUE)
      self$file_name <- NULL
      self$started   <- NULL
      invisible(self)
    }
  )
)
