# -----------------------------------------------------------------------------
# One shard of the coverage run, for .github/workflows/test-coverage.yaml.
#
# Coverage has to run the full NOT_CRAN suite serially: covr's line counters are
# per-process, so anything executed in a testthat parallel worker is never merged
# back (see the workflow header). Serial + covr-instrumented, the suite outgrew
# the job's 120-minute budget on PR #109 -- it was cancelled mid-run at ~116
# minutes of `test_check()` with nothing having failed.
#
# So the files are split across COV_SHARDS runners instead. Each shard installs
# and instruments the package itself, runs its own slice, and uploads its own
# cobertura report; Codecov merges uploads that carry the same commit SHA, so the
# per-PR diff and the badge are unchanged. Wall-clock becomes the slowest shard
# rather than the sum of all of them.
#
# Env:
#   COV_SHARD     1-based index of this shard           (required)
#   COV_SHARDS    how many shards the run is split into (required)
#   COV_PLAN_ONLY print the split and exit, running no tests -- for retuning the
#                 weights below, and the only mode that needs no covr
#   RUNNER_TEMP   where the instrumented package is installed (GitHub sets it)
#
# Run one shard locally with:
#   COV_SHARD=1 COV_SHARDS=4 Rscript tools/ci/coverage-shard.R
# -----------------------------------------------------------------------------

shard  <- as.integer(Sys.getenv("COV_SHARD",  "1"))
nshard <- as.integer(Sys.getenv("COV_SHARDS", "1"))
stopifnot(!is.na(shard), !is.na(nshard), nshard >= 1, shard >= 1, shard <= nshard)

test_dir <- normalizePath(file.path("tests", "testthat"), winslash = "/", mustWork = TRUE)
files    <- sort(list.files(test_dir, pattern = "^test-.*\\.[rR]$"))
stopifnot(length(files) > 0)


# --- Cost proxy --------------------------------------------------------------
# Balancing on file COUNT alone puts every heavy fit in one shard by accident, so
# each file is weighted by what it makes the model do. These weights are a static
# proxy read off the source, not measured seconds -- the ratios are what matter:
# an optimization dominates a build, and a diagnostic that refits in a loop
# dominates a single optimization.
#
# The workflow's reporter now prints a per-file elapsed time, so a future pass can
# replace this proxy with the real table.
COST_BASE  <- 1     # sourcing the file, fixtures, plain unit tests
COST_BUILD <- 1     # fit_mod(estimateMode = 3 or 4) -- MakeADFun, no optimization
COST_FIT   <- 12    # a real optimization (estimateMode 0/1/2, and the default)
COST_REFIT <- 40    # a .refit_like() caller: re-invokes fit_mod() in a loop

# The eight .refit_like() entry points (R/6-refit_like.R) plus sim_mod(), which
# draws through a built model.
REFIT_CALLS <- c("run_mse", "retrospective", "jitter", "self_test", "profile",
                 "remove_F", "sample_rec", "reweight_comps")

count_matches <- function(txt, pattern) {
  m <- gregexpr(pattern, txt, perl = TRUE)[[1]]
  if (identical(as.integer(m), -1L)) 0L else length(m)
}

file_cost <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")

  fits   <- count_matches(txt, "fit_mod\\s*\\(")
  builds <- count_matches(txt, "estimateMode\\s*=\\s*[34]\\b")
  refits <- sum(vapply(REFIT_CALLS,
                       function(f) count_matches(txt, paste0("\\b", f, "\\s*\\(")),
                       integer(1)))

  # Anything not explicitly built is assumed to optimize -- fit_mod()'s own
  # default estimateMode does.
  optimizations <- max(0L, fits - builds)

  COST_BASE +
    COST_BUILD * builds +
    COST_FIT   * optimizations +
    COST_REFIT * refits
}

costs <- vapply(file.path(test_dir, files), file_cost, numeric(1))
names(costs) <- files


# --- Longest-processing-time bin packing -------------------------------------
# Heaviest file first, each one onto whichever shard is currently lightest. For a
# spread like this suite's it lands within a few percent of a perfect split, and
# unlike round-robin it is stable: adding a light test file does not reshuffle
# the heavy ones onto different shards.
assign_shards <- function(costs, nshard) {
  ord    <- order(costs, decreasing = TRUE)
  load   <- numeric(nshard)
  bin    <- integer(length(costs))
  for (i in ord) {
    b        <- which.min(load)
    bin[i]   <- b
    load[b]  <- load[b] + costs[i]
  }
  list(bin = bin, load = load)
}

plan <- assign_shards(costs, nshard)
mine <- files[plan$bin == shard]
stopifnot(length(mine) > 0)


# --- Report the whole plan, so the log explains itself ------------------------
message(sprintf("Coverage shard %d of %d -- %d of %d test files",
                shard, nshard, length(mine), length(files)))
for (b in seq_len(nshard)) {
  message(sprintf("  shard %d: %3d files, weight %6.0f%s",
                  b, sum(plan$bin == b), plan$load[b],
                  if (b == shard) "   <- this runner" else ""))
}
message("\nFiles in this shard:")
for (f in mine) message(sprintf("  %-52s %5.0f", f, costs[[f]]))

if (nzchar(Sys.getenv("COV_PLAN_ONLY"))) {
  message("\nCOV_PLAN_ONLY set -- not running the tests.")
  quit(save = "no", status = 0)
}


# --- Run it under covr -------------------------------------------------------
# testthat's `filter` is matched against the file name with the `test-` prefix
# and the extension stripped, so anchor an alternation of exactly this shard's
# names. Every test file in the package is [a-z0-9_-] only, so none of them carry
# a regex metacharacter (test-schema-canonical.R would be the place to pin that
# if it ever changes).
context_names <- sub("\\.[rR]$", "", sub("^test-", "", mine))
filter <- paste0("^(", paste(context_names, collapse = "|"), ")$")

# Per-file timings, written as the run goes.
#
# covr::package_coverage(code = ) evaluates `code` in a subprocess and surfaces
# its output only when that returns, so on a shard killed by `timeout-minutes`
# the progress reporter's per-file times are never flushed -- the job log jumps
# from "* DONE (Rceattle)" straight to "The operation was canceled". The shard
# that overruns is the one that reports nothing, which is why the static weights
# above have never been checked against measured seconds.
#
# CoverageTimingReporter appends each file's time to this path as that file
# ends, so the table is on disk even when the job is cancelled; the workflow
# uploads it with `if: always()`. Feed the artifacts back into COST_* above.
timing_log <- Sys.getenv("COV_TIMING_LOG", unset = "")
if (!nzchar(timing_log)) {
  timing_log <- file.path(getwd(), sprintf("coverage-timing-shard-%d.tsv", shard))
}
timing_src <- normalizePath(file.path("tools", "ci", "coverage-timing.R"),
                            winslash = "/", mustWork = TRUE)
message(sprintf("\nPer-file timings -> %s", timing_log))

# `load_package = "installed"` picks up covr's instrumented build rather than the
# source tree, and `package =` makes the test environment's parent the Rceattle
# namespace so the internal helpers the suite calls resolve (see CLAUDE.md).
#
# stop_on_failure = FALSE: a failing test must not throw away the shard's
# coverage before it is written. R-CMD-check is the correctness gate for this
# suite -- test-coverage is deliberately non-blocking -- and the progress
# reporter prints any failure into this log.
#
# MultiReporter leaves the progress reporter's console output exactly as it was
# and writes the timing rows alongside it.
code <- sprintf(
  'source(%s)
   testthat::test_dir(%s, package = "Rceattle", load_package = "installed",
                      filter = %s, stop_on_failure = FALSE,
                      reporter = testthat::MultiReporter$new(reporters = list(
                        testthat::ProgressReporter$new(),
                        CoverageTimingReporter$new(%s))))',
  deparse1(timing_src), deparse1(test_dir), deparse1(filter),
  deparse1(timing_log)
)

install_path <- Sys.getenv("RUNNER_TEMP", unset = "")
install_path <- if (nzchar(install_path)) {
  file.path(normalizePath(install_path, winslash = "/"), "package")
} else {
  tempfile("Rceattle-covr-")
}
dir.create(install_path, recursive = TRUE, showWarnings = FALSE)

cov <- covr::package_coverage(
  quiet        = FALSE,
  clean        = FALSE,
  type         = "none",
  code         = code,
  install_path = install_path
)

covr::to_cobertura(cov)
message(sprintf("\nShard %d coverage: %s", shard, format(covr::percent_coverage(cov))))
