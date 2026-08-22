# The column schema is the source of truth for workbook columns, but it is not
# the only registry in the package, and the others are hand-synced. Each test
# here binds one of them to the thing it must agree with, so a value added in
# one place and forgotten in another fails here rather than in a fit.

test_that("write_template() emits every fleet_control column the schema defines", {
  # A template missing a column is how a user never learns the column exists:
  # they fill in what the template gave them and the feature stays off.
  testthat::skip_on_cran()   # writes and reads a workbook

  f <- tempfile(fileext = ".xlsx")
  suppressMessages(write_template(f, nages = 5, nyrs = 6, nprojyrs = 2, minage = 1))
  d <- suppressMessages(read_data(f))
  want <- .rce_schema_names("fleet_control")
  testthat::expect_length(setdiff(want, names(d$fleet_control)), 0)

  # And the template must be a valid model on its own: a starter workbook that
  # does not pass the data check is worse than none.
  d2 <- suppressMessages(suppressWarnings(switch_check(d)))
  testthat::expect_length(validate_switches(d2), 0)
})

test_that("every model-level switch's allowed set survives into the YAML comment", {
  # .rce_config_schema() projects from .rce_model_switch_schema(), so comparing
  # the two would be a tautology. What is worth pinning is the PROJECTION: that
  # every value of the governing map reaches the comment text a user reads in a
  # saved config, exactly and not as a substring of a longer name.
  dict <- .rce_config_schema()
  pairs <- list(estimateMode = estimateMode_map, initMode = initMode_map,
                suitMode = suitMode_map, msmMode = msmMode_map)

  for (nm in names(pairs)) {
    entry <- dict[[nm]]
    testthat::expect_false(is.null(entry), info = paste(nm, "absent from the config schema"))
    if (is.null(entry)) next
    # The comment is written as "A / B / C"; split it rather than substring-match,
    # or `MSVPA` is satisfied by `TypeIIIMSVPA`.
    got <- trimws(strsplit(paste(entry$allowed, collapse = " "), "/", fixed = TRUE)[[1]])
    testthat::expect_setequal(got[nzchar(got)], names(pairs[[nm]]))
  }
})


test_that("the hand-coded defaults still equal the schema's", {
  # Some defaults cannot be schema-driven: the bioenergetics scalars are
  # per-species (`rep(x, nspp)`) and are filled only in SINGLE-species mode
  # (`if (msmMode == 0)`), precisely because they are never read there -- which
  # the schema's flat `default` field cannot express. They
  # stay imperative in switch_check(), and the schema rows carry the same
  # values as documentation. Nothing kept the two in step; this does.
  #
  # A mismatch means the meta sheet and the field dictionary are telling users
  # one default while switch_check() applies another.
  sc <- .rce_column_schema()

  # Read the values out of switch_check() rather than restating them here --
  # a literal copy in the test is a third place for them to drift, and would
  # stay green when the source changed, which is the one thing this must catch.
  sw <- c("R/0-switches.R", testthat::test_path("..", "..", "R", "0-switches.R"))
  sw <- sw[file.exists(sw)]
  testthat::skip_if(length(sw) == 0, "R/0-switches.R not available")
  lines <- readLines(sw[1], warn = FALSE)
  from <- grep("bioenergetics_defaults <- list\\(", lines)[1]
  testthat::expect_false(is.na(from))
  to <- from + which(grepl("^\\s*\\)\\s*$", lines[(from + 1):(from + 30)]))[1]
  block <- lines[(from + 1):(to - 1)]

  hand_coded <- list()
  for (l in block) {
    m <- regmatches(l, regexec("^\\s*([A-Za-z0-9_]+)\\s*=\\s*rep\\(\\s*([0-9.eE+-]+)L?", l))[[1]]
    if (length(m) == 3) hand_coded[[m[2]]] <- as.numeric(m[3])
  }
  testthat::expect_gt(length(hand_coded), 8L)

  checked <- 0L
  for (nm in names(hand_coded)) {
    row <- sc[[nm]]
    if (is.null(row) || is.null(row$default) || all(is.na(row$default))) next
    checked <- checked + 1L
    testthat::expect_equal(
      as.numeric(row$default), as.numeric(hand_coded[[nm]]),
      info = paste0(nm, ": switch_check() applies ", hand_coded[[nm]],
                    " but the schema documents ", row$default))
  }
  # If the schema stops carrying these, the test must fail rather than pass on
  # an empty loop.
  testthat::expect_gt(checked, 8L)
})


test_that("the model-level switch table is complete and internally consistent", {
  # The switches that configure the model rather than a workbook column. They
  # lived in maps, were defaulted in switch_check(), documented in fit_mod()'s
  # roxygen and re-described in .rce_config_schema() -- four places, hand-kept.
  t <- .rce_model_switch_schema()
  testthat::expect_gt(length(t), 10)

  for (nm in names(t)) {
    row <- t[[nm]]
    testthat::expect_identical(row$name, nm)
    testthat::expect_true(nzchar(row$doc), info = paste(nm, "has no doc"))
    testthat::expect_true(row$scope %in% c("scalar", "per-species", "per-fleet"),
                          info = paste(nm, "has an unknown scope:", row$scope))

    # A named map must resolve, or the config comment silently loses its values.
    if (!is.na(row$allowed)) {
      m <- tryCatch(get(row$allowed, envir = asNamespace("Rceattle")),
                    error = function(e) NULL)
      testthat::expect_false(is.null(m),
                             info = paste0(nm, ": allowed = '", row$allowed,
                                           "' does not resolve"))
      testthat::expect_gt(length(m), 1)
      # The default must be a value the map actually defines.
      if (!is.na(row$default)) {
        testthat::expect_true(as.numeric(row$default) %in% unname(m),
                              info = paste0(nm, ": default ", row$default,
                                            " is not in ", row$allowed))
      }
    }
  }
})

# (The former "projects, not a hand-kept copy" test was removed: once the
#  projection landed it compared a value with itself. The test above pins the
#  property that actually matters -- that the values reach the user.)


test_that("scope and tmb_target on the switch table are true, not just present", {
  # The table is read by anyone deciding whether a switch takes one value or a
  # vector. `srr_fun` is a DATA_INTEGER scalar; calling it per-species would
  # invite a length error or a silent first-element truncation.
  t <- .rce_model_switch_schema()
  cpp <- c("src/TMB/ceattle.cpp", testthat::test_path("..", "..", "src", "TMB", "ceattle.cpp"))
  cpp <- cpp[file.exists(cpp)]
  testthat::skip_if(length(cpp) == 0, "src/TMB/ceattle.cpp not available")
  src <- paste(readLines(cpp[1], warn = FALSE), collapse = "\n")

  for (nm in names(t)) {
    tgt <- t[[nm]]$tmb_target
    if (is.na(tgt)) next
    # A declared target must exist, and its DATA_ kind must match the scope:
    # DATA_INTEGER / DATA_SCALAR is one value, DATA_IVECTOR / DATA_VECTOR many.
    m <- regmatches(src, regexpr(paste0("DATA_[A-Z]+\\(\\s*", tgt, "\\s*\\)"),
                                 src, perl = TRUE))
    testthat::expect_length(m, 1)
    if (!length(m)) next
    is_vector <- grepl("VECTOR", m, fixed = TRUE)
    testthat::expect_equal(
      is_vector, t[[nm]]$scope != "scalar",
      info = paste0(nm, ": declared scope '", t[[nm]]$scope, "' but the template says ", m))
  }
})


test_that("a schema column declaring `allowed` is either validated or knowingly not", {
  # `allowed` says which values are legal. A column that declares one and is
  # never checked documents a rule nothing enforces -- the state this work set
  # out to end. These four are known-unvalidated; the list must not grow
  # silently.
  sw <- c("R/0-switches.R", testthat::test_path("..", "..", "R", "0-switches.R"))
  sw <- sw[file.exists(sw)]
  testthat::skip_if(length(sw) == 0, "R/0-switches.R not available")
  lines <- readLines(sw[1], warn = FALSE)
  start <- grep("^validate_switches <- function", lines)[1]
  # Bound the function at the next top-level definition, or the whole tail is
  # swept in and any column named by a LATER helper looks validated.
  after <- grep("^[A-Za-z_.][A-Za-z0-9_.]* <- function", lines)
  end <- after[after > start]
  end <- if (length(end)) end[1] - 1L else length(lines)
  validator <- paste(lines[start:end], collapse = "\n")

  sc <- .rce_column_schema()
  declares <- names(Filter(function(r) !is.null(r$allowed) && !all(is.na(r$allowed)), sc))

  not_validated <- Filter(function(nm) !grepl(nm, validator, fixed = TRUE), declares)

  # Pinned, with the reason each is outstanding.
  known <- c("estDynamics",        # a control-sheet scalar, checked in model_config()
             "Estimate_index_sd",  # read numerically before validation runs
             "Estimate_catch_sd",
             "Diet_distribution")  # per-species, resolved to an integer early
  testthat::expect_setequal(not_validated, known)
})
