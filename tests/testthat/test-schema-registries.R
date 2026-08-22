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

test_that("every model-level switch in the config schema carries its full allowed set", {
  # .rce_config_schema() is a second, hand-written dictionary of the model-level
  # switches, used to comment the YAML that save_config() writes. It hardcodes
  # some allowed sets instead of reading the map, which is how estimateMode's
  # DebugOptimize went missing from the comments while being a real mode.
  dict <- .rce_config_schema()

  pairs <- list(
    estimateMode = estimateMode_map,
    initMode     = initMode_map,
    suitMode     = suitMode_map,
    msmMode      = msmMode_map
  )

  for (nm in names(pairs)) {
    entry <- dict[[nm]]
    testthat::expect_false(is.null(entry),
                           info = paste(nm, "absent from .rce_config_schema()"))
    if (is.null(entry)) next

    # The value set may be carried as `allowed`, or written into `doc`.
    txt <- paste(c(as.character(entry$allowed), as.character(entry$doc)), collapse = " ")
    for (canon in names(pairs[[nm]])) {
      testthat::expect_true(
        grepl(canon, txt, fixed = TRUE),
        info = paste0(nm, ": '", canon, "' is a valid value the YAML comment never mentions"))
    }
  }
})

test_that("the parameter dictionary covers every parameter the template declares", {
  # .PAR_INFO maps the TMB parameter names to what a user actually chose to
  # estimate; error messages, diagnostics and output summaries read it. A
  # parameter missing from it surfaces to the user under its raw internal name.
  cpp <- c("src/TMB/ceattle.cpp", testthat::test_path("..", "..", "src", "TMB", "ceattle.cpp"))
  cpp <- cpp[file.exists(cpp)]
  testthat::skip_if(length(cpp) == 0, "src/TMB/ceattle.cpp not available")
  src <- paste(readLines(cpp[1], warn = FALSE), collapse = "\n")

  # Declared parameters, ignoring anything inside a /* ... */ comment block --
  # ceattle.cpp keeps several retired predation parameters commented out.
  src <- gsub("/\\*.*?\\*/", "", src, perl = TRUE)
  m <- regmatches(src, gregexpr("PARAMETER[A-Z_]*\\(\\s*[A-Za-z0-9_]+\\s*\\)", src, perl = TRUE))[[1]]
  declared <- unique(gsub(".*\\(\\s*([A-Za-z0-9_]+)\\s*\\).*", "\\1", m))
  testthat::expect_gt(length(declared), 20)

  # Parameters belonging to features the model declares but does not implement.
  # They exist so the template compiles; no user can estimate one, so no error
  # message needs to name them. Pinned, so the list cannot quietly grow.
  stubbed <- c(
    "index_q_pow",                                  # Catchability = PowerEquation (rejected in data_check)
    "logH_1", "logH_1a", "logH_1b", "logH_2", "logH_3", "H_4"  # Kinzey-Punt predation, msmMode 3-9
  )

  known <- c(.PAR_INFO$internal, stubbed)
  missing_from_dict <- setdiff(declared, known)
  testthat::expect_length(missing_from_dict, 0)

  # And the exemption must stay honest: a stubbed parameter that gains a real
  # implementation should be documented, not left on this list.
  testthat::expect_length(intersect(.PAR_INFO$internal, stubbed), 0)
})


test_that("the hand-coded defaults still equal the schema's", {
  # Some defaults cannot be schema-driven: the bioenergetics scalars are
  # per-species (`rep(x, nspp)`) and only filled when the model is
  # multispecies, which the schema's flat `default` field cannot express. They
  # stay imperative in switch_check(), and the schema rows carry the same
  # values as documentation. Nothing kept the two in step; this does.
  #
  # A mismatch means the meta sheet and the field dictionary are telling users
  # one default while switch_check() applies another.
  sc <- .rce_column_schema()

  hand_coded <- list(
    Ceq = 1L, Cindex = 0L, Pvalue = 1, fday = 365,
    CA = 0, CB = 0, Qc = 0, Tco = 0, Tcm = 0, Tcl = 0, CK1 = 0, CK4 = 0
  )

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
