# Catchability = "AR1" was the QAR1 form (Rogers et al. 2024). It never worked:
# build_map() gates the log-q deviates on Time_varying_q %in% c("IID","AR1",
# "RandomWalk"), but under QAR1 that column holds an env_data COLUMN INDEX
# rather than a mode, so the deviates were never freed and q came back
# constant. Nothing errored -- the fit ran and reported a time-invariant q
# where the user asked for a time-varying one, and its objective carried the
# AR1 normalizing constant plus the environmental index fitted as noise about
# zero, so it was not comparable with any other model's.
#
# It is now an ERROR, not a warning: a warned fit still returns a summary()
# that looks ordinary, which is the failure mode this package cannot ship.
# That is the severity 'PowerEquation' already carries in the same function.
#
# The error must fire on that switch and ONLY that switch: `Time_varying_q =
# "AR1"` is a different setting sharing the same string (an AR1 structure on an
# ordinary "Estimated" q) and works correctly. An earlier draft of the
# developer notes named the wrong one.
#
# The replacement it offers must also RUN, and must be the QAR1 model rather
# than a free AR1 on q -- the last test lifts the snippet out of the message,
# fills in its placeholders and evaluates it, so advice that does not execute
# fails here rather than in someone's assessment. GOA pollock 2025 runs that
# same form in production.

# data_check() accumulates and stops once, so the QAR1 text arrives in the
# error. Collect warnings too, to assert it is NOT merely warned about.
check_messages <- function(d) {
  w <- character(0)
  err <- NA_character_
  withCallingHandlers(
    tryCatch(data_check(d), error = function(e) err <<- conditionMessage(e)),
    warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  list(error = err, warnings = w)
}

qar1_error <- function(d) {
  m <- check_messages(d)
  if (!is.na(m$error) && grepl("QAR1", m$error, fixed = TRUE)) m$error else NA_character_
}

base_data <- function() {
  suppressMessages(suppressWarnings(switch_check(make_test_data(nyrs = 4, nages = 3))))
}

test_that("Catchability = 'AR1' errors, and names the fleet", {
  d <- base_data()
  d$fleet_control$Catchability[1] <- "AR1"
  qe <- qar1_error(d)

  testthat::expect_false(is.na(qe))
  # A user must be able to act on it: which fleet, what was actually happening,
  # and what to write instead.
  testthat::expect_match(qe, d$fleet_control$Fleet_name[1], fixed = TRUE)
  testthat::expect_match(qe, "~ ar1(1 | Year)", fixed = TRUE)
  testthat::expect_match(qe, "qFun", fixed = TRUE)
  # `observe` is the whole difference between QAR1 and a free AR1 on q, so the
  # replacement is the wrong model without it.
  testthat::expect_match(qe, "observe", fixed = TRUE)
  testthat::expect_match(qe, "obs_sd", fixed = TRUE)
})

test_that("it is an error, not a warning -- and quiet_data_check cannot mute it", {
  d <- base_data()
  d$fleet_control$Catchability[1] <- "AR1"
  m <- check_messages(d)

  testthat::expect_false(any(grepl("QAR1", m$warnings, fixed = TRUE)))
  testthat::expect_true(grepl("QAR1", m$error, fixed = TRUE))
  # suppressWarnings is what fit_mod(quiet_data_check = TRUE) applies; an error
  # must survive it, or the switch stays silently usable.
  testthat::expect_error(suppressWarnings(data_check(d)), "QAR1")
})

test_that("the error names the env_data series the fleet points at", {
  d <- base_data()
  d$fleet_control$Catchability[1] <- "AR1"
  d$fleet_control$Time_varying_q[1] <- 2
  d$env_data <- data.frame(Year = d$styr:d$endyr, first_cov = 0, second_cov = 0)

  # Time_varying_q indexes the covariates, i.e. env_data without its Year
  # column, so 2 is `second_cov`.
  testthat::expect_match(qar1_error(d), "'second_cov'", fixed = TRUE)
})

test_that("the error falls back to a placeholder when the series cannot be named", {
  d <- base_data()
  d$fleet_control$Catchability[1] <- "AR1"
  d$fleet_control$Time_varying_q[1] <- 9   # out of range
  d$env_data <- data.frame(Year = d$styr:d$endyr, only_cov = 0)

  testthat::expect_match(qar1_error(d), "<env_data column>", fixed = TRUE)
})

test_that("two QAR1 fleets stay aligned with their series when one index is unusable", {
  # `covs[c(0, 1)]` returns ONE element, and a negative index drops elements
  # instead of returning one per fleet -- either way the fleet names would be
  # pasted against the wrong series. Every fleet must get its own label.
  d <- base_data()
  testthat::skip_if(nrow(d$fleet_control) < 2)
  d$fleet_control$Catchability[1:2] <- "AR1"
  d$fleet_control$Time_varying_q[1:2] <- c(0, 1)   # 0 is not a valid index
  d$env_data <- data.frame(Year = d$styr:d$endyr, only_cov = 0)
  qe <- qar1_error(d)

  testthat::expect_match(qe, d$fleet_control$Fleet_name[1], fixed = TRUE)
  testthat::expect_match(qe, d$fleet_control$Fleet_name[2], fixed = TRUE)
  testthat::expect_match(qe, "<env_data column>", fixed = TRUE)   # fleet 1
  testthat::expect_match(qe, "'only_cov'", fixed = TRUE)          # fleet 2
})

test_that("the integer spelling errors too", {
  d <- base_data()
  d$fleet_control$Catchability[1] <- 6
  testthat::expect_false(is.na(qar1_error(d)))
})

test_that("Time_varying_q = 'AR1' is untouched -- a different, working switch", {
  d <- base_data()
  d$fleet_control$Time_varying_q[1] <- "AR1"
  m <- check_messages(d)
  testthat::expect_true(is.na(qar1_error(d)))
  testthat::expect_false(any(grepl("QAR1", m$warnings, fixed = TRUE)))
})

test_that("a model using neither is silent about QAR1", {
  m <- check_messages(base_data())
  testthat::expect_true(is.na(qar1_error(base_data())))
  testthat::expect_false(any(grepl("QAR1", m$warnings, fixed = TRUE)))
})

test_that("the replacement the error prints parses and evaluates", {
  d <- base_data()
  d$fleet_control$Catchability[1] <- "AR1"
  qe <- qar1_error(d)

  # Lift the build_catchability() call straight out of the message, so this
  # checks the text a user actually copies rather than a restatement of it.
  code <- regmatches(qe, regexpr("(?s)build_catchability\\(.*?\\)\\)\\)", qe, perl = TRUE))
  testthat::expect_length(code, 1L)

  code <- sub("<Fleet_code>", "1", code, fixed = TRUE)
  code <- sub("\"<that fleet's env_data column>\"", "\"qcov\"", code, fixed = TRUE)
  code <- sub("<its measurement SD>", "0.1", code, fixed = TRUE)

  spec <- eval(parse(text = code), envir = asNamespace("Rceattle"))
  testthat::expect_type(spec, "list")

  # And it is the QAR1 model, not a free AR1: the deviates are observed against
  # the named series at a fixed measurement SD.
  q <- spec$linkages$q
  testthat::expect_equal(q$observe, "qcov")
  testthat::expect_equal(q$obs_sd, 0.1)
})
