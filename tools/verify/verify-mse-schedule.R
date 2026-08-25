# tools/verify/verify-mse-schedule.R
# End-to-end net for the two run_mse() arguments that describe an irregular
# assessment design: an explicit vector of assessment years, and a catch
# multiplier indexed by year and species.
#
# test-mse-schedule-and-catch-mult.R covers the two helpers in isolation. What
# it cannot show is that the values reach the loop unchanged -- that a schedule
# supplied as years produces the SAME trajectory as the equivalent period, and
# that a year-indexed multiplier moves catch in the years it names and in no
# others. Both need a real closed loop, so they live here rather than in the
# suite.
#
# Fixture: BS2017SS truncated to projyr 2020, the verify-mse-repro.R recipe.
# endyr is 2017, so assessment_period = 1 is exactly c(2018, 2019, 2020) and
# the two forms must agree to the last digit. Three projection years is the
# shortest horizon that leaves room to skip an assessment and still have one on
# either side of the gap.
#
# Determinism: cores = 1 and each sim self-seeds from seed + sim.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-mse-schedule.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# Proj_F_proportion must sum to 1 across each species' fishery fleets or
# run_mse() refuses to project. Same helper as verify-mse-repro.R.
activate_proj <- function(dl) {
  fc <- dl$fleet_control
  p <- rep(0, nrow(fc))
  is_fish <- fc$Fleet_type == 1
  for (s in unique(fc$Species[is_fish])) {
    idx <- which(is_fish & fc$Species == s)
    p[idx] <- 1 / length(idx)
  }
  dl$fleet_control$Proj_F_proportion <- p
  dl
}

data(BS2017SS)
BS2017SS$projyr <- 2020
BS2017SS <- activate_proj(BS2017SS)

om <- suppressWarnings(Rceattle::fit_mod(
  data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
  random_rec = FALSE, msmMode = 0,
  fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)))

em <- suppressWarnings(Rceattle::fit_mod(
  data_list = BS2017SS, inits = om$estimated_params, estimateMode = 2,
  HCR = build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35, Plimit = 0.2,
                  Alpha = 0.05),
  msmMode = 0, fit_control = fit_control(getsd = FALSE, verbose = 0)))

ENDYR <- om$data_list$endyr
NSPP  <- om$data_list$nspp
stopifnot(ENDYR == 2017)

run <- function(...) suppressMessages(suppressWarnings(Rceattle::run_mse(
  om = om, em = em, nsim = 1, sampling_period = 1, simulate_data = TRUE,
  sample_rec = TRUE, regenerate_past = FALSE, seed = 666, cores = 1, ...)))

# Realized catch in the operating model, by year and species. This is what a
# reduction is supposed to move, and the only thing that distinguishes the
# scenarios below.
om_catch <- function(mse) {
  do.call(rbind, lapply(names(mse), function(nm) {
    om_i <- mse[[nm]]$OM
    if (is.null(om_i)) return(NULL)
    cd <- om_i$data_list$catch_data
    keep <- cd$Year > ENDYR & cd$Year <= om_i$data_list$endyr
    data.frame(sim = nm, Year = cd$Year[keep], Species = cd$Species[keep],
               Catch = cd$Catch[keep])
  }))
}

report <- function(label, ok, detail = "") {
  cat(sprintf("%-58s %s%s\n", label, if (isTRUE(ok)) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("  [", detail, "]") else ""))
  invisible(isTRUE(ok))
}

results <- logical(0)

# --- 1. An explicit schedule equals the period it reproduces ---------------
# endyr 2017, projyr 2020: period 1 IS c(2018, 2019, 2020). Anything but an
# exact match means the vector path reaches the loop differently.
mse_period <- run(assessment_period = 1)
mse_vector <- run(assessment_period = c(2018, 2019, 2020))

c_period <- om_catch(mse_period)
c_vector <- om_catch(mse_vector)
results["equivalence"] <- report(
  "explicit years == the equivalent period",
  isTRUE(all.equal(c_period, c_vector, tolerance = 0)),
  paste0("max |diff| = ",
         format(max(abs(c_period$Catch - c_vector$Catch)), digits = 3)))

# The assessment nodes are named by year, so the schedule is visible in the
# returned object rather than only in the caller's script.
results["names"] <- report(
  "assessment nodes are named by their year",
  identical(names(mse_vector$Sim_1$EM),
            c("EM", paste0("OM_Sim_1. EM_yr_", c(2018, 2019, 2020)))))

# --- 2. A skipped assessment is a shorter schedule, not a longer period ----
# Assess in 2018 and 2020; the 2019 assessment is missed, so 2019 is fished on
# 2018's advice and the 2020 assessment covers a two-year interval.
mse_skip <- run(assessment_period = c(2018, 2020))
results["skip"] <- report(
  "a skipped assessment runs one fewer assessment, over the same horizon",
  length(mse_skip$Sim_1$EM) == 3L &&
    mse_skip$Sim_1$OM$data_list$endyr == 2020)

# It must not be the same trajectory as assessing in every year.
c_skip <- om_catch(mse_skip)
results["skip_differs"] <- report(
  "skipping an assessment changes the trajectory",
  !isTRUE(all.equal(c_skip$Catch, c_period$Catch)))

# WHICH years the gap moves. By the advice mapping alone it should be 2020 and
# 2020 only: the assessment in year Y sets catch for Y+1 onward, so dropping
# the 2019 assessment leaves 2020 on 2018's advice, while 2019 is set by the
# 2018 assessment in both schedules.
#
# It is not what happens. 2018 is identical, but 2019 moves as well -- 2880119
# against 2941211 on this fixture, 2.1% -- on advice that is identical by
# construction. The cause is the operating model's horizon shortening: at each
# assessment the OM is trimmed to the NEXT assessment year, so sim_mod() draws
# over a different number of projection rows depending on the schedule, and the
# observation draws diverge from the first assessment onward. Two scenarios
# differing only in their schedule are therefore NOT on common random numbers,
# and a paired comparison between them carries an observation-error difference
# it did not ask for.
#
# Pinned as it behaves today. If the draw is ever made independent of the
# schedule, this check fails and points at the reason.
i18 <- c_period$Year == 2018
i19 <- c_period$Year == 2019
i20 <- c_period$Year == 2020
results["gap_year_2018_clean"] <- report(
  "the year before the gap is bit-identical",
  isTRUE(all.equal(c_skip$Catch[i18], c_period$Catch[i18], tolerance = 0)))

results["gap_year_advice"] <- report(
  "the year the missed assessment would have set moves",
  !isTRUE(all.equal(c_skip$Catch[i20], c_period$Catch[i20])))

results["gap_crn_broken"] <- report(
  "KNOWN: a same-advice year moves too (draws are schedule-dependent)",
  !isTRUE(all.equal(c_skip$Catch[i19], c_period$Catch[i19])),
  paste0("2019 ratio = ",
         format(sum(c_skip$Catch[i19]) / sum(c_period$Catch[i19]),
                digits = 4)))

# It must also not be the biennial period it superficially resembles: a period
# of 2 from 2017 is c(2019), a different schedule over a different horizon.
results["skip_not_period"] <- report(
  "a skipped schedule is not the equivalent period",
  !identical(.mse_assess_years(c(2018, 2020), ENDYR, 2020),
             suppressWarnings(.mse_assess_years(2, ENDYR, 2020))))

# The trailing years of a schedule that stops short carry NA catch, which is
# what the short-schedule warning exists to report. Asserted against the run
# rather than the helper, because it is a property of the loop's fill window.
mse_short <- run(assessment_period = c(2018, 2019))
short_cd  <- mse_short$Sim_1$OM$data_list$catch_data
results["short_schedule_na"] <- report(
  "a schedule stopping short leaves the trailing years' catch NA",
  all(is.na(short_cd$Catch[short_cd$Year == 2020])) &&
    !any(is.na(short_cd$Catch[short_cd$Year %in% 2018:2019])))

warned <- character(0)
withCallingHandlers(
  .mse_assess_years(c(2018, 2019), ENDYR, 2020),
  warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
results["short_schedule_warns"] <- report(
  "and the short schedule is warned about",
  any(grepl("catch in 2020 is never set", warned)))

# --- 3. A year-indexed catch multiplier moves only the years it names ------
buffer <- expand.grid(Year = 2020, Species = seq_len(NSPP))
buffer$mult <- 0.5

mse_buffer <- run(assessment_period = 1, catch_mult = buffer)
c_buffer   <- om_catch(mse_buffer)

# 2018 and 2019 precede the reduction, and the operating model's dynamics in
# those years are settled before the 2020 catch is ever written, so they must
# be untouched.
i_before <- c_period$Year < 2020
results["year_before"] <- report(
  "the years before the reduction are bit-identical",
  isTRUE(all.equal(c_buffer$Catch[i_before], c_period$Catch[i_before],
                   tolerance = 0)))

# 2020 is halved, up to the exploitable-biomass limit -- which only ever binds
# downward, so the reduced catch can never come out above the baseline.
i_of <- c_period$Year == 2020
results["year_of"] <- report(
  "the reduction year is cut",
  all(c_buffer$Catch[i_of] <= c_period$Catch[i_of]) &&
    any(c_buffer$Catch[i_of] < c_period$Catch[i_of]),
  paste0("mean ratio = ",
         format(mean(c_buffer$Catch[i_of] / c_period$Catch[i_of]), digits = 4)))

# --- 4. A vector multiplier still cuts every year --------------------------
# The contrast that motivates the data.frame form: the same 0.5, supplied as a
# per-species vector, is a permanent harvest cut and takes 2018 with it.
mse_flat <- run(assessment_period = 1, catch_mult = rep(0.5, NSPP))
c_flat   <- om_catch(mse_flat)

results["vector_every_year"] <- report(
  "a vector multiplier cuts the earlier years too",
  all(c_flat$Catch[i_before] < c_period$Catch[i_before]),
  paste0("mean ratio pre-2020 = ",
         format(mean(c_flat$Catch[i_before] / c_period$Catch[i_before]),
                digits = 4)))

cat("\n")
if (all(results)) {
  cat("verify-mse-schedule: all", length(results), "checks passed\n")
} else {
  cat("verify-mse-schedule: FAILED --",
      paste(names(results)[!results], collapse = ", "), "\n")
  quit(status = 1)
}
