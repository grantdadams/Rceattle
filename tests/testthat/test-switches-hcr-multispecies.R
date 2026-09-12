# NPFMC Tier 3 (HCR 5) reads SBF, which under msmMode > 0 is built on the
# projection's realized M2 rather than an equilibrium M2. validate_switches()
# refuses it there, along with HCRs 4 and 7. These tests pin that refusal, so
# allowing HCR 5 in multispecies mode cannot happen without a failing test.

.msm_hcr_guard   <- "work in multi-species mode"
.msm_hcr_allowed <- c("NoFishing", "CMSY", "ConstantF", "ConstantFSSB", "PFMC")

test_that("validate_switches() refuses exactly the HCRs outside the multispecies set", {
  dl <- Rceattle::BS2017MS
  dl$msmMode <- 1
  for (h in names(hcr_map)) {
    dl$HCR <- h
    refused <- any(grepl(.msm_hcr_guard, validate_switches(dl), fixed = TRUE))
    expect_identical(refused, !h %in% .msm_hcr_allowed, info = h)
  }
})

test_that("fit_mod() stops on NPFMC in multispecies mode, by code or by name", {
  for (h in list(5L, "NPFMC")) {
    expect_error(
      suppressWarnings(suppressMessages(Rceattle::fit_mod(
        data_list = Rceattle::BS2017MS, msmMode = 1, estimateMode = 3,
        HCR = Rceattle::build_hcr(HCR = h),
        fit_control = Rceattle::fit_control(verbose = 0)))),
      .msm_hcr_guard, fixed = TRUE, info = as.character(h))
  }
})
