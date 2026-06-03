# Regression note for ceattle_v01_11.cpp:1989 (age_first_selected guard in
# marginal length-comp prediction loop):
#
# The marginal length-comp prediction at pred_CAAL must skip ages below
# age_first_selected to mirror SS3 pattern-10's sel_a[0] = 0 convention.
# Without the guard, age-0 N x sel_at_length_pop(small L) x ALK(age=0,
# small L) bleeds into the smallest data bin via the ALK tail. On the
# GOA Pcod 2024 model this inflated Srv 2023 bin 4.5 from 0.003 (SS3) to
# 0.067 (Rce) before the fix -- a 23x error in the predicted proportion.
#
# Reproducing this bug in a synthetic data unit test is fragile because
# it requires:
#   (a) length-based DoubleLogistic / DoubleNormal selectivity (not the
#       synthetic helper's default age-based logistic),
#   (b) length-comp rows (Age0_Length1 == 1) -- the helper produces only
#       age-comp rows (Age0_Length1 == 0),
#   (c) parametric growth so growth_matrix_pop is non-zero,
#   (d) populated lengths_pop and the pop-grid pipeline.
#
# This combination is not currently constructible via the synthetic
# helper. The fix is empirically verified on the GOA Pcod 2024 model
# (see ../../../Rceattle-models/GOA cod/SS3_Rceattle_parity_status.md):
# Srv LenComp NLL dropped from +44.84 to +0.09 vs SS3 after this fix.

testthat::test_that("age_first_selected guard exists in cpp length-comp loop", {
  # Lightweight grep-based check on the cpp source: the marginal length-
  # comp prediction loop must contain the age < amin_lcomp continue guard.
  # If someone removes it (regression), this fires.
  # test_dir() runs each test file with wd = file's directory (e.g.
  # tests/testthat/tests-Selectivity/). Source cpp is 3 levels up.
  cpp_candidates <- c(
    file.path("..", "..", "..", "src", "TMB", "ceattle_v01_11.cpp"),
    file.path("..", "..", "src", "TMB", "ceattle_v01_11.cpp"),
    file.path("src", "TMB", "ceattle_v01_11.cpp")
  )
  exists_mask <- file.exists(cpp_candidates)
  testthat::skip_if(!any(exists_mask), "cpp source not found")
  cpp <- readLines(cpp_candidates[which(exists_mask)[1]])

  # Find the pred_comp section (line range with both is_len_based and
  # age_first_selected references)
  has_amin_lcomp <- any(grepl("amin_lcomp", cpp))
  has_continue   <- any(grepl("if \\(age < amin_lcomp\\) continue", cpp))
  testthat::expect_true(has_amin_lcomp,
                        label = "amin_lcomp variable present in cpp")
  testthat::expect_true(has_continue,
                        label = "age < amin_lcomp continue guard present in cpp")
})
