# The DSEM is extended to growth and natural mortality by writing its latent
# states into the EXISTING linkage offset tensors, so growth and M gain no
# second time-variation mechanism. That rests entirely on the order three
# things appear in the template:
#
#   1. the linkage fill   zeroes growth_linkage_offset / M_linkage_offset and
#                         adds the linkage-table rows
#   2. section 5.5b       calculate_dsem(), which will ADD the latent states to
#                         those same tensors
#   3. the reads          growth_parameters (section 5.6) and M (inside the
#                         niter loop)
#
# A fill that moved below 5.5b would zero the latents away; a read that moved
# above it would use the pre-DSEM tensor. Either silently drops the growth or M
# process from the model while the fit converges and reports an ABC. Nothing in
# the build enforces the order, and the 5.64.0 merge of `dev` into this line put
# four conflict hunks into this file, so the next merge gets this guard instead
# of a re-derivation. It cannot be a numeric test until a growth latent exists.
#
# The REPORT assertions pin the placement between 5.5b and the reads. TMB binds
# the reported copy at the call site, so a REPORT above 5.5b hands back a
# pre-DSEM snapshot that disagrees with growth_parameters and M.

cpp_lines_no_comments <- function() {
  cpp <- c("src/TMB/ceattle.cpp",
           testthat::test_path("..", "..", "src", "TMB", "ceattle.cpp"))
  cpp <- cpp[file.exists(cpp)]
  testthat::skip_if(length(cpp) == 0, "src/TMB/ceattle.cpp not available")
  src <- readLines(cpp[1], warn = FALSE)
  # Blank `//` lines, keeping the line count so a reported position still
  # matches the file. Needed because "calculate_dsem()" and "5.5b" both appear
  # in the surrounding prose; without it the scan would match a comment and
  # could pass while the code order was wrong.
  src[grepl("^[[:space:]]*//", src)] <- ""
  src
}

# One match per anchor, except the two fill calls (log-link and natural-scale
# tensors). Asserted so a future duplicate fails here rather than letting the
# position checks silently read the wrong occurrence.
anchor <- function(src, pattern, expected_n = 1L) {
  hits <- grep(pattern, src)
  testthat::expect_length(hits, expected_n)
  hits
}

test_that("the linkage offset fill precedes the DSEM section, which precedes every read", {
  src <- cpp_lines_no_comments()

  dsem <- anchor(src, "calculate_dsem[(]")

  g_zero <- anchor(src, "growth_linkage_offset[.]setZero[(][)]")
  g_fill <- anchor(src, "rceattle_apply_growth_linkages[(]", 2L)
  g_read <- anchor(src, "growth_linkage_offset[(]sp, sex, yr, par[)]")

  m_zero <- anchor(src, "M_linkage_offset[.]setZero[(][)]")
  m_fill <- anchor(src, "rceattle_apply_M_linkages[(]", 2L)
  m_read <- anchor(src, "M_linkage_offset[(]sp, sex, age, yr[)]")

  # Fill before 5.5b: a latent written into the tensor must not be zeroed or
  # overwritten afterwards.
  expect_lt(g_zero, dsem)
  expect_lt(max(g_fill), dsem)
  expect_lt(m_zero, dsem)
  expect_lt(max(m_fill), dsem)

  # 5.5b before every read: growth_parameters and M must see the tensor the
  # DSEM has already contributed to.
  expect_lt(dsem, g_read)
  expect_lt(dsem, m_read)

  # setZero before its own fill, or the linkage rows are discarded.
  expect_lt(g_zero, min(g_fill))
  expect_lt(m_zero, min(m_fill))
})

test_that("the offset REPORTs sit below the DSEM section, above every read", {
  src <- cpp_lines_no_comments()

  dsem   <- anchor(src, "calculate_dsem[(]")
  g_read <- anchor(src, "growth_linkage_offset[(]sp, sex, yr, par[)]")
  m_read <- anchor(src, "M_linkage_offset[(]sp, sex, age, yr[)]")

  for (p in c("REPORT[(]growth_linkage_offset[)]",
              "REPORT[(]growth_linkage_offset_nat[)]",
              "REPORT[(]M_linkage_offset[)]",
              "REPORT[(]M_linkage_offset_nat[)]")) {
    expect_gt(anchor(src, p), dsem)
  }

  # Above the reads as well, so the reported tensor is the one the dynamics
  # used -- a REPORT below a read would be the same value today but stops
  # being a snapshot of the input to that read.
  expect_lt(anchor(src, "REPORT[(]growth_linkage_offset[)]"), g_read)
  expect_lt(anchor(src, "REPORT[(]M_linkage_offset[)]"), m_read)
})
