# Regenerate inst/extdata/meta_data_names.xlsx from the canonical column schema.
#
# The embedded meta_data documentation sheet is generated from
# R/0-column_schema.R (`.rce_build_meta_data_df()`), so it can never drift from
# the schema. Re-run this after any change to the schema's docs or column set:
#
#   export PATH=/usr/bin:$PATH
#   Rscript data-raw/regenerate_meta_xlsx.R
#
# The drift-guard test in tests/testthat/test-schema-canonical.R asserts the
# committed xlsx matches the generator, so CI fails if this was not re-run.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

df <- Rceattle:::.rce_build_meta_data_df()
out <- file.path("inst", "extdata", "meta_data_names.xlsx")
writexl::write_xlsx(stats::setNames(list(df), "meta_data"), out)
cat("Wrote", nrow(df), "rows to", out, "\n")
