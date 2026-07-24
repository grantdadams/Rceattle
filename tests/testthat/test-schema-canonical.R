# Canonical workbook-column schema (R/0-column_schema.R): faithful-reproduction
# and integrity guards. These prove the single-source schema reproduces the
# hand-maintained copies BEFORE any consumer is switched over to it.

# Compare descriptions on their meaningful text: collapse CR/LF and runs of
# spaces so a drift guard flags wording changes, not cosmetic whitespace.
norm_desc <- function(x) {
  x <- gsub("[\r\n]+", "\n", x)
  x <- gsub("[ \t]+", " ", x)
  x <- gsub(" ?\n ?", "\n", x)
  trimws(x)
}

test_that("schema is structurally sound", {
  schema <- .rce_column_schema()
  nms <- vapply(schema, function(r) r$name, character(1))

  # Unique canonical names; names() keyed correctly.
  expect_equal(anyDuplicated(nms), 0L)
  expect_equal(names(schema), unname(nms))

  # Valid sheets.
  sheets <- vapply(schema, function(r) r$sheet, character(1))
  expect_true(all(sheets %in% c("control", "fleet_control", "bioenergetics_control")))

  # No alias collides with a canonical name (would make a rename ambiguous).
  aliases <- unlist(lapply(schema, function(r) r$aliases))
  expect_length(intersect(aliases, nms), 0)

  # No alias is shared by two columns (would make an upgrade ambiguous).
  expect_equal(anyDuplicated(aliases), 0L)

  # The orphan set is pinned: these are documented historically but read
  # nowhere and slated for removal. Promoting one to "live" (or adding a new
  # orphan) must be a deliberate, visible change -- pinning catches it, so the
  # file's "explicit, testable record of the drop" claim actually holds.
  orphan <- vapply(Filter(function(r) r$status == "orphan", schema),
                   function(r) r$name, character(1))
  expect_setequal(unname(orphan), c("R_sexr", "Catch_units", "Sex"))

  # Orphans never leak into the live accessors used by the generators.
  expect_length(intersect(.rce_schema_names("fleet_control"), orphan), 0)
  expect_length(intersect(.rce_schema_names("control"), orphan), 0)
})

test_that("meta_data generator reproduces the bundled meta_data_names.xlsx", {
  schema <- .rce_column_schema()
  gen <- .rce_build_meta_data_df()

  xlsx <- as.data.frame(readxl::read_xlsx(
    system.file("extdata", "meta_data_names.xlsx", package = "Rceattle")))
  colnames(xlsx) <- c("sheet", "col", "desc")

  # --- Column rows: map each side to a canonical key and compare descriptions.
  alias_to_canon <- list()
  for (r in schema) for (a in r$aliases) alias_to_canon[[a]] <- r$name
  to_canon <- function(col) {
    if (is.na(col)) return(NA_character_)
    if (!is.null(alias_to_canon[[col]])) alias_to_canon[[col]] else col
  }

  gen_cols <- gen[!is.na(gen$`Column/row name`), ]
  xlsx_cols <- xlsx[!is.na(xlsx$col) & !is.na(xlsx$desc), ]

  gen_map <- setNames(gen_cols$Description,
                      vapply(gen_cols$`Column/row name`, to_canon, character(1)))
  xlsx_map <- setNames(xlsx_cols$desc,
                       vapply(xlsx_cols$col, to_canon, character(1)))

  # Same set of documented columns (by canonical name).
  expect_setequal(names(gen_map), names(xlsx_map))

  # Identical descriptions for every documented column.
  for (nm in names(xlsx_map))
    expect_equal(norm_desc(gen_map[[nm]]), norm_desc(xlsx_map[[nm]]),
                 info = paste("meta description for", nm))

  # --- Sheet-header rows (Column/row name is NA, Description names the sheet).
  gen_hdr <- gen[is.na(gen$`Column/row name`) & !is.na(gen$Description), ]
  xlsx_hdr <- xlsx[is.na(xlsx$col) & !is.na(xlsx$desc), ]
  expect_setequal(gen_hdr$`Sheet name`, xlsx_hdr$sheet)
  gh <- setNames(gen_hdr$Description, gen_hdr$`Sheet name`)
  xh <- setNames(xlsx_hdr$desc, xlsx_hdr$sheet)
  for (nm in names(xh))
    expect_equal(norm_desc(gh[[nm]]), norm_desc(xh[[nm]]),
                 info = paste("sheet header for", nm))
})

test_that("every bundled meta_data column maps to a schema name or alias", {
  schema <- .rce_column_schema()
  allkeys <- unlist(lapply(schema, function(r) c(r$name, r$aliases)))
  xlsx <- as.data.frame(readxl::read_xlsx(
    system.file("extdata", "meta_data_names.xlsx", package = "Rceattle")))
  cols <- xlsx[[2]]
  cols <- cols[!is.na(cols)]
  # The lone "Sex" legend row (col holds the legend text, not a column) is not a
  # real column; every other documented column must resolve.
  cols <- cols[!grepl("combined|female|male", cols)]
  expect_length(setdiff(cols, allkeys), 0)
})
