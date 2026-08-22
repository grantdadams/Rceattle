# Every vignette is `eval = FALSE`, so R CMD check renders them without
# executing a line. An API break in a vignette is therefore invisible until a
# user copies the code and it fails -- and the vignettes are the only
# user-facing documentation of the switch codes and units. This is the cheap
# half of the guard: parse the code and check that every Rceattle call names an
# exported function and passes arguments that function actually has.
#
# It cannot catch return-shape drift (mse_summary() changing from a data frame
# to a list). That needs the chunks executed; see the vignette CI job.

vignette_files <- function() {
  dirs <- c("vignettes", testthat::test_path("..", "..", "vignettes"))
  d <- dirs[dir.exists(dirs)]
  if (!length(d)) return(character(0))
  c(list.files(d[1], pattern = "\\.Rmd$", full.names = TRUE),
    list.files(file.path(d[1], "articles"), pattern = "\\.Rmd$", full.names = TRUE))
}

# Every call in an expression tree, as (function name, argument names).
calls_in <- function(expr, acc = list()) {
  if (is.call(expr)) {
    fn <- expr[[1]]
    nm <- if (is.name(fn)) as.character(fn)
          else if (is.call(fn) && identical(as.character(fn[[1]]), "::"))
            as.character(fn[[3]]) else NA_character_
    if (!is.na(nm)) {
      args <- names(as.list(expr))[-1]
      acc[[length(acc) + 1L]] <- list(fn = nm, args = if (is.null(args)) character(0) else args)
    }
    for (el in as.list(expr)[-1]) if (!missing(el)) acc <- calls_in(el, acc)
  }
  acc
}

test_that("vignette code parses", {
  files <- vignette_files()
  testthat::skip_if(length(files) == 0, "vignettes not available")
  testthat::skip_if_not_installed("knitr")
  for (f in files) {
    code <- knitr::purl(f, output = tempfile(fileext = ".R"), quiet = TRUE, documentation = 0L)
    testthat::expect_no_error(parse(code))
  }
})

test_that("vignettes call only exported functions, with arguments those functions have", {
  files <- vignette_files()
  testthat::skip_if(length(files) == 0, "vignettes not available")
  testthat::skip_if_not_installed("knitr")

  ns       <- asNamespace("Rceattle")
  exported <- getNamespaceExports(ns)
  internal <- setdiff(ls(ns, all.names = TRUE), exported)

  bad_fn <- character(0)
  bad_arg <- character(0)

  for (f in files) {
    code <- knitr::purl(f, output = tempfile(fileext = ".R"), quiet = TRUE, documentation = 0L)
    exprs <- tryCatch(parse(code), error = function(e) NULL)
    if (is.null(exprs)) next
    seen <- unlist(lapply(exprs, calls_in), recursive = FALSE)

    for (cl in seen) {
      # An internal function a vignette tells users to call is a promise they
      # cannot keep without :::.
      if (cl$fn %in% internal && !cl$fn %in% exported) {
        bad_fn <- c(bad_fn, paste0(basename(f), ": ", cl$fn, "()"))
        next
      }
      if (!cl$fn %in% exported) next          # base/other packages: not ours

      fm <- names(formals(get(cl$fn, envir = ns)))
      if (!length(fm) || "..." %in% fm) next  # dots absorb anything
      supplied <- cl$args[nzchar(cl$args)]
      unknown  <- setdiff(supplied, fm)
      if (length(unknown)) {
        bad_arg <- c(bad_arg,
                     paste0(basename(f), ": ", cl$fn, "(", paste(unknown, collapse = ", "), ")"))
      }
    }
  }

  testthat::expect_length(bad_fn, 0)
  testthat::expect_length(bad_arg, 0)
})
