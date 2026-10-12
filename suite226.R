suppressMessages(pkgload::load_all(".", quiet = TRUE))
r <- devtools::test(reporter = "summary", stop_on_failure = FALSE)
d <- as.data.frame(r)
cat(sprintf("SUITE226 files=%d pass=%d fail=%d err=%d skip=%d warn=%d\n",
            length(unique(d$file)), sum(d$passed), sum(d$failed), sum(d$error),
            sum(d$skipped), sum(d$warning)))
bad <- d[d$failed > 0 | d$error > 0, c("file", "test", "failed", "error")]
if (nrow(bad)) print(bad) else cat("SUITE226 no failures\n")
