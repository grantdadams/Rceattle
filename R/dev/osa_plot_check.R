# Render-check the OSA / process residual plots. Run via:
#   Rscript R/dev/osa_plot_check.R
suppressMessages(devtools::load_all(".", compile = FALSE, quiet = TRUE))

ok <- TRUE
try_plot <- function(label, expr) {
  f <- tempfile(fileext = ".png")
  res <- tryCatch({ grDevices::png(f); print(force(expr)); grDevices::dev.off(); TRUE },
                  error = function(e) { try(grDevices::dev.off(), silent = TRUE)
                    cat("  ERROR:", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("[%s] %s\n", if (res) "PASS" else "FAIL", label))
  if (!res) ok <<- FALSE
}

data("GOApollock")
fit <- Rceattle::fit_mod(
  data_list = GOApollock, inits = NULL, file = NULL, estimateMode = 1,
  random_rec = FALSE, msmMode = 0,
  fit_control = fit_control(verbose = 0, phase = TRUE, getsd = TRUE))

try_plot("plot(osa aggregate)", plot(osa_residuals(fit, source = c("index", "catch"))))
try_plot("plot_indexresidual(residual_type='osa')",
         plot_indexresidual(fit, residual_type = "osa"))
try_plot("plot_comp(residual_type='osa')", plot_comp(fit, residual_type = "osa"))
try_plot("plot(process residuals)", plot(process_residuals(fit, "recruitment")))

cat(sprintf("\n==== OVERALL: %s ====\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1)
