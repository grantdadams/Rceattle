# dev/bench-fit-decompose.R
# Decompose a fit_mod() call into taping (MakeADFun) vs nlminb evaluation vs
# sdreport, to settle whether the AD taping or the repeated fn/gr evals dominate.
# Read-only. Isolate:
#   estimateMode = 3                       -> data-prep + ONE MakeADFun (taping), no nlminb/sd
#   estimateMode = 1, phase=F, loopnum=1, getsd=F  -> taping + nlminb (one pass)
#   + getsd=T (getJointPrecision default T)         -> + full sdreport
#   + getsd=T, getJointPrecision=F                  -> + sdreport w/o joint precision
# Differences give taping / optimize / sdreport / joint-precision seconds.
#
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript dev/bench-fit-decompose.R

Sys.setenv(OMP_NUM_THREADS = "1")
suppressMessages(pkgload::load_all(".", quiet = TRUE, compile = FALSE))
data(BS2017SS); data(BS2017MS); data(GOA2018SS)

timeit <- function(expr) { t <- proc.time()[["elapsed"]]; v <- force(expr)
  list(sec = proc.time()[["elapsed"]] - t, val = v) }

decompose <- function(label, dl, msm, niter, inits = NULL) {
  ff <- function(...) suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = dl, inits = inits, file = NULL, msmMode = msm, niter = niter,
    suitMode = 0, random_rec = FALSE, ...)))
  # 1. taping: data-prep + 1 MakeADFun, no optimize/sd (obj kept alive at mode 3)
  b3 <- timeit(ff(estimateMode = 3, fit_control = fit_control(getsd = FALSE, verbose = 0)))
  obj <- b3$val$obj; p <- obj$par
  eps <- 1e-6 * (abs(p) + 1); i <- 0L
  pert <- function() { i <<- i + 1L; p + (as.numeric(i %% 2L) - 0.5) * eps }
  ev_fn <- median(replicate(15, timeit(obj$fn(pert()))$sec))
  ev_gr <- median(replicate(15, timeit(obj$gr(pert()))$sec))
  # 2. hindcast-only, single nlminb pass, no phase/sd
  h0 <- timeit(ff(estimateMode = 1, fit_control = fit_control(phase = FALSE, getsd = FALSE, loopnum = 1, verbose = 0)))
  it  <- tryCatch(h0$val$opt$iterations, error = function(e) NA)
  evs <- tryCatch(paste(h0$val$opt$evaluations, collapse = "/"), error = function(e) NA)
  # 3. + full sdreport (getJointPrecision default TRUE)
  hs <- timeit(ff(estimateMode = 1, fit_control = fit_control(phase = FALSE, getsd = TRUE, loopnum = 1, verbose = 0)))
  # 3b. + sdreport without joint precision
  hj <- timeit(ff(estimateMode = 1, fit_control = fit_control(phase = FALSE, getsd = TRUE, getJointPrecision = FALSE, loopnum = 1, verbose = 0)))

  taping   <- b3$sec
  optimize <- h0$sec - b3$sec
  sd_full  <- hs$sec - h0$sec
  sd_noj   <- hj$sec - h0$sec
  jointp   <- hs$sec - hj$sec
  cat(sprintf(
    "%-7s taping=%5.1fs  optimize=%6.1fs (%s iters, f/g evals=%s)  sdreport=%5.1fs (no-joint=%5.1fs, jointprec=%5.1fs)\n         per-eval fn=%5.0fms gr=%5.0fms  nfree=%d niter=%g\n",
    label, taping, optimize, it, evs, sd_full, sd_noj, jointp,
    1000*ev_fn, 1000*ev_gr, length(p), niter))
  invisible(list(label = label, taping = taping, optimize = optimize, sd_full = sd_full,
                 sd_noj = sd_noj, jointp = jointp, ev_fn = ev_fn, ev_gr = ev_gr,
                 nlm_iters = it, nfree = length(p), niter = niter))
}

cat("Decomposing fit_mod wall-clock (hindcast-only, isolating taping vs eval vs sdreport)\n\n")
res <- list()
res$BS_SS  <- decompose("BS_SS",  BS2017SS, 0, 1)
ss_bs  <- suppressWarnings(suppressMessages(Rceattle::fit_mod(BS2017SS, estimateMode = 1, msmMode = 0,
  fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0))))
res$BS_MS  <- decompose("BS_MS",  BS2017MS, 1, 5, inits = ss_bs$estimated_params)
res$GOA_SS <- decompose("GOA_SS", GOA2018SS, 0, 1)
ss_goa <- suppressWarnings(suppressMessages(Rceattle::fit_mod(GOA2018SS, estimateMode = 1, msmMode = 0,
  fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0))))
res$GOA_MS <- decompose("GOA_MS", GOA2018SS, 1, 3, inits = ss_goa$estimated_params)

saveRDS(res, "dev/bench-fit-decompose.rds")
cat("\nwrote -> dev/bench-fit-decompose.rds\n")
