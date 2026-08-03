# Fast dev builds for Rceattle.
#
# pkgload::load_all() / devtools::load_all() compile the TMB model in pkgbuild's
# DEBUG mode, whose flags end in -O0 (unoptimized). That makes fit_mod() roughly
# 10x slower than a normal `R CMD INSTALL .`, which compiles at -O2 (see
# dev/PERF-findings.md for the measured 10x, bit-identical). Setting
# pkg.build_extra_flags = FALSE drops the debug flags so load_all() compiles the
# model at the same -O2 as a production install -- dev fits run at production
# speed.
#
# Trade-off: an -O2 build is not gdb/lldb-friendly (no line-level C++ debugging).
# When you need to step through src/TMB/*.cpp, start R with RCEATTLE_DEBUG_CPP=1
# (or run `options(pkg.build_extra_flags = TRUE)` then re-load) to restore the
# -O0 debuggable build.
if (!nzchar(Sys.getenv("RCEATTLE_DEBUG_CPP"))) {
  options(pkg.build_extra_flags = FALSE)
}
