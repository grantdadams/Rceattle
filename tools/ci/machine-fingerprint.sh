#!/usr/bin/env bash
# What machine produced this number?
#
# Two jobs here pin numeric results and can go red for reasons that have nothing
# to do with the change under test: `golden` compares four objectives against
# literals, and `test-coverage` runs regression tests that pin fitted values.
# Both have done so. `inst/dev/TRAPS.md` records golden red on `main` for weeks
# over a summation order, and the 5.45.0 red was a gradient of 3.5e-4 on ubuntu
# against 2.5e-5 on macOS.
#
# The case that prompted this: a mirrored-block linkage fit returned THREE
# different objectives on ubuntu-latest (13807.6, 13818.2, 13817.0067) with the
# fitted coefficient ranging -0.0304 to -0.240043, on byte-identical model code,
# against 12879.997040 reproducibly on macOS arm64. Two fits in one process were
# bitwise identical, so the model is deterministic per machine and the variation
# is between machines -- but NONE of the data needed to say WHICH difference
# mattered was being recorded. The runner's own logs carry only
# "Platform: x86_64-pc-linux-gnu".
#
# So this prints the things that plausibly move a last bit, cheaply (~1 s), into
# every log that pins a number:
#
#   * CPU model and the FMA/AVX flags. glibc's libm dispatches exp/log/pow by
#     CPU feature through ifunc resolvers, and this model calls exp() 246 times
#     and log() 105 times per objective evaluation -- so a different kernel is a
#     different last bit, many times over. The runner fleet is heterogeneous.
#   * glibc, which is that libm.
#   * the compiler, which chooses whether to contract a*b+c into an FMA.
#   * R's BLAS/LAPACK, for completeness. Note this model does NOT route Eigen
#     through them (no EIGEN_USE_BLAS, and one BLAS-level op in the template),
#     so they are the least likely culprit here -- recorded so the next person
#     does not have to re-establish that, as happened once already.
#   * any threading override, since a multithreaded BLAS reduces in a
#     nondeterministic order.
#
# It is a REPORT, not a check: nothing here fails a job. The point is that a
# cross-platform disagreement becomes attributable instead of a mystery.
set -u

echo "=== machine fingerprint ==="
echo "  runner image : ${ImageOS:-unknown} ${ImageVersion:-unknown}"
echo "  arch         : $(uname -m)  $(uname -s)"

if command -v lscpu >/dev/null 2>&1; then
  echo "  cpu          : $(lscpu | sed -n 's/^Model name: *//p' | head -1)"
  echo "  cpu id       : family $(lscpu | sed -n 's/^CPU family: *//p' | head -1)"\
" model $(lscpu | sed -n 's/^Model: *//p' | head -1)"\
" stepping $(lscpu | sed -n 's/^Stepping: *//p' | head -1)"
elif [ "$(uname -s)" = "Darwin" ]; then
  echo "  cpu          : $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)"
fi

# The flags that change which libm kernel an ifunc resolver picks.
if [ -r /proc/cpuinfo ]; then
  flags=""
  for f in fma avx avx2 avx512f avx512dq; do
    if grep -qw "$f" /proc/cpuinfo 2>/dev/null; then flags="$flags $f"; fi
  done
  echo "  cpu features :${flags:- none detected}"
fi

if command -v ldd >/dev/null 2>&1; then
  echo "  glibc (libm) : $(ldd --version 2>/dev/null | head -1)"
fi
echo "  compiler     : $(${CC:-cc} --version 2>/dev/null | head -1 || echo unknown)"

# Threading overrides, if a caller has pinned them.
thr=""
for v in OPENBLAS_NUM_THREADS OMP_NUM_THREADS MKL_NUM_THREADS; do
  eval "val=\${$v:-}"
  if [ -n "${val}" ]; then thr="$thr $v=$val"; fi
done
echo "  thread limits:${thr:- unset (BLAS picks its own)}"

# R-side: version, and which BLAS/LAPACK it actually linked.
if command -v Rscript >/dev/null 2>&1; then
  Rscript -e 'si <- sessionInfo()
    cat("  R            :", R.version.string, "\n")
    cat("  R platform   :", R.version$platform, "\n")
    cat("  BLAS         :", if (is.null(si$BLAS)) "unknown" else si$BLAS, "\n")
    cat("  LAPACK       :", if (is.null(si$LAPACK)) "unknown" else si$LAPACK, "\n")
    for (p in c("TMB", "Matrix", "TMBhelper")) {
      v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
      cat(sprintf("  %-13s %s\n", p, v))
    }' 2>/dev/null || echo "  R            : present but sessionInfo failed"
fi
echo "=== end fingerprint ==="
