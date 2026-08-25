---
description: Recompile the TMB C++ model and reload Rceattle
---

Recompile the TMB model and reload the package, then report the result cleanly.

Run from the repo root:

```
export PATH=/usr/bin:$PATH && Rscript -e 'pkgload::load_all(".", quiet = TRUE)'
```

Notes:
- The `export PATH=/usr/bin:$PATH` prefix puts the system toolchain first — this repo's
  sessions do it consistently so the TMB DLL compiles cleanly (avoids a Homebrew
  clang/gfortran shadowing the build). Keep it on every compile/check command.
- `load_all()` recompiles `src/TMB/*.cpp` + `*.hpp` (slow, ~tens of seconds).
- If `$ARGUMENTS` contains `r-only` or `no-compile`, skip the rebuild by passing
  `compile = FALSE` to `load_all()` (use this when you only changed `R/`).
- If compilation fails, surface the compiler error with file:line and STOP — do not
  continue as if it built.
