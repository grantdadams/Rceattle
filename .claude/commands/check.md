---
description: Run R CMD check the way CI does (backgrounded, slow)
---

Run a package check similar to CI. This is slow — run it in the **background** and report
when it finishes.

```
export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), error_on = "warning")'
```

- If `$ARGUMENTS` includes `no-vignettes`, also add `"--no-build-vignettes"` to `args`
  (this repo's sessions often skip vignettes for speed).
- When it completes, summarize ERRORs / WARNINGs / NOTEs. Don't try to fix them unless I
  ask.
