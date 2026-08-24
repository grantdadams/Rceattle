---
description: Regenerate roxygen man/*.Rd + NAMESPACE
---

Regenerate package documentation after roxygen (`@…`) changes:

```
export PATH=/usr/bin:$PATH && Rscript -e 'suppressMessages(devtools::document(quiet = TRUE))'
```

Then show `git status --short man NAMESPACE` so I can see exactly what changed. Don't
commit. If `man/` or `NAMESPACE` changed unexpectedly (beyond what I edited), flag it.
