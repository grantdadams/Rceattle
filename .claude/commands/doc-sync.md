---
description: Check the current diff updated NEWS, DESCRIPTION, the vignettes and _pkgdown.yml
---

CLAUDE.md's rule: a behaviour, API, or doc change means **NEWS.md + `DESCRIPTION` `Version:` +
the affected vignette, in the same commit**. This reports whether that held for the current
diff. It reports; it does not fix.

## Invocation

```
git diff --name-only <merge-base>...HEAD
git diff <merge-base>...HEAD -- NEWS.md DESCRIPTION | head -40
```

## What to check

1. **Did `R/`, `src/` or `man/` change without a `NEWS.md` bullet?** Repo/tooling-only changes
   (`.claude/`, `tools/`, `.github/`) are exempt; say which exemption you applied.
2. **Is the `DESCRIPTION` `Version:` bumped, and at the right level?** Patch for bug fixes and
   docs, minor for new features, major for breaking. **"Breaking" here means no back-compat
   path** — a removal that ships with a deprecation message and keeps old fits working is a
   *minor* bump even when NEWS files it under `## Breaking changes`.
3. **Does NEWS open a section for a version `DESCRIPTION` never carries?** Intermediate headings
   written mid-branch get folded into the shipped section before merge. An unshipped section
   also breaks any `(x.y.z)` cross-reference pointing at it.
4. **Did an exported signature change without its vignette?** Vignettes are `eval = FALSE`, so
   nothing else will catch it.
5. **Did a documented topic appear or disappear?** Then `/pkgdown-check`.
6. **Does a NEWS bullet claim something the code does not do?** Read the bullet against the
   diff, not against the commit message.

## After running

A short table — NEWS / DESCRIPTION / vignettes / pkgdown, each ✓, ✗, or n/a with one line of
why — then only the problems. If everything is in order, say so in one sentence.
