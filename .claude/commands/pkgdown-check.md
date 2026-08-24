---
description: Build the pkgdown reference index — the only pkgdown check a PR to dev gets
---

`pkgdown.yaml` runs on push to `main` and on PRs targeting `main` — **not `dev`**. So a
pkgdown break on a dev-targeted PR is invisible until the release PR opens. This is that check,
and it takes seconds.

## Before running

Run it after any change that adds, renames, or removes a documented topic — which includes
adding `@export`, adding `@rdname`, or changing `@keywords internal`.

## Invocation

```
export PATH=/usr/bin:$PATH && Rscript -e 'pkgdown::build_reference_index(pkgdown::as_pkgdown("."))'
```

## After running

Two failure modes, both silent in CI on this branch:

- **A documented topic missing from `_pkgdown.yml` stops the build.** `simulate.Rceattle` did
  exactly this from 5.9.0. Add it to the right section.
- **`@keywords internal` drops a topic from the site silently**, which is wrong for one the
  help or a vignette tells users to read.

Note that a function documented under another topic via `@rdname` needs no entry of its own —
`rearrange_dat()` shares `rearrange_data`'s topic, so it is not missing.

Report the exit status and any topic it named. Don't add entries to `_pkgdown.yml` for topics
that should not be on the site — make them internal instead, and say which you chose.
