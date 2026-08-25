# `.claude/` — shared, and what that means

`settings.json`, `commands/` and `agents/` are **committed**. They arrive by `git pull`, so
treat a change to them as you would a change to a CI workflow: it takes effect on everyone's
machine without anyone approving it locally.

`settings.local.json` is gitignored. **Anything personal or machine-specific belongs there** —
extra `WebFetch` domains, broader `Bash`/`Edit` permissions, model preferences.

`.Rbuildignore` keeps this whole directory out of the package tarball, so nothing here ships.

## Why `settings.json` grants so little

It allows exactly eight documentation domains: the TMB, WHAM, FIMS, sdmTMB, PBS-assess and
glmmTMB sites, plus R's own docs and CRAN. These are the references `CLAUDE.md` tells you to
consult when documenting a switch or naming a process argument.

It deliberately does **not** pre-approve `github.com`, `raw.githubusercontent.com` or the
publisher sites. Those are not "documentation for a handful of model families" — they are most
of the internet, and a committed allow-list is the wrong place to widen reach for everyone. Add
them to your own `settings.local.json` if you want them.

There is no `deny` block. The rules that matter here are enforced by review and by the
verification gates in `CLAUDE.md`, not by permissions.

## The commands

Each wraps an operation that spans several files and is easy to leave half-done. They state
their preconditions, and — more importantly — what they do **not** cover:

| Command | Covers |
|---|---|
| `/recompile`, `/test`, `/document`, `/check` | the dev loop |
| `/golden-check` | four-model numeric equivalence |
| `/verify` | the `tools/verify/` harnesses golden-check cannot reach |
| `/new-column` | a schema column, end to end |
| `/doc-sync` | NEWS / DESCRIPTION / vignettes / pkgdown |
| `/pkgdown-check` | the pkgdown build a PR to `dev` otherwise never gets |
| `/ecosystem-sweep` | the sibling assessment repos |
| `/handoff` | `inst/dev/SESSION_HANDOFF.md` |

`agents/density-reviewer.md` is an adversarial statistical reviewer for likelihood and
AD-taping changes. It is read-only by construction.
