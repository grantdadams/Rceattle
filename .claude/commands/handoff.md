---
description: Update inst/dev/SESSION_HANDOFF.md so the next session can resume without re-deriving
---

Keep `inst/dev/SESSION_HANDOFF.md` current: what is in flight, what is verified, what is
blocked, and where to pick up. CLAUDE.md is policy and does not change often; this is state and
changes every session. Keeping them apart is why CLAUDE.md can be trusted.

## Before running

Read the existing file first and **edit it** — do not regenerate it. Losing a blocker someone
recorded is worse than an out-of-date line.

## The shape (keep these sections)

```
# Session handoff
## Now            — the one thing in flight, and its branch
## Done & verified — what shipped, with the evidence: which gate, what it showed
## Known flags     — things that are true and surprising, with file:line
## Blocked         — what is waiting, and on whom or what
## Resume here     — the next concrete action
```

## Rules

- **"Verified" means a gate was run and its output read.** Name the gate and the number. "Should
  be fine" belongs under Known flags, not Done.
- **Absolute dates**, never "last week".
- **A number that matters gets written down.** A golden objective, an SSB delta, a max gradient.
  The next session cannot re-derive what you did not record.
- **`[[TODO]]`, `[[VERIFY]]` and `[[author]]` markers are Grant's calls.** Surface them; never
  resolve one yourself.
- Planning documents and scratch stay in the untracked `dev/`. `inst/dev/` is committed and
  ships, so it holds things that stay true.

## After running

Show the diff to the handoff file and nothing else.
