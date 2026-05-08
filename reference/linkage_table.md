# Linkage table: long-format coefficient registry for Rceattle processes

A *linkage table* is the unified data structure used to describe every
estimated coefficient that connects a process (recruitment, natural
mortality, growth, ...) to a covariate or stratum. Each row of the table
corresponds to exactly one scalar coefficient (a `beta`) and identifies
(a) which process and parameter it modifies, (b) which species/sex/age
subset it applies to, (c) which column of the shared design matrix `X`
it multiplies, and (d) its prior, bounds, link, and estimation phase.

## Details

The motivation is that the historical Rceattle pattern of one wide
parameter array per process (e.g. `beta_rec_pars[nspp, n_env]`,
`M1_beta[nspp, nsex, n_env]`) does not generalize cleanly across
processes that need different stratifications. A long-format table grows
only with the linkages the user actually requests, regardless of the
underlying species/sex/age dimensionality.

This file defines the schema and basic helpers; user-facing
formula-driven constructors live in `R/0-build_linkage.R`.
