# Encode an Rceattle linkage table into TMB-friendly inputs

The TMB template consumes the pooled linkage table as a set of parallel
`IVECTOR` / `VECTOR` inputs plus a dense design matrix. R encodes
string-valued columns (`process`, `param`, `link`, `prior_family`) as
0-based integer codes (TMB-friendly); converts `NA` stratum ids
(`species`, `sex`, `age_bin`) to a sentinel `0` meaning "applies to all
levels" (TMB-side dispatch expands the shared rows over the relevant
1-based levels); converts the 1-based `X_col` to 0-based.

## Details

This is mechanism only; no model behavior depends on the encoding until
the TMB template wires up the corresponding inputs.
