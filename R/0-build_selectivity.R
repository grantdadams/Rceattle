#' Selectivity parameters that accept a linkage
#'
#' Slot names shared across the parametric forms, plus DoubleNormal aliases.
#' @keywords internal
#' @noRd
SEL_LINKAGE_PARAMS <- c("slp_asc", "slp_desc", "inf_asc", "inf_desc", "coff",
                        "sigma_asc", "sigma_desc", "peak", "right_floor")


#' @keywords internal
#' @noRd
.validate_sel_linkages <- function(linkages) {
  .validate_process_linkages(linkages, SEL_LINKAGE_PARAMS, "sel")
}


#' Selectivity specification
#'
#' @description
#' Carries environmental linkages on selectivity parameters. The effect on a
#' parameter is written as a formula and composes additively with any
#' `Time_varying_sel` process error on the same fleet (the two are separate
#' mechanisms: a covariate effect versus a deviation).
#'
#' The parameter names are the shape parameters of the parametric selectivity
#' forms:
#' \describe{
#'   \item{`slp_asc`, `slp_desc`}{ascending / descending logistic slope (log
#'     scale); for a double-normal the ascending / descending width, aliased
#'     `sigma_asc` / `sigma_desc`.}
#'   \item{`inf_asc`, `inf_desc`}{ascending / descending inflection age/length
#'     (natural scale); for a double-normal the peak and the logit right-floor,
#'     aliased `peak` / `right_floor`.}
#'   \item{`coff`}{non-parametric selectivity-at-bin coefficients.}
#' }
#'
#' Every parameter accepts `link = "log"` (multiplicative on the natural
#' parameter) or `link = "identity"` (additive), like the other processes.
#'
#' **Priors on a selectivity parameter.** An intercept-only formula (`~ 1`) with
#' a `priors` entry places a prior on the selectivity parameter itself (no
#' year-to-year offset is added). Read the prior on the parameter's own scale:
#' the slopes (`slp_asc` / `slp_desc`) are on the log scale (use `lognormal()`),
#' the inflections (`inf_asc` / `inf_desc`) on the natural scale (use `normal()`).
#' See Examples for a normal prior on the ascending inflection. This mirrors the
#' prior-only [build_composition()] path.
#'
#' A selectivity prior targets one parameter, so in a two-sex model an
#' unstratified `~ 1` prior constrains sex 1 only -- use `by = ~ sex` for a
#' per-sex prior. An `init` on a selectivity intercept has no effect (the
#' starting value comes from the data), and a prior on the double-normal
#' `right_floor` is not supported.
#' For a fleet that mirrors another fleet's selectivity (shared
#' `Selectivity_index`), place the prior on the lead fleet so the shared
#' parameter block is not penalized more than once.
#'
#' @param linkages Optional named list of [linkage_spec()] objects keyed by
#'   selectivity parameter. Coefficients are per fleet by default
#'   (`by = ~ fleet`); use the `fleet` argument of [linkage_spec()] to restrict
#'   a spec to particular fleets.
#'
#' @return A list of selectivity settings for [fit_mod()].
#'
#' @examples
#' \donttest{
#' # A cold-pool effect on the ascending inflection of a logistic fleet
#' build_selectivity(linkages = list(
#'   inf_asc = linkage_spec(~ cold_pool, by = ~ fleet)))
#'
#' # A normal prior on the ascending inflection (intercept-only formula)
#' build_selectivity(linkages = list(
#'   inf_asc = linkage_spec(~ 1, priors = list(`(Intercept)` = normal(0, 3)))))
#' }
#'
#' @export
build_selectivity <- function(linkages = NULL) {
  list(linkages = .validate_sel_linkages(linkages))
}


# Map a selectivity linkage param name to (array, 1-based slot). log_sel_slp
# and sel_inf are [2, fleet, sex]; sel_coff is [fleet, sex, bin].
.SEL_PARAM_TO_SLOT <- list(
  slp_asc     = list(arr = "log_sel_slp", slot = 1L),
  slp_desc    = list(arr = "log_sel_slp", slot = 2L),
  sigma_asc   = list(arr = "log_sel_slp", slot = 1L),
  sigma_desc  = list(arr = "log_sel_slp", slot = 2L),
  inf_asc     = list(arr = "sel_inf",     slot = 1L),
  inf_desc    = list(arr = "sel_inf",     slot = 2L),
  peak        = list(arr = "sel_inf",     slot = 1L),
  right_floor = list(arr = "sel_inf",     slot = 2L),
  coff        = list(arr = "sel_coff",    slot = NA_integer_)
)


# Selectivity forms whose consume site reads the linkage offset: every
# PARAMETRIC form. DoubleNormal reuses the slp/inf slots (peak/sigma/floor
# aliases) and LogisticPM's multiplicative deviates carry the offset inside
# their exp. The non-parametric forms are excluded on purpose -- see the
# `coff` note in .check_sel_linkage_support().
.SEL_LINKAGE_WIRED_FORMS <- c("Logistic", "DoubleLogistic", "DescendingLogistic",
                              "DoubleNormal", "LogisticPM")
.SEL_LINKAGE_WIRED_PARAMS <- c("slp_asc", "slp_desc", "inf_asc", "inf_desc",
                               "sigma_asc", "sigma_desc", "peak", "right_floor")


#' Reject selectivity linkages the template does not yet consume
#'
#' @param linkage_table pooled linkage table (may be NULL / empty).
#' @param fleet_control the fleet control table.
#' @return invisibly NULL; errors on an unsupported sel linkage.
#' @keywords internal
#' @noRd
.check_sel_linkage_support <- function(linkage_table, fleet_control) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(invisible())
  sel <- linkage_table[linkage_table$process == "sel", , drop = FALSE]
  if (nrow(sel) == 0L) return(invisible())

  bad_param <- setdiff(unique(sel$param), .SEL_LINKAGE_WIRED_PARAMS)
  if (length(bad_param) > 0) {
    extra <- if ("coff" %in% bad_param) paste0(
      "\n  `coff` (non-parametric selectivity) cannot carry a linkage: those ",
      "forms mean-centre their coefficients each year, so a per-year offset ",
      "applied across all bins cancels exactly. A meaningful effect would need ",
      "a per-bin covariate, which the formula grammar does not express.") else ""
    stop(sprintf(
      paste0("selectivity linkage parameter(s) not supported: %s.\n",
             "  Supported: %s.%s"),
      paste(bad_param, collapse = ", "),
      paste(.SEL_LINKAGE_WIRED_PARAMS, collapse = ", "), extra), call. = FALSE)
  }

  # Every fleet a sel row targets must use a wired selectivity form.
  flts <- unique(sel$fleet)
  flts <- flts[!is.na(flts)]
  if (length(flts) == 0L) flts <- seq_len(nrow(fleet_control))  # NA = all fleets
  forms <- as.character(fleet_control$Selectivity[flts])
  bad_flt <- flts[!forms %in% .SEL_LINKAGE_WIRED_FORMS]
  if (length(bad_flt) > 0) {
    stop(sprintf(
      paste0("selectivity linkage on fleet(s) %s whose form (%s) is not yet ",
             "wired for linkages.\n  Wired forms: %s."),
      paste(fleet_control$Fleet_name[bad_flt], collapse = ", "),
      paste(unique(as.character(fleet_control$Selectivity[bad_flt])),
            collapse = ", "),
      paste(.SEL_LINKAGE_WIRED_FORMS, collapse = ", ")), call. = FALSE)
  }

  # A PRIOR on a selectivity intercept re-targets the base parameter, whose scale
  # depends on the fleet's form and whose ownership depends on Selectivity_index.
  # Two cases the re-target cannot yet express correctly are rejected up front
  # (a covariate linkage on the same slot is fine -- only priors are affected).
  prior_rows <- sel[!is.na(sel$prior_family) & sel$prior_family != "none", ,
                    drop = FALSE]
  if (nrow(prior_rows) > 0L) {
    row_flt <- function(f) if (is.na(f)) 1L else as.integer(f)   # NA fleet = cell 1 (cpp default)

    # (a) DoubleNormal stores sel_inf(1) as logit(right_floor), so a natural-scale
    # prior on `inf_desc` / `right_floor` would be evaluated on the logit scale.
    # Reject until the logit transform is wired -- the ascending peak (`inf_asc`)
    # and the sigmas/slopes (log scale) are unaffected.
    dn <- prior_rows[prior_rows$param %in% c("inf_desc", "right_floor"), , drop = FALSE]
    dn_flt <- unique(vapply(dn$fleet, row_flt, integer(1)))
    dn_flt <- dn_flt[as.character(fleet_control$Selectivity[dn_flt]) == "DoubleNormal"]
    if (length(dn_flt) > 0L) {
      stop(sprintf(paste0(
        "prior on `inf_desc` / `right_floor` for DoubleNormal fleet(s) %s is not ",
        "supported: that slot holds logit(right_floor), so a natural-scale prior ",
        "would be applied on the logit scale. Prior the ascending peak / sigmas ",
        "instead."),
        paste(fleet_control$Fleet_name[dn_flt], collapse = ", ")), call. = FALSE)
    }

    # (b) Fleets that mirror another fleet's selectivity (Selectivity_index != own
    # Fleet_code) share one parameter block; a prior on the mirror double-counts
    # the block (cf. the shared-block penalty trap). Require the prior on the lead
    # fleet (Selectivity_index == Fleet_code).
    sidx <- fleet_control$Selectivity_index
    mir_flt <- unique(vapply(prior_rows$fleet, row_flt, integer(1)))
    mir_flt <- mir_flt[!is.na(sidx[mir_flt]) & sidx[mir_flt] != mir_flt]
    if (length(mir_flt) > 0L) {
      stop(sprintf(paste0(
        "selectivity prior on fleet(s) %s that mirror another fleet's ",
        "selectivity (Selectivity_index != Fleet_code): the shared block would be ",
        "penalized once per sharing fleet. Place the prior on the lead fleet ",
        "(the one whose Selectivity_index equals its Fleet_code)."),
        paste(fleet_control$Fleet_name[mir_flt], collapse = ", ")), call. = FALSE)
    }

    # (c) A prior on a limb the fleet's own curve never uses. Logistic reads only
    # the ascending slots, DescendingLogistic only the descending ones; the other
    # pair stays at its build default and never enters selectivity-at-age. The
    # prior would still be added to the objective -- a constant that shifts the
    # reported likelihood and moves with an unrelated default, while doing
    # nothing to the fit. Silently accepting it is how a reconciliation against
    # another model picks up an unexplained offset.
    used <- list(Logistic           = c("slp_asc", "inf_asc"),
                 DescendingLogistic = c("slp_desc", "inf_desc"))
    for (form in names(used)) {
      f_rows <- prior_rows[vapply(prior_rows$fleet, function(f)
        as.character(fleet_control$Selectivity[row_flt(f)]) == form,
        logical(1)), , drop = FALSE]
      unused <- f_rows[!f_rows$param %in% used[[form]], , drop = FALSE]
      if (nrow(unused) > 0L) {
        stop(sprintf(paste0(
          "selectivity prior on `%s` for %s fleet(s) %s: that %s curve does not ",
          "use those parameters, so the prior would add a constant to the ",
          "objective without affecting the fit. Prior %s instead, or drop the ",
          "fleet from this prior's fleet list."),
          paste(unique(unused$param), collapse = "`, `"), form,
          paste(unique(fleet_control$Fleet_name[
            vapply(unused$fleet, row_flt, integer(1))]), collapse = ", "),
          form, paste(used[[form]], collapse = " / ")), call. = FALSE)
      }
    }
  }
  invisible()
}
