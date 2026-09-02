# Which model coordinate each estimated parameter belongs to.
#
# TMB names every element of `obj$par` after its parameter BLOCK, so a
# diagnostic naming the offender can only say `sel_coff_dev` -- 392 times, for
# one fleet's selectivity deviations. This turns that into the fleet, sex, bin
# and year the element actually is.

# Dimension token (parameter_dictionary()$dims) -> the axis it indexes. The
# dictionary is the source of truth for what an axis MEANS; the extent is read
# off the built array, because a declared extent can describe the estimated
# portion rather than the allocation (`init_dev` is [nspp, nages-1] in a
# [nspp, nages] array).
.PAR_AXIS <- c(
  nspp             = "species",
  nsex             = "sex",
  n_flt            = "fleet",
  n_sel            = "fleet",
  n_fsh            = "fleet",
  nyrs             = "year",
  nyrs_hind        = "year",
  nages            = "age",
  `nages-1`        = "age",
  n_sel_bins       = "bin",
  n_env            = "covariate",
  n_linkage        = "linkage",
  n_linkage_re     = "linkage",
  n_linkage_re_pen = "linkage",
  n_re_group       = "linkage",
  n_re_ar1_group   = "linkage",
  n_re_obs_group   = "linkage"
)

# Axes a bare integer dimension stands for, per block. A literal in the
# dictionary is a fixed-width slot whose meaning is specific to its parameter:
# the two limbs of a double-logistic, the three AR1 correlations. Without this
# the slot prints as a number, which is the problem being fixed.
.PAR_SLOT_LABELS <- list(
  log_sel_slp     = c("ascending", "descending"),
  sel_inf         = c("ascending", "descending"),
  log_sel_slp_dev = c("ascending", "descending"),
  sel_inf_dev     = c("ascending", "descending"),
  sel_curve_pen   = c("shape", "curvature", "dev magnitude"),
  M1_rho          = c("age", "year"),
  growth_log_sd   = c("mean", "plus group"),
  rec_pars        = c("R0/mean", "steepness", "beta"),
  weight_length_pars = c("a", "b")
)

#' Axis names for one parameter block
#'
#' Reads `parameter_dictionary()$dims` and returns one axis name per dimension,
#' in array order. A dimension the dictionary gives as a bare integer becomes
#' `"slot"`; one it does not name becomes `NA`.
#' @noRd
.rce_par_axes <- function(block, dict = parameter_dictionary()) {
  dims <- dict$dims[match(block, dict$internal)]
  if (is.na(dims)) return(character(0))
  toks <- trimws(strsplit(gsub("]", "", gsub("[", "", dims, fixed = TRUE),
                               fixed = TRUE), ",")[[1]])
  toks <- toks[nzchar(toks)]
  out <- unname(.PAR_AXIS[toks])
  out[is.na(out) & grepl("^[0-9]+$", toks)] <- "slot"
  out
}

#' Labels along one axis
#'
#' `n` comes from the built array, so an axis allocated wider than the
#' dictionary describes still labels every position it has.
#' @noRd
.rce_axis_labels <- function(axis, n, data_list, block = NULL, species = NA) {
  idx <- seq_len(n)
  pad <- function(x) if (length(x) >= n) as.character(x[idx]) else
    c(as.character(x), rep(NA_character_, n - length(x)))[idx]

  switch(axis,
    fleet = {
      nm <- data_list$fleet_control$Fleet_name
      if (is.null(nm)) as.character(idx) else pad(nm)
    },
    species = {
      nm <- data_list$spnames
      if (is.null(nm)) as.character(idx) else pad(nm)
    },
    # Ages run minage .. minage + nages - 1, so age `a` sits at index
    # a - minage + 1. minage is per species; without one the index is reported
    # rather than a guessed age.
    age = {
      ma <- data_list$minage
      if (is.null(ma) || is.na(species) || species > length(ma)) as.character(idx)
      else as.character(idx - 1L + ma[species])
    },
    # A selectivity bin is an index on the fleet's own Selectivity_dimension,
    # which may be length rather than age, so it is not converted to an age.
    bin  = as.character(idx),
    year = {
      if (is.null(data_list$styr)) as.character(idx)
      else as.character(data_list$styr + idx - 1L)
    },
    # Reported only where the model distinguishes sexes. Naming a sex on a
    # single-sex species implies a split that is not in the model, and the
    # array is allocated to the widest species, so the species decides.
    sex = {
      ns <- data_list$nsex
      single <- if (!is.na(species) && !is.null(ns) && species <= length(ns)) {
        ns[species] == 1L
      } else {
        n == 1L
      }
      if (single) rep(NA_character_, n) else c("female", "male")[idx]
    },
    slot = {
      lab <- .PAR_SLOT_LABELS[[block]]
      if (is.null(lab)) as.character(idx) else pad(lab)
    },
    as.character(idx)
  )
}

#' Locate every estimated parameter in the model's own coordinates
#'
#' @description
#' TMB names each element of `obj$par` after its parameter block, so a
#' convergence diagnostic can report that `sel_coff_dev` is non-identifiable but
#' not which fleet, sex, bin or year. This returns one row per estimated
#' (fixed-effect) parameter, carrying the coordinates it occupies.
#'
#' @details
#' The mapping is recovered by pushing a tagged parameter vector back through
#' TMB's own `parList()`, so it reflects `build_map()` exactly: parameters
#' mapped off are absent, and fleets sharing a `Selectivity_index` or
#' `Catchability_index` appear once, as the single parameter they are, with
#' `n_cells` recording how many array cells they drive.
#'
#' Coordinate columns are `NA` where the axis does not apply to a block.
#'
#' @param object A fitted Rceattle model from [fit_mod()].
#'
#' @return A data frame with one row per estimated parameter and columns
#'   `par_index` (position in `obj$par`), `block`, `species`, `fleet`, `sex`,
#'   `age`, `bin`, `year`, `slot`, `n_cells` and `label`.
#'
#' @seealso [parameter_dictionary()] for what each block means.
#' @export
parameter_index <- function(object) {
  obj <- object$obj
  if (is.null(obj) || is.null(obj$env$parList)) {
    stop("`object` does not carry a TMB object; refit with fit_mod().",
         call. = FALSE)
  }
  npar <- length(obj$par)
  if (npar == 0) return(.rce_par_index_empty())

  # A sentinel offset, so a tag cannot be confused with a parameter value left
  # in place for a mapped-off or random-effect cell.
  off <- 1e6
  pl  <- obj$env$parList(x = seq_len(npar) + off)
  dl  <- object$data_list
  dict <- parameter_dictionary()

  AX <- c("species", "fleet", "sex", "age", "bin", "year", "slot")
  rows <- list()

  for (block in names(pl)) {
    arr <- pl[[block]]
    v   <- as.vector(arr)
    hit <- which(v > off / 2)
    if (length(hit) == 0) next

    d <- dim(arr); if (is.null(d)) d <- length(arr)
    axes <- .rce_par_axes(block, dict)
    if (length(axes) != length(d)) axes <- rep(NA_character_, length(d))

    coord <- arrayInd(hit, .dim = d)          # cell coordinates, one row per hit
    par_i <- v[hit] - off

    # Species is resolved first, because minage and nsex are per species. A
    # block indexed by fleet rather than species still has one: fleet_control
    # says which species each fleet fishes or surveys.
    sp_col  <- match("species", axes)
    flt_col <- match("fleet", axes)
    sp_idx <- if (!is.na(sp_col)) {
      coord[, sp_col]
    } else if (!is.na(flt_col) && !is.null(dl$fleet_control$Species)) {
      as.integer(dl$fleet_control$Species)[coord[, flt_col]]
    } else {
      rep(NA_integer_, length(hit))
    }

    lab <- matrix(NA_character_, nrow = length(hit), ncol = length(AX),
                  dimnames = list(NULL, AX))
    # Axes whose labels depend on the species holding the row.
    per_sp <- c("age", "sex")
    for (k in seq_along(axes)) {
      ax <- axes[k]
      if (is.na(ax) || !ax %in% AX) next
      if (ax %in% per_sp) {
        lab[, ax] <- vapply(seq_along(hit), function(r)
          .rce_axis_labels(ax, d[k], dl, block, sp_idx[r])[coord[r, k]],
          character(1))
      } else {
        lab[, ax] <- .rce_axis_labels(ax, d[k], dl, block)[coord[, k]]
      }
    }
    # Name the species a fleet-indexed block belongs to, so a multispecies fit
    # says which stock a flagged parameter is on.
    if (is.na(sp_col) && !all(is.na(sp_idx))) {
      lab[, "species"] <- .rce_axis_labels("species", length(dl$spnames %||% 0),
                                           dl, block)[sp_idx]
    }

    df <- data.frame(par_index = par_i, block = block,
                     as.data.frame(lab, stringsAsFactors = FALSE),
                     stringsAsFactors = FALSE)
    rows[[block]] <- df
  }

  if (length(rows) == 0) return(.rce_par_index_empty())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # One row per PARAMETER. A mirrored parameter drives several cells; collapsing
  # keeps the count of unidentified things equal to the count of parameters, and
  # names every fleet that shares the block.
  out <- .rce_collapse_mirrored(out, AX)
  out <- out[order(out$par_index), , drop = FALSE]
  rownames(out) <- NULL
  out$label <- .rce_par_label(out, .rce_varying_axes(out, AX))
  out
}

# Axes worth printing: one that takes a single value across the whole model
# distinguishes nothing, and naming it on every line buries the axes that do.
# A single-species fit does not need "Pollock" on all 525 rows; a two-species
# fit does. The structured columns keep the value either way -- this governs
# only what the rendered label and the summary say.
#' @noRd
.rce_varying_axes <- function(df, AX) {
  keep <- vapply(AX, function(a)
    length(unique(stats::na.omit(df[[a]]))) > 1L, logical(1))
  AX[keep]
}

#' @noRd
.rce_par_index_empty <- function() {
  data.frame(par_index = integer(0), block = character(0), species = character(0),
             fleet = character(0), sex = character(0), age = character(0),
             bin = character(0), year = character(0), slot = character(0),
             n_cells = integer(0), label = character(0), stringsAsFactors = FALSE)
}

# Collapse the cells a mirrored parameter drives into one row, listing the
# distinct values on each axis.
#' @noRd
.rce_collapse_mirrored <- function(df, AX) {
  sp <- split(seq_len(nrow(df)), df$par_index)
  if (all(lengths(sp) == 1L)) { df$n_cells <- 1L; return(df) }
  one <- function(i) {
    r <- df[i, , drop = FALSE]
    o <- r[1, , drop = FALSE]
    for (ax in AX) {
      u <- unique(stats::na.omit(r[[ax]]))
      o[[ax]] <- if (length(u) == 0) NA_character_ else paste(u, collapse = " + ")
    }
    o$n_cells <- length(i)
    o
  }
  out <- do.call(rbind, lapply(sp, one))
  rownames(out) <- NULL
  out
}

# Axes that read as an ordered sequence collapse to a range; the rest name every
# distinct value, because "fleets 2-5" would imply an order fleets do not have.
.PAR_ORDINAL <- c("age", "bin", "year")

#' Summarise a set of flagged parameters by coordinate
#'
#' One line per block and per distinct combination of its categorical axes, with
#' the ordinal axes collapsed to ranges. 49 rows of `sel_coff_dev` become the
#' fleet, sex, bins and years they occupy.
#' @noRd
.rce_par_summary <- function(par_idx, index, max_lines = 8L) {
  if (is.null(index) || nrow(index) == 0 || length(par_idx) == 0) {
    return(character(0))
  }
  df <- index[match(par_idx, index$par_index), , drop = FALSE]
  df <- df[!is.na(df$block), , drop = FALSE]
  if (nrow(df) == 0) return(character(0))

  # Axes that vary across the MODEL, not across the flagged subset: an axis
  # constant everywhere distinguishes nothing wherever it is printed.
  AX  <- .rce_varying_axes(index, c("species", "fleet", "sex", "age", "bin",
                                    "year", "slot"))
  if (length(AX) == 0) AX <- "fleet"
  cat_ax <- setdiff(AX, .PAR_ORDINAL)
  key <- do.call(paste, c(list(df$block), lapply(cat_ax, function(a) df[[a]]),
                          list(sep = "\r")))

  rng <- function(v) {
    v <- stats::na.omit(v)
    if (length(v) == 0) return(NA_character_)
    u <- unique(v)
    n <- suppressWarnings(as.numeric(u))
    if (anyNA(n) || length(u) == 1L) return(paste(u, collapse = ", "))
    if (length(u) == 2L) return(paste(sort(u), collapse = ", "))
    paste0(min(n), "-", max(n))
  }

  lines <- vapply(split(seq_len(nrow(df)), key), function(i) {
    r <- df[i, , drop = FALSE]
    parts <- character(0)
    for (a in AX) {
      v <- if (a %in% .PAR_ORDINAL) rng(r[[a]]) else {
        u <- unique(stats::na.omit(r[[a]])); if (length(u) == 0) NA_character_ else paste(u, collapse = ", ")
      }
      if (is.na(v) || !nzchar(v)) next
      pre <- switch(a, age = if (grepl("-|,", v)) "ages " else "age ",
                       bin = if (grepl("-|,", v)) "bins " else "bin ", "")
      parts <- c(parts, paste0(pre, v))
    }
    # A block indexed by something with no model coordinate (a linkage group,
    # say) still has to say which elements it means.
    if (length(parts) == 0) {
      rg <- range(r$par_index)
      parts <- paste0("element ",
                      if (rg[1] == rg[2]) rg[1] else paste(rg, collapse = "-"))
    }
    sprintf("  %-16s %s  (%d)", r$block[1], paste(parts, collapse = ", "), nrow(r))
  }, character(1), USE.NAMES = FALSE)

  lines <- lines[order(-as.numeric(sub(".*\\((\\d+)\\)$", "\\1", lines)))]
  if (length(lines) > max_lines) {
    lines <- c(lines[seq_len(max_lines)],
               sprintf("  ... and %d more group(s)", length(lines) - max_lines))
  }
  lines
}

# The convergence battery must never fail because a label could not be built, so
# an index that cannot be recovered yields no coordinates rather than an error.
#' @noRd
.conv_par_index <- function(object) {
  tryCatch(parameter_index(object), error = function(e) NULL)
}

# Append the coordinate summary under a check's own message.
#' @noRd
.conv_with_coords <- function(msg, par_idx, index) {
  lines <- .rce_par_summary(par_idx, index)
  if (length(lines) == 0) return(msg)
  paste0(msg, "\n", paste(lines, collapse = "\n"))
}

# Put the coordinates beside the numbers in a check's own table, so
# `fit$convergence` carries them for a caller that wants to filter.
#' @noRd
.conv_attach_label <- function(tab, par_idx, index) {
  if (is.null(index) || nrow(index) == 0) return(tab)
  m <- match(par_idx, index$par_index)
  tab$where <- index$label[m]
  tab
}

# "GOA_pollock_fishery, female, age 3, 1985" -- the axes that apply, in reading
# order, so a diagnostic can print one string.
#' @noRd
.rce_par_label <- function(df, AX) {
  pre <- c(species = "", fleet = "", sex = "", age = "age ", bin = "bin ",
           year = "", slot = "")
  vapply(seq_len(nrow(df)), function(i) {
    parts <- character(0)
    for (ax in AX) {
      v <- df[[ax]][i]
      if (!is.na(v) && nzchar(v)) parts <- c(parts, paste0(pre[[ax]], v))
    }
    if (length(parts) == 0) "" else paste(parts, collapse = ", ")
  }, character(1))
}
