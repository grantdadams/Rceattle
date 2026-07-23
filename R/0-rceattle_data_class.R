# =============================================================================
# Rceattle_data S3 class + spec-tree printer
# =============================================================================
#
# build_data() tags its result with the class "Rceattle_data" so that printing a
# data object shows an indented spec tree (dimensions -> fleets -> per-process
# form -> active linkages -> configuration) instead of dumping the ~40-element
# list. The same renderer feeds print.Rceattle() for a fitted model. The class
# is a thin tag: the object is still a plain list to every function that consumes
# it (clean_data(), data_check(), fit_mod(), write_data() all treat it as a
# list), so tagging cannot change a fit.

# Readable labels for the small fixed set of integer-coded fleet columns, used
# only when switch_check() could not normalise a (partial) object to strings.
.RCE_FLEET_TYPE_LABEL <- c("0" = "off", "1" = "fishery", "2" = "survey")

#' Render an Rceattle data list as an indented spec tree
#'
#' @param dl A data list (from build_data()/read_data(), or a fitted model's
#'   `$data_list`).
#' @return A character vector, one element per line (invisibly printable).
#' @keywords internal
#' @noRd
.rce_spec_tree <- function(dl) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # Normalise switches to their canonical string forms so the tree is readable
  # whether it is fed a raw (integer-coded) build_data() object or an already
  # switch-checked fitted data_list. Fall back to the raw object if a partial
  # object cannot be normalised.
  dl2 <- tryCatch(
    suppressWarnings(suppressMessages(Rceattle::switch_check(dl))),
    error = function(e) dl
  )

  L <- character(0)
  add <- function(...) L[[length(L) + 1L]] <<- paste0(...)

  # ---- header + dimensions --------------------------------------------------
  spn <- dl2$spnames %||% (if (!is.null(dl2$nspp)) seq_len(dl2$nspp) else NULL)
  nspp <- dl2$nspp %||% length(spn)
  add("  dimensions")
  add("  ├─ species  : ", nspp %||% "?",
      if (!is.null(spn)) paste0("  (", paste(spn, collapse = ", "), ")") else "")
  yr <- paste0(dl2$styr %||% "?", "–", dl2$endyr %||% "?")
  if (!is.null(dl2$projyr) && !identical(dl2$projyr, dl2$endyr)) {
    yr <- paste0(yr, "  (projection → ", dl2$projyr, ")")
  }
  add("  ├─ years    : ", yr)
  agr <- if (!is.null(dl2$minage) && !is.null(dl2$nages)) {
    paste0(dl2$minage, "–", dl2$minage + dl2$nages - 1)
  } else as.character(dl2$nages %||% "?")
  add("  └─ ages     : ", paste(agr, collapse = ", "),
      "  |  lengths ", paste(dl2$nlengths %||% "?", collapse = ", "),
      "  |  sexes ", paste(dl2$nsex %||% "?", collapse = ", "))

  # ---- fleets ---------------------------------------------------------------
  # Coerce to a data.frame: a fitted model$data_list holds the original
  # data.frame, but harden against a rearranged (matrix) fleet_control so a
  # printer can never crash auto-printing a model.
  fc <- tryCatch(as.data.frame(dl2$fleet_control, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (!is.null(fc) && nrow(fc) > 0) {
    add("  fleets (", nrow(fc), ")")
    lab <- function(col, i) {
      v <- if (col %in% colnames(fc)) fc[[col]][i] else NA
      if (is.na(v)) return(NA_character_)
      if (col == "Fleet_type" && v %in% names(.RCE_FLEET_TYPE_LABEL)) {
        return(unname(.RCE_FLEET_TYPE_LABEL[as.character(v)]))
      }
      as.character(v)
    }
    # Which Selectivity_index / Q_index values are shared (mirrored)?
    shared <- function(col) {
      if (!col %in% colnames(fc)) return(integer(0))
      t <- table(fc[[col]][!is.na(fc[[col]])])
      as.integer(names(t[t > 1]))
    }
    sel_mir <- shared("Selectivity_index"); q_mir <- shared("Q_index")
    for (i in seq_len(nrow(fc))) {
      last <- i == nrow(fc)
      elbow <- if (last) "  └─ " else "  ├─ "
      ftype <- lab("Fleet_type", i)
      sel   <- lab("Selectivity", i)
      qform <- lab("Catchability", i)
      parts <- c(
        paste0("[", fc$Fleet_code[i] %||% i, "] ", fc$Fleet_name[i]),
        if (!is.na(ftype)) ftype,
        if (!is.na(sel))   paste0("sel ", sel),
        if (!is.na(qform)) paste0("q ", qform)
      )
      # Mirror annotations.
      mir <- c()
      si <- if ("Selectivity_index" %in% colnames(fc)) fc$Selectivity_index[i] else NA
      qi <- if ("Q_index" %in% colnames(fc)) fc$Q_index[i] else NA
      if (!is.na(si) && si %in% sel_mir) mir <- c(mir, paste0("sel↔", si))
      if (!is.na(qi) && qi %in% q_mir)   mir <- c(mir, paste0("q↔", qi))
      if (length(mir)) parts <- c(parts, paste0("[", paste(mir, collapse = ", "), "]"))
      add(elbow, paste(parts, collapse = "   "))
    }
  }

  # ---- processes (present on a configured / fitted object) ------------------
  procs <- c()
  if (!is.null(dl2$srr_fun))      procs <- c(procs, paste0("recruitment  srr_fun = ", paste(dl2$srr_fun, collapse = ", ")))
  if (!is.null(dl2$M1_model))     procs <- c(procs, paste0("M1           M1_model = ", paste(dl2$M1_model, collapse = ", ")))
  if (!is.null(dl2$growth_model)) procs <- c(procs, paste0("growth       growth_model = ", paste(dl2$growth_model, collapse = ", ")))
  if (length(procs)) {
    add("  processes")
    for (k in seq_along(procs)) {
      add(if (k == length(procs)) "  └─ " else "  ├─ ", procs[k])
    }
  }

  # ---- active linkages ------------------------------------------------------
  link_slots <- c(q = "q_linkages", sel = "sel_linkages", growth = "growth_linkages",
                  M1 = "M1_linkages", srr = "srr_linkages", comp = "comp_linkages")
  link_lines <- c()
  for (nm in names(link_slots)) {
    lk <- dl2[[link_slots[[nm]]]]
    if (!is.null(lk) && length(lk) > 0) {
      forms <- vapply(seq_along(lk), function(j) {
        f <- lk[[j]]
        ff <- tryCatch(deparse(f$formula %||% attr(f, "formula")), error = function(e) "")
        paste0(names(lk)[j] %||% nm, " ", ff[1])
      }, character(1))
      link_lines <- c(link_lines, paste0(nm, ": ", forms))
    }
  }
  if (length(link_lines)) {
    add("  linkages (", length(link_lines), ")")
    for (k in seq_along(link_lines)) {
      add(if (k == length(link_lines)) "  └─ " else "  ├─ ", link_lines[k])
    }
  }

  # ---- model_config ---------------------------------------------------------
  cfg <- dl2$model_config
  if (!is.null(cfg)) {
    add("  model_config : msmMode ", cfg$msmMode %||% "?",
        " · initMode ", cfg$initMode %||% "?",
        if (!is.null(cfg$HCR$HCR)) paste0(" · HCR ", cfg$HCR$HCR) else "")
  }

  L
}

#' Print method for an Rceattle data list
#'
#' Shows the model specification as an indented tree -- dimensions, fleets and
#' their selectivity / catchability forms, configured processes, active
#' linkages, and any attached [model_config()] -- rather than dumping the full
#' data list.
#'
#' @param x An `"Rceattle_data"` object from [build_data()].
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @export
print.Rceattle_data <- function(x, ...) {
  cat("<Rceattle data>\n")
  tree <- tryCatch(.rce_spec_tree(x),
                   error = function(e) paste0("  (spec tree unavailable: ",
                                              conditionMessage(e), ")"))
  cat(paste(tree, collapse = "\n"), "\n", sep = "")
  invisible(x)
}

#' @rdname print.Rceattle_data
#' @param object An `"Rceattle_data"` object.
#' @export
summary.Rceattle_data <- function(object, ...) {
  print(object, ...)
  invisible(object)
}
