#' Compare Two-Way Fixed Effects and Saturated-Model R-Squared Values
#'
#' Computes goodness-of-fit summaries for a two-way fixed effects model and for
#' a saturated worker-firm interaction model on the same sample. The function is
#' intended as a diagnostic companion to the leave-out decomposition routines and
#' follows the same basic data-preparation conventions.
#'
#' @param y Numeric outcome vector.
#' @param id Worker identifier vector.
#' @param firmid Firm identifier vector.
#' @param controls Optional matrix or vector of additional controls.
#' @param txt_file Optional path for exporting a text summary of the comparison.
#' @param progress Logical scalar indicating whether stage progress messages
#'   should be emitted.
#'
#' @return An object of class `"rsquared_comp_result"` containing a summary
#'   table for the two fitted models and the elapsed time.
#'
#' @details
#' The two-way fixed effects model includes worker effects, firm effects, and
#' optional controls. The saturated model replaces separate worker and firm
#' effects with worker-firm interaction indicators. Comparing the two summaries
#' can be useful when evaluating how much additional fit is obtained by moving
#' from the standard Abowd, Kramarz, and Margolis (1999; AKM) specification to
#' a fully saturated match design.
#'
#' @seealso [leave_out_KSS()], [leave_out_KSS_fe()], [fast_fe_est()]
#'
#' @examples
#' path <- system.file("extdata", "test.csv", package = "LeaveOutKSS")
#' dt <- data.table::fread(path, header = FALSE)
#'
#' res <- rsquared_comp(
#'   y = dt[[4]],
#'   id = dt[[1]],
#'   firmid = dt[[2]],
#'   progress = FALSE
#' )
#'
#' print(res)
#'
#' @export
rsquared_comp <- function(y, id, firmid,
                          controls = NULL,
                          txt_file = NULL,
                          progress = FALSE) {
  txt_file <- .kss_validate_optional_path(txt_file, "txt_file")
  no_controls <- 0
  start_time <- proc.time()[["elapsed"]]

  # Check for minimum required arguments
  if (missing(y) || missing(id) || missing(firmid)) {
    stop("`y`, `id`, and `firmid` must all be supplied.", call. = FALSE)
  }
  .kss_validate_flag(progress, "progress")
  .kss_assert_same_length(y = y, id = id, firmid = firmid)
  .kss_stop_if_missing_or_nonfinite(y, "y")
  .kss_stop_if_missing_or_nonfinite(id, "id")
  .kss_stop_if_missing_or_nonfinite(firmid, "firmid")

  controls <- .kss_validate_row_aligned_input(controls, length(y), "controls")
  if (!is.null(controls)) {
    .kss_stop_if_missing_or_nonfinite(controls, "controls")
  }
  
  if (is.null(controls)) {
    controls <- matrix(1, nrow=length(y), ncol=1)
    no_controls <- 1
  } else {
    controls <- cbind(1, controls)
  }

  .kss_progress(progress, "rsquared_comp: preparing fixed-effects comparison")
  # Assume that y, id, firmid, and controls are available and controls is a matrix
  DT <- data.table(y = y, id = id, old_id = id, firmid = firmid, old_firmid = firmid)
  DT_controls <- as.data.table(controls)

  # Drop stayers with a single person year observation
  DT[, count := .N, by = id]
  DT_controls <- DT_controls[DT$count > 1]
  DT <- DT[count > 1]
  .kss_stop_if_no_rows(nrow(DT), "rsquared_comp")

  # Resetting ids
  DT[, firmid := match(firmid, unique(firmid))]
  DT[, id := match(id, unique(id))]
  DT[, interact := paste(id, firmid, sep = "-")]
  DT[, int_id := match(interact, unique(interact))]

  NT <- nrow(DT)
  D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1)
  N <- ncol(D)
  F <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1)
  J <- ncol(F)

  .kss_progress(progress, "rsquared_comp: fitting two-way fixed-effects benchmark")
  S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x = 0, dims = c(1, J-1)))
  X <- cbind2(cbind2(D, F %*% S), as.matrix(DT_controls))
  K <- ncol(X)
  K_twfe <- K

  # Calculating xx (X' * X) and xy (X' * y)
  xx <- t(X) %*% X
  xy <- t(X) %*% DT$y

  # Solving for beta
  b <- sanic::solve_cg(a = xx, b = xy, type = "CG", iter = min(1000L, max(dim(xx))), tol = 1e-10, precond = 2, verbose = FALSE)
  b <- as.vector(b)

  DT[, y_hat := as.vector(X %*% b)]
  mean_y <- mean(DT$y)
  DT[, residuals := (y - y_hat)^2]
  DT[, total := (y - mean_y)^2]
  r_squared_twfe <- 1 - sum(DT$residuals) / sum(DT$total)
  adj_r_squared_twfe <- 1 - (1 - r_squared_twfe) * (NT - 1) / (NT - K - 1)

  .kss_progress(progress, "rsquared_comp: fitting saturated interaction benchmark")
  flag <- DT$int_id != 1
  D <- sparseMatrix(
    i = which(DT$int_id != 1),
    j = as.integer(factor(DT[flag, ]$int_id)),
    x = 1,
    dims = c(nrow(DT), length(unique(DT$int_id)) - 1)
  )
  # Constructing the X matrix
  X <- cbind2(D, as.matrix(DT_controls))
  K <- ncol(X)

  # Calculating xx (X' * X) and xy (X' * y)
  xx <- t(X) %*% X
  xy <- t(X) %*% DT$y

  # Solving for beta
  b <- sanic::solve_cg(a = xx, b = xy, type = "CG", iter = min(1000L, max(dim(xx))), tol = 1e-10, precond = 2, verbose = FALSE)
  b <- as.vector(b)

  DT[, y_hat := as.vector(X %*% b)]
  mean_y <- mean(DT$y)
  DT[, residuals := (y - y_hat)^2]
  DT[, total := (y - mean_y)^2]
  r_squared_saturated <- 1 - sum(DT$residuals) / sum(DT$total)
  adj_r_squared_saturated <- 1 - (1 - r_squared_saturated) * (NT - 1) / (NT - K - 1)

  summary_tbl <- data.frame(
    model = c("TWFE", "Saturated"),
    r_squared = c(r_squared_twfe, r_squared_saturated),
    adjusted_r_squared = c(adj_r_squared_twfe, adj_r_squared_saturated),
    n_parameters = c(K_twfe, K),
    check.names = FALSE
  )

  elapsed_time <- proc.time()[["elapsed"]] - start_time
  result <- structure(
    list(
      call = match.call(),
      controls_included = !no_controls,
      summary = summary_tbl,
      elapsed_seconds = elapsed_time
    ),
    class = "rsquared_comp_result"
  )

  if (!is.null(txt_file)) {
    .kss_progress(progress, sprintf("rsquared_comp: writing summary to %s", txt_file))
    .kss_write_text(.kss_format_rsquared_lines(result), txt_file)
  }

  result
}
