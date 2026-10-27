#' Fit a One-Way or Two-Way Fixed Effects Model
#'
#' Solves a fixed effects model using conjugate gradients and returns fitted
#' values and adjusted outcomes as an object. When `firmid` is omitted, the
#' routine estimates a one-way worker fixed effects model. When `firmid` is
#' supplied, it estimates a two-way worker-firm fixed effects model.
#'
#' @param y Numeric outcome vector.
#' @param id Worker identifier vector.
#' @param firmid Optional firm identifier vector. If `NULL`, a one-way model is
#'   fitted.
#' @param controls Optional matrix or vector of controls.
#' @param csv_file Optional path for exporting the fitted values table as a
#'   `.csv` file.
#' @param progress Logical scalar indicating whether stage progress messages
#'   should be emitted.
#'
#' @return An object of class `"fast_fe_est_result"` containing the fitted
#'   values table, model metadata, and elapsed time.
#'
#' @details
#' This helper is useful when the goal is to recover fitted values and
#' residualized outcomes rather than the leave-out variance decomposition. The
#' returned fitted-values table includes `y_hat`, `y_adj`, and the original
#' identifiers. When `csv_file` is supplied, that table is also written to disk.
#'
#' @seealso [leave_out_KSS()], [rsquared_comp()]
#'
#' @examples
#' path <- system.file("extdata", "test.csv", package = "LeaveOutKSS")
#' dt <- data.table::fread(path, header = FALSE)
#'
#' res <- fast_fe_est(
#'   y = dt[[4]],
#'   id = dt[[1]],
#'   firmid = dt[[2]],
#'   controls = cbind(year = dt[[3]])
#' )
#'
#' print(res)
#'
#' @export
fast_fe_est <- function(y, id, firmid = NULL,
                        controls = NULL,
                        csv_file = NULL,
                        progress = FALSE) {
  csv_file <- .kss_validate_optional_path(csv_file, "csv_file")
  no_controls <- 0
  start_time <- proc.time()[["elapsed"]]

  # Check for minimum required arguments
  if (missing(y) || missing(id)) {
    stop("`y` and `id` must both be supplied.", call. = FALSE)
  }
  .kss_validate_flag(progress, "progress")
  .kss_assert_same_length(y = y, id = id)
  .kss_stop_if_missing_or_nonfinite(y, "y")
  .kss_stop_if_missing_or_nonfinite(id, "id")
  if (!is.null(firmid)) {
    .kss_assert_same_length(y = y, id = id, firmid = firmid)
    .kss_stop_if_missing_or_nonfinite(firmid, "firmid")
  }
  model_type <- if (is.null(firmid)) "one_way" else "two_way"
  .kss_progress(progress, sprintf("fast_fe_est: fitting a %s fixed-effects model", gsub("_", "-", model_type)))

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
  
  .kss_progress(progress, "fast_fe_est: preparing sparse design matrices")
  if (is.null(firmid)) {
    DT <- data.table(y = y, id = id, old_id=id)
  } else {
    DT <- data.table(y = y, id = id,old_id = id, firmid = firmid, old_firmid=firmid)
  }
  DT_controls <- as.data.table(controls)
  

  
  
  # Drop stayers with a single person year observation
  DT[, count := .N, by = id]
  DT_controls <- DT_controls[DT$count > 1]
  DT <- DT[count > 1]
  .kss_stop_if_no_rows(nrow(DT), "fast_fe_est")
  
  # Resetting ids
  if (!is.null(firmid)){
    DT[, firmid := match(firmid, unique(firmid))]
  }
  DT[, id := match(id, unique(id))]
  
  NT <- nrow(DT)
  
  D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1)
  N <- ncol(D)
  if (!is.null(firmid)){
    F <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1)
    J <- ncol(F)
  }
 
  # STEP 3: Residualizing
  .kss_progress(progress, "fast_fe_est: residualizing outcome")
  
  
  if (!is.null(firmid)) {
    # Creating the restriction matrix S
    S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(1, J-1)))
    
    # Constructing the X matrix
    X <- cbind2(cbind2(D, F %*% S), as.matrix(DT_controls))
    n_fe_cols <- N + J - 1L
  } else {
    X <- cbind2(D, as.matrix(DT_controls))
    n_fe_cols <- N
  }
    # Calculating xx (X' * X) and xy (X' * y)
    xx <- t(X) %*% X
    xy <- t(X) %*% DT$y
    
    
    
    # Solving for beta
    b <- sanic::solve_cg(a = xx, b = xy, type = "CG", iter = min(1000L, max(dim(xx))), tol = 1e-10, precond = 2, verbose = FALSE)
    b <- as.vector(b)
  control_cols <- seq.int(n_fe_cols + 1L, ncol(X))
  DT[, y_adj := y - as.vector(X[, control_cols, drop = FALSE] %*% b[control_cols])]
  DT[, y_hat := as.vector(X %*% b)]

  if (!is.null(firmid)){
    output_dt <- data.table(
      y_hat = DT$y_hat,
      y_adj = DT$y_adj,
      personid = DT$old_id,
      firm_id = DT$old_firmid
    )
  } else {
    output_dt <- data.table(
      y_hat = DT$y_hat,
      y_adj = DT$y_adj,
      personid = DT$old_id)
  }

  elapsed_time <- proc.time()[["elapsed"]] - start_time
  result <- structure(
    list(
      call = match.call(),
      model = model_type,
      controls_included = !no_controls,
      fitted_values = output_dt,
      coefficients = b,
      elapsed_seconds = elapsed_time
    ),
    class = "fast_fe_est_result"
  )

  if (!is.null(csv_file)) {
    .kss_progress(progress, sprintf("fast_fe_est: writing fitted values to %s", csv_file))
    .kss_write_csv(output_dt, csv_file)
  }

  result
}
