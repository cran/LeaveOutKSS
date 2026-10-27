#' Linear Projections of Estimated Firm Effects with Kline, Saggio, and Solvsten (KSS) Standard Errors
#'
#' Regresses transformed fixed effects on observables and reports both naive and
#' Kline, Saggio, and Solvsten (KSS)-corrected standard errors. This
#' corresponds to the "lincom" discussion in the original vignette on
#' regressing firm effects on observables.
#'
#' @param y Outcome vector used to estimate the original model.
#' @param X Design matrix used to estimate the fixed effects model.
#' @param Z Matrix of observables used in the linear projection.
#' @param Transform Matrix that maps model coefficients into the fixed effect of
#'   interest, typically firm effects.
#' @param sigma_i Observation-specific leave-out variance estimates.
#' @param labels Optional labels for the columns of `Z`.
#'
#' @return An object of class `"lincom_kss_result"` containing a results table
#'   with coefficient estimates, naive standard errors, KSS-corrected standard
#'   errors, and t statistics.
#'
#' @seealso [leave_out_KSS()], [kss_quadratic_form()]
#'
#' @export
lincom_KSS <- function(y, X, Z, Transform, sigma_i, labels = NULL) {
  if (missing(y) || missing(X) || missing(Z) || missing(Transform) || missing(sigma_i)) {
    stop("`y`, `X`, `Z`, `Transform`, and `sigma_i` must all be supplied.", call. = FALSE)
  }
  Z <- .kss_validate_row_aligned_input(Z, nrow(X), "Z")
  .kss_assert_same_length(y = y, sigma_i = sigma_i)
  if (nrow(Transform) != nrow(X)) {
    stop("`Transform` must have the same number of rows as `X`.", call. = FALSE)
  }
  .kss_stop_if_missing_or_nonfinite(y, "y")
  .kss_stop_if_missing_or_nonfinite(Z, "Z")
  .kss_stop_if_missing_or_nonfinite(sigma_i, "sigma_i")

   # SET DIMENSIONS
  n <- nrow(X)
  r <- ncol(Z) + 1

  # Add Constant
  Z <- cbind(rep(1, nrow(Z)), Z)

  # COMPUTE
  xx <- t(X) %*% X # Equivalent to X'X in 'MATLAB'
  xy <- t(X) %*% y # Equivalent to X'y in 'MATLAB'
  beta <- sanic::solve_cg(a = xx, b = as.vector(xy), type = "CG", iter = min(1000L, max(dim(xx))), tol = 1e-10, precond = 2, verbose = FALSE)
  wy <- Transform %*% beta
  zz <- t(Z) %*% Z # Equivalent to Z'Z in 'MATLAB'
  if (anyNA(wy) || any(!is.finite(wy))) {
    stop("The transformed fixed effects contain missing or non-finite values.", call. = FALSE)
  }
  numerator <- qr.solve(Z, wy)
  sigma_i <- sparseMatrix(i = 1:n, j = 1:n, x = as.vector(sigma_i),dims = c(n,n))
  error_vec <- as.vector((y - X %*% beta)^2)
  sigma_i_naive <- sparseMatrix(i = 1:n, j = 1:n, x = error_vec,dims = c(n,n))
  denominator <- numeric(r)
  denominator_naive <- numeric(r)

  for (q in 1:r) {
    v <- sparseMatrix(i = q, j = 1, x = 1, dims = c(r, 1))
    v <- qr.solve(zz, v)
    v <- Z %*% v
    v <- t(Transform) %*% v
    right <- sanic::solve_cg(a = xx, b = as.vector(v), type = "CG", iter = min(1000L, max(dim(xx))), tol = 1e-6, precond = 2, verbose = FALSE)
    left <- t(right)
    denominator[q] <- left %*% (t(X) %*% sigma_i %*% X) %*% right
    denominator_naive[q] <- left %*% (t(X) %*% sigma_i_naive %*% X %*% right)
  }

  kss_se <- rep(NaN, r)
  naive_se <- rep(NaN, r)
  test_statistic <- rep(NaN, r)

  if (r >= 2) {
    kss_se[2:r] <- sqrt(denominator[2:r])
    naive_se[2:r] <- sqrt(denominator_naive[2:r])
    test_statistic[2:r] <- numerator[2:r] / kss_se[2:r]
  }

  if (is.null(labels)) {
    terms <- paste("Column", seq_len(r - 1), "of Z")
  } else {
    if (length(labels) != (r - 1)) {
      stop("`labels` must have one entry for each column of `Z`.", call. = FALSE)
    }
    terms <- as.character(labels)
  }

  test_statistic <- test_statistic[-1]
  linear_combination <- numerator[-1]
  SE_linear_combination_KSS <- kss_se[-1]
  white_se <- naive_se[-1]

  results <- data.frame(
    term = terms,
    estimate = linear_combination,
    white_se = white_se,
    kss_se = SE_linear_combination_KSS,
    t_statistic = test_statistic,
    check.names = FALSE
  )

  structure(
    list(
      call = match.call(),
      results = results,
      test_statistic = test_statistic,
      linear_combination = linear_combination,
      SE_linear_combination_KSS = SE_linear_combination_KSS
    ),
    class = "lincom_kss_result"
  )
}
