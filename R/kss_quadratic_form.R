#' Evaluate Plug-In and Kline, Saggio, and Solvsten (KSS)-Corrected Quadratic Forms
#'
#' Computes a covariance-like quadratic form from transformed coefficient
#' vectors and subtracts the Kline, Saggio, and Solvsten (KSS) bias adjustment
#' based on observation-specific variances and `Bii` weights.
#'
#' @param sigma_i Vector of leave-out variance estimates.
#' @param A_1 Matrix used to transform the coefficient vector on the left side of
#'   the quadratic form.
#' @param A_2 Matrix used to transform the coefficient vector on the right side
#'   of the quadratic form.
#' @param beta Estimated coefficient vector.
#' @param Bii Vector of observation-specific bias terms for the target
#'   variance component.
#'
#' @return A named list with `theta`, the plug-in estimate, and `theta_KSS`, the
#'   bias-corrected estimate.
#'
#' @seealso [leave_out_KSS()], [lincom_KSS()]
#'
#' @examples
#' A <- diag(2)
#' kss_quadratic_form(
#'   sigma_i = c(1, 2),
#'   A_1 = A,
#'   A_2 = A,
#'   beta = c(0.5, 1),
#'   Bii = c(0.1, 0.2)
#' )
#'
#' @export
kss_quadratic_form <- function(sigma_i, A_1, A_2, beta, Bii) {
  # Calculate the right and left terms
  right <- A_2 %*% beta
  left <- A_1 %*% beta

  # Calculate means
  mean_left <- mean(left)
  mean_right <- mean(right)

  # Calculate the covariance manually to avoid memory issues
  theta <- sum((left - mean_left) * (right - mean_right)) / (length(left) - 1)

  # Calculate the degrees of freedom
  dof <- length(left) - 1

  # Calculate the KSS adjustment
  theta_KSS <- theta - (1 / dof) * sum(Bii * sigma_i)

  # Return the results as a list
  return(list(theta = theta, theta_KSS = theta_KSS))
}
