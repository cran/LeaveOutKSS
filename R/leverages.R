#' Compute Statistical Leverages and Bias Terms
#'
#' Computes the observation-level leverage quantities used in the Kline,
#' Saggio, and Solvsten (KSS) bias correction, either exactly or with a
#' Johnson-Lindenstrauss approximation (JLA).
#'
#' @param X_fe Matrix used for the firm-effect variance component.
#' @param X_pe Matrix used for the person-effect variance component.
#' @param X Main design matrix.
#' @param xx Crossproduct matrix `t(X) %*% X`.
#' @param type_algorithm Character scalar, either `"exact"` or `"JLA"`.
#' @param scale Number of random projections when `type_algorithm = "JLA"`.
#' @param progress Logical scalar indicating whether leverage progress should be
#'   displayed.
#'
#' @return A list with elements `Pii`, `Mii`, `correction_JLA`, `Bii_fe`,
#'   `Bii_cov`, and `Bii_pe`.
#'
#' @details
#' The exact branch solves one linear system per observation. The
#' Johnson-Lindenstrauss approximation (JLA) branch follows the randomized
#' projection logic described in the original vignette to approximate the same
#' quantities at lower computational cost on large panels.
#'
#' @seealso [leverages_parallel()], [leave_out_KSS()]
#'
#' @export
leverages <- function(X_fe, X_pe, X, xx, type_algorithm, scale,
                      progress = FALSE) {
  .kss_validate_flag(progress, "progress")
  .kss_validate_scalar_number(scale, "scale", integerish = TRUE, min_value = 1)
  NT <- nrow(X)
  tolProb=0.5;
  Pii <- rep(0, NT)
  Bii_fe <- rep(0, NT)
  Bii_cov <- rep(0, NT)
  Bii_pe <- rep(0, NT)
  Mii <- rep(0, NT)
  Pii_sq <- rep(0, NT)
  Mii_sq <- rep(0, NT)
  Pii_Mii <- rep(0, NT)

  numIterations <- 300
  tol <- 1e-6

  if (type_algorithm == 'exact') {
    results <- matrix(0, nrow = NT, ncol = 4)
    pb <- NULL
    if (isTRUE(progress)) {
      pb <- txtProgressBar(min = 0, max = NT, style = 3)
    }
        for (i in 1:NT) {
          if (!is.null(pb)) {
            setTxtProgressBar(pb, i)
            if (i %% 10 == 0) flush.console()
          }

          xtilde <- sanic::solve_cg(a = xx, b = as.vector(X[i, , drop=FALSE]), type = "CG", iter = numIterations, tol = tol, precond = 2, verbose = FALSE)

          Pii_i <- as.numeric(X[i,] %*% xtilde)
          Bii_fe_i <- as.numeric(crossprod(X_fe %*% xtilde - mean(X_fe %*% xtilde)) / (nrow(X_fe) - 1))
          Bii_pe_i <- as.numeric(crossprod(X_pe %*% xtilde - mean(X_pe %*% xtilde)) / (nrow(X_pe) - 1))
          Bii_cov_i <- as.numeric(crossprod(X_pe %*% xtilde - mean(X_fe %*% xtilde)) / (nrow(X_pe) - 1))

          results[i,] <- c(Pii_i, Bii_fe_i, Bii_pe_i, Bii_cov_i)
        }
        if (!is.null(pb)) {
          close(pb)
        }

        Pii <- as.numeric(results[, 1])
        Bii_fe <- as.numeric(results[, 2])
        Bii_pe <- as.numeric(results[, 3])
        Bii_cov <- as.numeric(results[, 4])

        correction_JLA <- 1
        Mii <- 1 - Pii
      } else {
        pb <- NULL
        if (isTRUE(progress)) {
          pb <- txtProgressBar(min = 0, max = scale, style = 3)
        }
        for (i in 1:scale) {
      if (!is.null(pb)) {
        setTxtProgressBar(pb, i)
        if (i %% 10 == 0) flush.console()
      }
      # Random Projection Matrix (Rademacher)
      ons <- ifelse(runif(NT) > tolProb, 1, -1)

      # Get me the row
      Z <- sanic::solve_cg(a = xx, b = as.vector(ons %*% X), type = "CG", iter = numIterations, tol = tol, precond = 2, verbose = FALSE)
      Z <- X %*% Z

      # Collect (construct augmented estimator that combines Mii, Pii)
      Pii <- Pii + (Z^2) / scale
      Pii_sq <- Pii_sq + (Z^4) / scale
      Mii <- Mii + ((ons - Z)^2) / scale
      Mii_sq <- Mii_sq + ((ons - Z)^4) / scale
      Pii_Mii <- Pii_Mii + ((Z^2) * ((ons - Z)^2)) / scale

      # Random Projection Matrix (Rademacher)
      ons <- ifelse(runif(nrow(X_fe)) > tolProb, 1/sqrt(scale), -1/sqrt(scale))
      ons <- ons - mean(ons)

      # Bii for Variance of Firm Effects
      Z <- sanic::solve_cg(a = xx, b = as.vector(ons %*% X_fe), type = "CG", iter = numIterations, tol = tol, precond = 2, verbose = FALSE)
      Z_fe <- X %*% Z
      Bii_fe <- Bii_fe + (Z_fe * Z_fe)

      # Bii for Variance of Person Effects
      Z <- sanic::solve_cg(a = xx, b = as.vector(ons %*% X_pe), type = "CG", iter = numIterations, tol = tol, precond = 2, verbose = FALSE)
      Z_pe <- X %*% Z
      Bii_pe <- Bii_pe + (Z_pe * Z_pe)

      # Bii for CoVariance of Person, Firm Effects
      Z_cov <- (Z_pe * Z_fe)
      Bii_cov <- Bii_cov + Z_cov
    }

  # Element-wise operations
    Pii <- Pii / (Pii + Mii)
    Mii <- 1 - Pii
    Vi <- (1 / scale) * ((Mii^2) * Pii_sq + (Pii^2) * Mii_sq - 2 * Mii * Pii * Pii_Mii)
    Bi <- (1 / scale) * (Mii * Pii_sq - Pii * Mii_sq + 2 * (Mii - Pii) * Pii_Mii)
    correction_JLA <- 1 - Vi / (Mii^2) + Bi / Mii
    if (!is.null(pb)) {
      close(pb)
    }

  }

  return(list(Pii = Pii, Mii = Mii, correction_JLA = correction_JLA, Bii_fe = Bii_fe, Bii_cov = Bii_cov, Bii_pe = Bii_pe))
}
