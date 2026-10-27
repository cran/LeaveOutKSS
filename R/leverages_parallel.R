# Register parallel backend
#' Parallel Computation of Statistical Leverages and Bias Terms
#'
#' Parallel version of [leverages()] using `foreach` and `doParallel`.
#'
#' @inheritParams leverages
#'
#' @return A list with the same elements returned by [leverages()].
#'
#' @details
#' The exact and Johnson-Lindenstrauss approximation (JLA) branches mirror
#' [leverages()], but the repeated linear solves are distributed across worker
#' processes. This routine is intended for larger problems where the leverage
#' stage dominates runtime.
#'
#' @seealso [leverages()], [leave_out_KSS()]
#'
#' @export
leverages_parallel <- function(X_fe, X_pe, X, xx, type_algorithm, scale,
                               progress = FALSE) {
  NT <- nrow(X)
  tolProb=0.5;

  numIterations <- min(300L, max(dim(xx)))
  tol <- 1e-6

    .kss_progress(progress, "Parallel leverages: registering worker processes")
    registerDoParallel(cores = detectCores()-1)

    if (type_algorithm == 'exact') {

    results <- foreach(i = 1:NT,
                      .combine = rbind,
                      .packages = c("sanic", "Matrix")) %dopar%
    {
      xtilde <- sanic::solve_cg(a = xx, b = as.vector(X[i, , drop = FALSE]),
                                type = "CG", iter = numIterations, tol = tol, precond = 2,
                                verbose = FALSE)
      Pii_i <- as.numeric(X[i, ] %*% xtilde)
      v_fe <- as.vector(X_fe %*% xtilde)
      v_pe <- as.vector(X_pe %*% xtilde)
      Bii_fe_i <- sum((v_fe - mean(v_fe))^2)
      Bii_pe_i <- sum((v_pe - mean(v_pe))^2)
      Bii_cov_i <- sum((v_pe - mean(v_pe)) * (v_fe - mean(v_fe)))

      c(Pii_i, Bii_fe_i, Bii_pe_i, Bii_cov_i)
    }
    # Stop the parallel backend after the loop is done
    stopImplicitCluster()

  # Extract results
  Pii <- as.numeric(results[, 1])
  Bii_fe <- as.numeric(results[, 2])
  Bii_pe <- as.numeric(results[, 3])
  Bii_cov <- as.numeric(results[, 4])
  correction_JLA <- 1
  Mii <- 1 - Pii
    } else {

      # Custom combine function to sum corresponding elements of lists containing sparse matrices
      combine_func <- function(x, y) {
          # Iterate over the lists and add corresponding matrices element-wise
            mapply(function(a, b) a + b, x, y, SIMPLIFY = FALSE)
        }

        # The rest of the setup remains the same


      # Initialize with sparse matrices
      init_list <- replicate(8, Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(nrow(X), 1)), simplify = FALSE)

      # Set up the parallel backend
      .kss_progress(progress, "Parallel leverages: starting JLA worker tasks")
      registerDoParallel(cores = detectCores()-1 )

      # Parallel loop
      results <- foreach(i = 1:scale,
                   .combine = 'combine_func',
                   .init = init_list,
                   .multicombine = FALSE,
                   .packages = c("sanic", "Matrix")) %dopar% {
    # Random Projection Matrix (Rademacher)
    #ons <- ifelse(runif(nrow(X)) > tolProb, 1, -1)
    ons <- rbinom(NT, 1, 1 - tolProb) * 2 - 1


    # Get me the row
    Z <- sanic::solve_cg(a = xx, b = as.vector(ons %*% X), type = "CG", iter = numIterations, tol = tol, precond = 2, verbose = FALSE)
    Z <- X %*% Z

    # Collect
    Pii_i <- (Z^2) / scale
    Pii_sq_i <- (Z^4) / scale
    ons_Z <- ons - Z
    Mii_i <- (ons_Z^2) / scale
    Mii_sq_i <- (ons_Z^4) / scale
    Pii_Mii_i <- (Z^2 * ons_Z^2) / scale

    # Random Projection Matrix (Rademacher)
    # Adjust ons in place for Variance of Firm and Person Effects
    ons <- (rbinom(nrow(X_fe), 1, 1 - tolProb) *2 - 1)/sqrt(scale)
    ons <- ons - mean(ons)

    # Bii for Variance of Firm Effects
    Z_fe <- sanic::solve_cg(a = xx, b = as.vector(ons %*% X_fe), type = "CG", iter = numIterations, tol = tol, precond = 2, verbose = FALSE)
    Bii_fe_i <- (X %*% Z_fe)^2

    # Bii for Variance of Person Effects
    Z_pe <- sanic::solve_cg(a = xx, b = as.vector(ons %*% X_pe), type = "CG", iter = numIterations, tol = tol, precond = 2, verbose = FALSE)
    Bii_pe_i <- (X %*% Z_pe)^2

    # Bii for CoVariance of Person, Firm Effects
    Bii_cov_i <- ((X %*% Z_pe) * (X %*% Z_fe))

    list(Pii_i, Pii_sq_i, Mii_i, Mii_sq_i, Pii_Mii_i, Bii_fe_i, Bii_pe_i, Bii_cov_i)
  }
  # Stop the parallel backend after the loop is done
    stopImplicitCluster()

  # Extract results
  Pii <- results[[1]]
  Pii_sq <- results[[2]]
  Mii <- results[[3]]
  Mii_sq <- results[[4]]
  Pii_Mii <- results[[5]]
  Bii_fe <- results[[6]]
  Bii_pe <- results[[7]]
  Bii_cov <- results[[8]]

  # Element-wise operations
  Pii <- Pii / (Pii + Mii)
  Mii <- 1 - Pii
  Vi <- (1 / scale) * ((Mii^2) * Pii_sq + (Pii^2) * Mii_sq - 2 * Mii * Pii * Pii_Mii)
  Bi <- (1 / scale) * (Mii * Pii_sq - Pii * Mii_sq + 2 * (Mii - Pii) * Pii_Mii)
  correction_JLA <- 1 - Vi / (Mii^2) + Bi / Mii
    }
  # Return results
  return(list(Pii = Pii, Mii = Mii, correction_JLA = correction_JLA, Bii_fe = Bii_fe, Bii_cov = Bii_cov, Bii_pe = Bii_pe))
}
