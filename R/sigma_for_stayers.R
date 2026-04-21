#' Approximate Leave-Out Variance Terms for Stayers
#'
#' Computes the stayer-specific adjustment used when the main decomposition is
#' performed at the match level. In that case, the current implementation uses a
#' leave-one-observation-out style adjustment for stayers, following the
#' approximation discussed in the original vignette.
#'
#' @param y Outcome vector in person-year space.
#' @param id Worker identifier vector in collapsed match space.
#' @param firmid Firm identifier vector in collapsed match space.
#' @param peso Match weights used to expand back to person-year space.
#' @param b Estimated coefficient vector from the worker-firm fixed effects
#'   regression.
#'
#' @return A vector of averaged stayer variance adjustments at the match level.
#'
#' @seealso [leave_out_KSS()], [leave_out_KSS_fe()]
#'
#' @export
sigma_for_stayers <- function(y, id, firmid, peso, b) {
  # Computes heteroskedastic robust variances for stayers

  # Get back to person-year space
  # Expanding the 'id' and 'firmid' vectors to match the expanded 'y' vector
  id <- rep(id, times = peso)
  firmid <- rep(firmid, times = peso)

  # Compute Pii for stayers (1/T_i)
  # Count the number of observations per 'id'
  T <- tapply(rep(1, length(id)), id, sum)
  # Inverse of T to find Pii, using the expanded 'id' vector
  Pii <- 1 / T[id]

  # Compute OLS residuals
  NT <- length(y)
  # Using the expanded vectors to create the design matrices
  D <- sparseMatrix(i = 1:NT, j = id, x = 1)
  F <- sparseMatrix(i = 1:NT, j = firmid, x = 1)

  # Create a restriction matrix 'S' if needed, assuming 'J' is the number of firms
  # Assuming you're applying some kind of restriction like in your 'MATLAB' code
  J <- dim(F)[2]
  S <- Diagonal(x = 1, n = J - 1)
  S <- rbind2(S, -rowSums(S))  # Assuming a similar restriction as in 'MATLAB' code

  # Construct the model matrix 'X'
  X <- cbind2(D, F %*% S)
  # Calculate the OLS residuals
  eta <- y - as.vector(X %*% b)

  # Compute mean of 'y' once
  y_mean <- mean(y)

  # Adjust residuals for heteroskedasticity
  # Ensure 'Pii' is a vector here
  eta_h <- eta / (1 - Pii)

  # Compute sigma_stayers using heteroskedasticity adjusted residuals
  # Reuse 'y_mean' instead of calling 'mean(y)' repeatedly
  sigma_stayers <- (y - y_mean) * eta_h

  # Collapse back to match level by averaging over matches
  # Using 'data.table' for more efficient aggregation
  DT <- data.table(id = id, firmid = firmid, sigma_stayers = sigma_stayers)
  sigma_stayers_collapsed <- DT[, .(mean_sigma_stayers = mean(sigma_stayers)), by = .(id, firmid)]


  return(sigma_stayers_collapsed$mean_sigma_stayers)
}
