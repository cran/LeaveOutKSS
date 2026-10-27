# Regression test for the exact leverage branch.
#
# Ground truth is computed with a dense inverse of the design with one firm
# column dropped. Within a connected set the person and firm effects are only
# identified up to a shared constant, and the Bii terms demean each projection,
# so they (and the leverages Pii) do not depend on that normalisation. This
# avoids a Moore-Penrose inverse and keeps the test base-R only.

make_panel <- function(N = 150, J = 10, Tp = 3, seed = 1) {
  set.seed(seed)
  id <- rep(seq_len(N), each = Tp)
  firmid <- sample(J, N * Tp, replace = TRUE)
  NT <- length(id)
  D <- Matrix::sparseMatrix(i = seq_len(NT), j = id, x = 1, dims = c(NT, N))
  Fm <- Matrix::sparseMatrix(i = seq_len(NT), j = firmid, x = 1, dims = c(NT, J))
  X <- cbind(D, -Fm)
  X_pe <- cbind(D, Matrix::Matrix(0, NT, J, sparse = TRUE))
  X_fe <- cbind(Matrix::Matrix(0, NT, N, sparse = TRUE), -Fm)
  xx <- methods::as(Matrix::crossprod(X), "generalMatrix")

  # Dense ground truth: drop the last firm column to get a full-rank design.
  keep <- seq_len(N + J - 1)
  Xr <- as.matrix(X[, keep])
  XA <- Xr %*% solve(crossprod(Xr)) # row i = xtilde_i' in reduced coordinates
  FE <- as.matrix(X_fe[, keep]) %*% t(XA) # column i = X_fe %*% xtilde_i
  PE <- as.matrix(X_pe[, keep]) %*% t(XA)
  FEc <- sweep(FE, 2, colMeans(FE))
  PEc <- sweep(PE, 2, colMeans(PE))

  list(
    X = X, X_pe = X_pe, X_fe = X_fe, xx = xx,
    truth = list(
      Pii = rowSums(XA * Xr),
      Bii_fe = colSums(FEc^2),
      Bii_pe = colSums(PEc^2),
      Bii_cov = colSums(PEc * FEc)
    )
  )
}

expect_leverages_match <- function(res, truth, tol = 1e-6) {
  expect_equal(as.vector(res$Pii), truth$Pii, tolerance = tol)
  expect_equal(as.vector(res$Mii), 1 - truth$Pii, tolerance = tol)
  expect_equal(as.vector(res$Bii_fe), truth$Bii_fe, tolerance = tol)
  expect_equal(as.vector(res$Bii_pe), truth$Bii_pe, tolerance = tol)
  expect_equal(as.vector(res$Bii_cov), truth$Bii_cov, tolerance = tol)
}

test_that("exact leverages() matches a dense ground truth", {
  p <- make_panel()
  res <- leverages(p$X_fe, p$X_pe, p$X, p$xx,
                   type_algorithm = "exact", scale = 1, progress = FALSE)
  expect_leverages_match(res, p$truth)
  expect_identical(res$correction_JLA, 1)
})

test_that("exact leverages_parallel() matches a dense ground truth", {
  skip_on_cran()
  p <- make_panel()
  res <- leverages_parallel(p$X_fe, p$X_pe, p$X, p$xx,
                            type_algorithm = "exact", scale = 1, progress = FALSE)
  expect_leverages_match(res, p$truth)
})
