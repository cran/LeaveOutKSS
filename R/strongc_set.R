#' Restrict a Panel to Firms Above a Minimum Graph Degree Threshold
#'
#' Graph-based trimming helper that keeps firms whose degree in the mobility
#' graph is at least `min_degree`. This is a stronger restriction than the basic
#' connected-set filter and can be useful when the analyst wants a denser firm
#' network.
#'
#' @param y Numeric outcome vector.
#' @param id Worker identifier vector.
#' @param firmid Firm identifier vector.
#' @param controls Matrix of controls aligned with the observations.
#' @param min_degree Minimum graph degree required for a firm to remain in the
#'   sample.
#' @param progress Logical scalar indicating whether graph summary messages
#'   should be emitted.
#'
#' @return A list with `DT` and `DT_controls`, analogous to [connected_set()].
#'
#' @seealso [connected_set()]
#'
#' @export
strongc_set <- function(y, id, firmid, controls, min_degree = 1,
                        progress = FALSE) {
  .kss_validate_flag(progress, "progress")
  .kss_validate_scalar_number(min_degree, "min_degree", integerish = TRUE, min_value = 1)
  .kss_assert_same_length(y = y, id = id, firmid = firmid)
  controls <- .kss_validate_row_aligned_input(controls, length(y), "controls")
  .kss_stop_if_missing_or_nonfinite(y, "y")
  .kss_stop_if_missing_or_nonfinite(id, "id")
  .kss_stop_if_missing_or_nonfinite(firmid, "firmid")
  .kss_stop_if_missing_or_nonfinite(controls, "controls")
  
  # Assume that y, id, firmid, and controls are available and controls is a matrix
  DT <- data.table(y = y, id = id, firmid)
  DT_controls <- as.data.table(controls)
  
  # Lagfirmid
  DT[, lagfirmid := shift(firmid), by = id]
  DT[is.na(id) | shift(is.na(id), type = "lead", fill = TRUE), lagfirmid := NA]
  
  
  # Find connected set
  y <- DT$y
  id <- DT$id
  firmid <- DT$firmid
  lagfirmid <- DT$lagfirmid
  controls <- as.matrix(DT_controls)
  
  
  # Save old IDs
  firmid_old <- firmid
  id_old <- id

  # Create a selection vector for non-NA lagfirmid
  sel <- !is.na(lagfirmid)

  # Relabel firms
  firms <- unique(c(firmid, lagfirmid[sel]))
  ids <- unique(id)


  firmid <- match(firmid, firms)
  lagfirmid[sel] <- match(lagfirmid[sel], firms)
  id <- match(id, ids)
  
  
  .kss_stop_if_no_rows(sum(sel), "strongc_set")

  # Create adjacency matrix: sparseMatrix is already efficient
  A <- Matrix::sparseMatrix(i = lagfirmid[sel], j = firmid[sel], x = 1)
  
  
  # Make the matrix square by extending the smaller dimension
  max_dim <- max(dim(A))
  A <- Matrix::sparseMatrix(i = c(lagfirmid[sel], max_dim), j = c(firmid[sel], max_dim), x = c(rep(1, length(lagfirmid[sel])), 0))
  
  .kss_progress(
    progress,
    sprintf(
      "Strong connected set: adjusted adjacency dimensions %s x %s",
      dim(A)[1],
      dim(A)[2]
    )
  )
  
  
  # Remove disconnected vertices: using sparse matrix operation
  degrees <- Matrix::rowSums(A) + Matrix::colSums(A)
  connected_firms <- which(degrees > 1)
  .kss_stop_if_no_rows(length(connected_firms), "strongc_set")
  A <- A[connected_firms, connected_firms]  # subsetting sparse matrix is efficient
  
  # Create igraph object and find connected components
  g <- kss_graph_from_adjacency(A, diag = FALSE)
  d_g <- kss_degree_from_adjacency(A, diag = FALSE)
  # Filtering based on largest connected component's membership
  degree_sel <- d_g >= min_degree
  A_degree <- A[degree_sel, degree_sel]
  firmlst <- connected_firms[degree_sel]
  
  
  # Select only observations in the largest connected set
  sel <- firmid %in% firmlst
  .kss_stop_if_no_rows(sum(sel), "strongc_set")
  # Directly create data.tables within the function
  DT_result <- data.table(
    y = y[sel],
    id = id[sel],
    firmid = firmid[sel],
    id_old = id_old[sel],
    firmid_old = firmid_old[sel]
  )
  
  DT_controls_result <- as.data.table(controls[sel, , drop = FALSE])
  .kss_progress(
    progress,
    sprintf(
      "Strong connected set: mean graph degree %.3f",
      mean(kss_degree_from_adjacency(A_degree, diag = FALSE))
    )
  )
  
  # Return a list containing both data.tables
  return(list(DT = DT_result, DT_controls = DT_controls_result))
  
}
