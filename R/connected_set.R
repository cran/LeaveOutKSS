#' Restrict a Panel to Its Largest Connected Set of Firms
#'
#' Builds a mobility graph from worker moves across firms and keeps only the
#' largest connected component of firms. This is the first graph-based trimming
#' step used by the leave-out routines before leave-one-worker-out pruning.
#'
#' @param y Numeric outcome vector.
#' @param id Worker identifier vector.
#' @param firmid Firm identifier vector.
#' @param lagfirmid Lagged firm identifier vector, typically constructed within
#'   worker.
#' @param controls Matrix of controls aligned with the observations.
#' @param prov_indicator Optional provider indicator carried along for interface
#'   compatibility.
#' @param progress Logical scalar indicating whether stage messages should be
#'   emitted.
#'
#' @return A list with two elements: `DT`, a `data.table` containing the
#'   restricted sample and original identifiers, and `DT_controls`, the
#'   correspondingly restricted controls.
#'
#' @details
#' The graph is built from observed worker transitions between lagged and
#' current firms. Firms not connected to the largest component are removed. The
#' function relabels worker and firm identifiers internally but preserves the
#' originals in the returned table.
#'
#' @seealso [pruning_unbal_v3()], [strongc_set()], [build_adj()]
#'
#' @export
connected_set <- function(y, id, firmid, lagfirmid, controls,
                          prov_indicator = rep(1, length(y)),
                          progress = FALSE) {
  .kss_validate_flag(progress, "progress")
  .kss_assert_same_length(
    y = y,
    id = id,
    firmid = firmid,
    lagfirmid = lagfirmid,
    prov_indicator = prov_indicator
  )
  controls <- .kss_validate_row_aligned_input(controls, length(y), "controls")
  .kss_stop_if_missing_or_nonfinite(y, "y")
  .kss_stop_if_missing_or_nonfinite(id, "id")
  .kss_stop_if_missing_or_nonfinite(firmid, "firmid")
  .kss_stop_if_missing_or_nonfinite(prov_indicator, "prov_indicator")
  .kss_stop_if_missing_or_nonfinite(controls, "controls")

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



    .kss_stop_if_no_rows(sum(sel), "connected_set")

    # Create adjacency matrix: sparseMatrix is already efficient
    A <- Matrix::sparseMatrix(i = lagfirmid[sel], j = firmid[sel], x = 1)

    # Make the matrix square by extending the smaller dimension
    max_dim <- max(dim(A))
    A <- Matrix::sparseMatrix(i = c(lagfirmid[sel], max_dim), j = c(firmid[sel], max_dim), x = c(rep(1, length(lagfirmid[sel])), 0))



  # Remove disconnected vertices: using sparse matrix operation
  degrees <- Matrix::rowSums(A) + Matrix::colSums(A)
  connected_firms <- which(degrees > 0)
  .kss_stop_if_no_rows(length(connected_firms), "connected_set")
  A <- A[connected_firms, connected_firms]  # subsetting sparse matrix is efficient

  # Create igraph object and find connected components
  g <- kss_graph_from_adjacency(A, diag = FALSE)
  # Directly find the largest connected component without intermediate steps
  components <- igraph::components(g)
  largest_component <- which.max(components$csize)

  # Filtering based on largest connected component's membership
  firmlst <- connected_firms[components$membership == largest_component]
  component_sel <- components$membership == largest_component
  A_component <- A[component_sel, component_sel]


  # Select only observations in the largest connected set
  sel <- firmid %in% firmlst
  .kss_stop_if_no_rows(sum(sel), "connected_set")
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
      "Connected set: mean graph degree %.3f",
      mean(kss_degree_from_adjacency(A_component, diag = FALSE))
    )
  )

  # Return a list containing both data.tables
  return(list(DT = DT_result, DT_controls = DT_controls_result))
}
