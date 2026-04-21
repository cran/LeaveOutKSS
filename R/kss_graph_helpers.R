#' Internal Helpers for Building Undirected Graphs from Sparse Adjacency Matrices
#'
#' These helpers avoid `igraph::graph_from_adjacency_matrix()` for very large
#' sparse mobility graphs by extracting nonzero edges directly from a sparse
#' adjacency matrix and then constructing an undirected `igraph` object from the
#' resulting edge list.
#'
#' They are used internally by the connected-set and pruning routines.
#'
#' @keywords internal
kss_sparse_undirected_edges <- function(A, diag = FALSE) {
  if (!inherits(A, "Matrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  }

  if (nrow(A) != ncol(A)) {
    stop("Adjacency matrix must be square.")
  }

  edge_dt <- data.table::as.data.table(Matrix::summary(A))

  if (nrow(edge_dt) == 0L) {
    return(data.table::data.table(u = integer(), v = integer(), w = numeric()))
  }

  if (!diag) {
    edge_dt <- edge_dt[i != j]
  }

  if (nrow(edge_dt) == 0L) {
    return(data.table::data.table(u = integer(), v = integer(), w = numeric()))
  }

  edge_dt[, `:=`(u = pmin(i, j), v = pmax(i, j))]
  edge_dt <- edge_dt[, .(w = max(x)), by = .(u, v)]
  edge_dt[]
}

#' @keywords internal
kss_graph_from_adjacency <- function(A, diag = FALSE) {
  vc <- nrow(A)
  edge_dt <- kss_sparse_undirected_edges(A, diag = diag)
  g <- igraph::make_empty_graph(n = vc, directed = FALSE)

  if (nrow(edge_dt) == 0L) {
    return(g)
  }

  edge_mat <- as.matrix(edge_dt[, .(u, v)])
  igraph::add_edges(g, edges = c(t(edge_mat)))
}

#' @keywords internal
kss_degree_from_adjacency <- function(A, diag = FALSE) {
  vc <- nrow(A)
  edge_dt <- kss_sparse_undirected_edges(A, diag = diag)

  if (nrow(edge_dt) == 0L) {
    return(numeric(vc))
  }

  d <- numeric(vc)
  deg_u <- edge_dt[, .(w = sum(w)), by = u]
  deg_v <- edge_dt[, .(w = sum(w)), by = v]
  d[deg_u$u] <- d[deg_u$u] + deg_u$w
  d[deg_v$v] <- d[deg_v$v] + deg_v$w
  d
}
