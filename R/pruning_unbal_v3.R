#' Prune to a Leave-One-Worker-Out Connected Set
#'
#' Iteratively removes articulation workers from the worker-firm mobility graph
#' until the remaining sample stays connected after dropping any single worker.
#' This implements the leave-one-worker-out connectivity requirement used by the
#' main Kline, Saggio, and Solvsten (KSS) routines.
#'
#' @param y Numeric outcome vector.
#' @param firmid Firm identifier vector.
#' @param id Worker identifier vector.
#' @param id_old Original worker identifiers.
#' @param firmid_old Original firm identifiers.
#' @param controls Matrix of controls aligned with the observations.
#' @param prov_indicator Optional provider indicator carried along with the
#'   sample.
#' @param progress Logical scalar indicating whether iterative pruning progress
#'   should be emitted.
#'
#' @return A list containing the pruned outcome, identifiers, controls, and
#'   provider indicator.
#'
#' @details
#' The routine constructs a bipartite worker-firm graph for movers, identifies
#' articulation workers, removes them, and recomputes the largest connected
#' component until no articulation worker remains.
#'
#' @seealso [connected_set()], [build_adj()], [leave_out_KSS()]
#'
#' @export
pruning_unbal_v3 <- function(y, firmid, id, id_old, firmid_old, controls,
                             prov_indicator = rep(1, length(y)),
                             progress = FALSE) {
  .kss_validate_flag(progress, "progress")
  .kss_assert_same_length(
    y = y,
    firmid = firmid,
    id = id,
    id_old = id_old,
    firmid_old = firmid_old,
    prov_indicator = prov_indicator
  )
  controls <- .kss_validate_row_aligned_input(controls, length(y), "controls")
  .kss_stop_if_missing_or_nonfinite(y, "y")
  .kss_stop_if_missing_or_nonfinite(firmid, "firmid")
  .kss_stop_if_missing_or_nonfinite(id, "id")
  .kss_stop_if_missing_or_nonfinite(id_old, "id_old")
  .kss_stop_if_missing_or_nonfinite(firmid_old, "firmid_old")
  .kss_stop_if_missing_or_nonfinite(prov_indicator, "prov_indicator")
  .kss_stop_if_missing_or_nonfinite(controls, "controls")
  n_of_bad_workers <- 1

  while (n_of_bad_workers >= 1) {
    # find movers
J <- max(firmid, na.rm = TRUE)

# Directly calculate move without multiple assignments
gcs <- c(NA, head(id, -1))
lagfirmid <- c(NA, head(firmid, -1))
move <- ifelse(is.na(firmid) | is.na(lagfirmid) | id != gcs, 0, firmid != lagfirmid)
move <- as.integer(ave(move, id, FUN = sum, na.rm = TRUE)) > 0
move <- move[!is.na(id)]

# Keep original id_movers before normalization
id_movers_orig <- id[move]
.kss_stop_if_no_rows(length(id_movers_orig), "pruning_unbal_v3")

# Normalizing id_movers and firmid_movers
id_movers <- as.integer(factor(id_movers_orig))
firmid_movers <- firmid[move]
n = max(id_movers)
m = max(firmid_movers, na.rm = TRUE)
# Direct sparse matrix creation for bipartite graph
# Create the zero blocks without specifying 'i', 'j', or 'x'
      zero_block_n <- sparseMatrix(
        i = n,      # Row indices for non-zero entries
        j = n,  # Column indices for non-zero entries
        x = 0,  # Values at the non-zero entries (all 1s)
        dims = c(n, n)  # Dimensions of the matrix
      )
      zero_block_J <- sparseMatrix(
        i = m,      # Row indices for non-zero entries
        j = m,  # Column indices for non-zero entries
        x = 0,  # Values at the non-zero entries (all 1s)
        dims = c(m, m)  # Dimensions of the matrix
      )



B <- sparseMatrix(i = id_movers, j = firmid_movers, x = 1, dims = c(n, m))

# Construct the block matrix G using rbind and cbind to bind the matrices by rows and columns
      G <- rbind2(cbind2(zero_block_n, B), cbind2(t(B), zero_block_J))


# Find articulation points
g <- kss_graph_from_adjacency(G, diag = FALSE)
bad_workers <- articulation_points(g)[articulation_points(g) <= n]

   # Translate back to original worker IDs
sel <- id_movers %in% bad_workers
bad_workers <- unique(id_movers_orig[sel])
n_of_bad_workers <- length(bad_workers)
.kss_progress(progress, sprintf("Pruning: identified %s articulation workers", n_of_bad_workers))


# Remove bad workers and select largest connected set
sel <- !id %in% bad_workers
# Define a function to apply the selection based on data type
apply_selection <- function(data, sel) {
  if(is.matrix(data) || is.array(data)) {
    return(data[sel, , drop = FALSE])
  } else {
    return(data[sel])
  }
}

# Update all data vectors/matrices in one go
data_list <- list(y, firmid, id, firmid_old, id_old, controls, prov_indicator)
data_list <- lapply(data_list, apply_selection, sel = sel)

# Reset IDs
data_list[[2]] <- match(data_list[[2]], unique(data_list[[2]]))  # firmid
data_list[[3]] <- match(data_list[[3]], unique(data_list[[3]]))  # id

# Largest connected set once removed bad workers
A <- build_adj(data_list[[3]], data_list[[2]])
g <- kss_graph_from_adjacency(A, diag = FALSE)
firm_list <- which(components(g)$membership == which.max(components(g)$csize))
# Update all data vectors/matrices in one go
sel <- data_list[[2]] %in% firm_list
data_list <- lapply(data_list, apply_selection, sel = sel)

# Reset IDs one last time
data_list[[2]] <- match(data_list[[2]], unique(data_list[[2]]))  # firmid
data_list[[3]] <- match(data_list[[3]], unique(data_list[[3]]))  # id
# Extract back the data from data_list
y <- data_list[[1]]
firmid <- data_list[[2]]
id <- data_list[[3]]
firmid_old <- data_list[[4]]
id_old <- data_list[[5]]
controls <- data_list[[6]]
prov_indicator <- data_list[[7]]
.kss_stop_if_no_rows(length(y), "pruning_unbal_v3")
gc()
.kss_progress(progress, "Pruning: completed one leave-one-worker-out refinement step")
  }

  return(list(y = y, firmid = firmid, id = id, id_old = id_old, firmid_old = firmid_old, controls = controls, prov_indicator = prov_indicator))
}
