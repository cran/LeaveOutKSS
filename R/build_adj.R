
#' Build a Firm-to-Firm Mobility Adjacency Matrix
#'
#' Constructs a symmetric sparse adjacency matrix of firm mobility links using
#' worker transitions. Only movers contribute edges.
#'
#' @param id Worker identifier vector.
#' @param firmid Firm identifier vector.
#'
#' @return A sparse square adjacency matrix whose nonzero entries count observed
#'   worker moves between firms.
#'
#' @seealso [connected_set()], [pruning_unbal_v3()]
#'
#' @examples
#' build_adj(
#'   id = c(1, 1, 2, 2, 3, 3),
#'   firmid = c(1, 2, 2, 3, 3, 3)
#' )
#'
#' @export
build_adj <- function(id, firmid) {

  # Initialize data.table
  DT <- data.table(id = id, firmid = firmid)

  # Generate gcs, a logical vector indicating the start of a new id
  DT[, gcs := c(NA, head(id, -1)) != id]

  # Generate lagfirmid, the previous firm id for each observation
  DT[, lagfirmid := shift(firmid)]
  DT[gcs == TRUE, lagfirmid := NA]

  # Identify stayer periods, where firm id doesn't change
  DT[, stayer := (!is.na(firmid) & !is.na(lagfirmid) & firmid == lagfirmid)]
  DT[gcs == TRUE, stayer := TRUE]

  # Summarize the number of stayers and total periods per id
  DT[, `:=`(stayer_sum = sum(stayer, na.rm = TRUE), T = .N), by = id]
  DT[, stayer := stayer_sum == T]
  DT[, movers := !stayer]

  # Subset data for movers
  DT_movers <- DT[movers == TRUE]

  # Calculate lagfirmid and gcs for movers
  DT_movers[, lagfirmid := shift(firmid)]
  DT_movers[gcs == TRUE, lagfirmid := NA]

  # Filter for valid moves
  list_data <- DT_movers[!is.na(lagfirmid) & lagfirmid != firmid, .(lagfirmid, firmid, id)]

  # Calculate edge weights
  edge_weights <- table(paste0(list_data$lagfirmid, "-", list_data$firmid))

  # Build sparse adjacency matrix
  e <- unique(list_data[, .(lagfirmid, firmid)])
  J <- max(DT$firmid, na.rm = TRUE)
  A <- sparseMatrix(i = c(e$lagfirmid, e$firmid),
                    j = c(e$firmid, e$lagfirmid),
                    x = c(as.numeric(edge_weights), as.numeric(edge_weights)),
                    dims = c(J, J))

  return(A)
}
