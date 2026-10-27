#' Leave-Out Bias-Corrected Variance Decomposition in a Two-Way Fixed Effects Model
#'
#' Estimates plug-in and leave-out bias-corrected variance components for a
#' two-way fixed effects model as part of the R translation of the original
#' 'MATLAB' package of Kline, Saggio, and Solvsten (2020). The function starts
#' from worker identifiers, firm identifiers, and an outcome, constructs the
#' leave-one-worker-out connected set, optionally partials out controls,
#' computes statistical leverages either exactly or via the
#' Johnson-Lindenstrauss approximation (JLA), and returns decomposition
#' summaries together with estimated worker and firm effects.
#'
#' Relative to the original 'MATLAB' package, this implementation follows the
#' same broad sequence: connected-set construction, leave-out pruning,
#' optional residualization of controls, leverage computation, and bias
#' correction of the variance of firm effects, the covariance of worker and
#' firm effects, and the variance of worker effects.
#'
#' @param y Numeric outcome vector.
#' @param id Worker identifier vector.
#' @param firmid Firm identifier vector.
#' @param controls Optional matrix or vector of controls. When supplied, the
#'   function prepends an intercept internally and residualizes the outcome with
#'   respect to worker, firm, and control regressors before computing variance
#'   components.
#' @param leave_out_level Character scalar. Use `"matches"` to leave out entire
#'   worker-firm matches or `"obs"` to leave out person-year observations.
#' @param type_algorithm Character scalar. Use the randomized
#'   Johnson-Lindenstrauss approximation (`"JLA"`) to the leverages or
#'   `"exact"` for the exact algorithm.
#' @param simulations_JLA Integer number of random projections when
#'   `type_algorithm = "JLA"`.
#' @param lincom_do Integer flag equal to `0` or `1`. When `1`, the function
#'   also calls [lincom_KSS()] to regress estimated firm effects on user-supplied
#'   observables.
#' @param Z_lincom Optional matrix of observables used by [lincom_KSS()] when
#'   `lincom_do = 1`. The main decomposition may collapse the estimation sample
#'   to match means, but the optional lincom step is still run on the
#'   post-pruning observation-level sample so observation weights are preserved.
#' @param labels_lincom Optional labels for the columns of `Z_lincom`.
#' @param csv_file Optional path for exporting the estimated effects table as a
#'   `.csv` file.
#' @param txt_file Optional path for exporting a text summary of the
#'   decomposition.
#' @param paral Logical scalar indicating whether leverage computation should use
#'   the parallel routine [leverages_parallel()].
#' @param Cd Integer random seed passed to [base::set.seed()].
#' @param progress Logical scalar indicating whether stage progress messages
#'   should be emitted.
#'
#' @return An object of class `"leave_out_kss_result"` containing biased and
#'   bias-corrected estimates, estimated worker and firm effects, optional
#'   lincom results, sample summaries, and elapsed time.
#'
#' @details
#' The decomposition is based on an Abowd, Kramarz, and Margolis (1999;
#' AKM)-style model with worker effects, firm effects, and optional controls.
#' By default, the function leaves out matches, which corresponds to allowing
#' unrestricted heteroskedasticity and arbitrary serial correlation within
#' worker-firm matches, in line with the discussion in the original vignette.
#' When `leave_out_level = "obs"`, the correction is based on leaving out one
#' person-year observation at a time.
#'
#' When controls are supplied, the function first estimates their coefficients
#' in the leave-out connected set and then works with the residualized outcome.
#' When `lincom_do = 1`, the function additionally reports linear projections of
#' firm effects on observables using [lincom_KSS()].
#'
#' The input vectors must be sorted by worker identifier and, within worker,
#' from earlier to later time periods before calling the function. When
#' `controls` or `Z_lincom` are supplied, they must follow that same sorted row
#' order.
#'
#' The returned object is the primary estimation record. It stores the
#' decomposition summaries, estimated worker and firm effects, and optional
#' lincom output. When `csv_file` or `txt_file` are supplied, those summaries
#' are also written to disk.
#'
#' @references
#' Kline, P., Saggio, R., and Solvsten, M. (2020). Leave-out estimation of
#' variance components. *Econometrica*, 88(5), 1859-1898.
#'
#' Abowd, J. M., Kramarz, F., and Margolis, D. N. (1999). High wage workers and
#' high wage firms. *Econometrica*, 67(2), 251-333.
#'
#' @seealso [leave_out_KSS_fe()], [rsquared_comp()], [lincom_KSS()],
#'   [leverages()], [leverages_parallel()]
#'
#' @examples
#' path <- system.file("extdata", "test.csv", package = "LeaveOutKSS")
#' dt <- data.table::fread(path, header = FALSE)
#' data.table::setorder(dt, V1, V3)
#'
#' res <- leave_out_KSS(
#'   y = dt[[4]],
#'   id = dt[[1]],
#'   firmid = dt[[2]],
#'   simulations_JLA = 5,
#'   paral = FALSE,
#'   progress = FALSE
#' )
#'
#' print(res)
#'
#' @export
leave_out_KSS <- function(y, id, firmid,
                          controls = NULL,
                          leave_out_level = 'matches',
                          type_algorithm = 'JLA',
                          simulations_JLA = 200,
                          lincom_do = 0,
                          Z_lincom = NULL,
                          labels_lincom = NULL,
                          csv_file = NULL,
                          txt_file = NULL,
                          paral = TRUE,
                          Cd = 12345,
                          progress = FALSE) {
  csv_file <- .kss_validate_optional_path(csv_file, "csv_file")
  txt_file <- .kss_validate_optional_path(txt_file, "txt_file")
  no_controls <- 0
  start_time <- proc.time()[["elapsed"]]
  set.seed(Cd)

  # Check for minimum required arguments
  if (missing(y) || missing(id) || missing(firmid)) {
    stop("`y`, `id`, and `firmid` must all be supplied.", call. = FALSE)
  }
  .kss_validate_flag(paral, "paral")
  .kss_validate_flag(progress, "progress")
  .kss_validate_scalar_number(simulations_JLA, "simulations_JLA", integerish = TRUE, min_value = 1)
  .kss_validate_scalar_number(Cd, "Cd")
  .kss_assert_same_length(y = y, id = id, firmid = firmid)
  .kss_stop_if_missing_or_nonfinite(y, "y")
  .kss_stop_if_missing_or_nonfinite(id, "id")
  .kss_stop_if_missing_or_nonfinite(firmid, "firmid")
  leave_out_level <- match.arg(leave_out_level, c("matches", "obs"))
  type_algorithm <- match.arg(type_algorithm, c("JLA", "exact"))

  if (!lincom_do %in% c(0, 1)) {
    stop("`lincom_do` must be either 0 or 1.", call. = FALSE)
  }
  if (lincom_do == 1 && is.null(Z_lincom)) {
    stop("`Z_lincom` must be supplied when `lincom_do = 1`.", call. = FALSE)
  }
  if (!is.null(labels_lincom)) {
    labels_lincom <- as.character(unlist(labels_lincom, use.names = FALSE))
  }

  controls <- .kss_validate_row_aligned_input(controls, length(y), "controls")
  if (!is.null(controls)) {
    .kss_stop_if_missing_or_nonfinite(controls, "controls")
  }
  Z_lincom <- .kss_validate_row_aligned_input(Z_lincom, length(y), "Z_lincom")
  if (!is.null(Z_lincom)) {
    .kss_stop_if_missing_or_nonfinite(Z_lincom, "Z_lincom")
  }

  if (is.null(controls)) {
    controls <- matrix(1, nrow = length(y), ncol = 1)
    no_controls <- 1
  } else {
    controls <- cbind(1, controls)
  }

  .kss_progress(progress, sprintf("leave_out_KSS: starting %s leave-out decomposition", leave_out_level))
  .kss_progress(progress, sprintf("leave_out_KSS: leverage method %s", type_algorithm))

  # Relabel potentially large IDs to dense indices for sparse matrix routines.
  id_original <- id
  firmid_original <- firmid
  id <- match(id, unique(id))
  firmid <- match(firmid, unique(firmid))

  # Assume that y, id, firmid, and controls are available and controls is a matrix
  DT <- data.table(y = y, id = id, firmid = firmid)
  DT_controls <- as.data.table(controls)

# STEP 1: FIND CONNECTED SET

# Lagfirmid
DT[, lagfirmid := shift(firmid), by = id]
DT[is.na(id) | shift(is.na(id), type = "lead", fill = TRUE), lagfirmid := NA]

# Update controls if lincom_do is 1
if (lincom_do == 1) {
  K=ncol(controls)
  DT_controls <- cbind(DT_controls, Z_lincom)
}

  .kss_progress(progress, "leave_out_KSS: restricting to the largest connected set")
  connected_set_result <- connected_set(
    y = DT$y,
    id = DT$id,
    firmid = DT$firmid,
    lagfirmid = DT$lagfirmid,
    controls = as.matrix(DT_controls),
    progress = FALSE
  )

# Directly updating DT and DT_controls
DT <- connected_set_result$DT
DT_controls <- connected_set_result$DT_controls
rm(connected_set_result)
gc()
.kss_stop_if_no_rows(nrow(DT), "leave_out_KSS")

  # STEP 2: LEAVE ONE OUT CONNECTED SET
  .kss_progress(progress, "leave_out_KSS: pruning to a leave-one-worker-out connected set")
  output <- pruning_unbal_v3(
    y = DT$y,
    firmid = DT$firmid,
    id = DT$id,
    id_old = DT$id_old,
    firmid_old = DT$firmid_old,
    controls = DT_controls,
    progress = progress
  )

# Create new data.tables
DT <- data.table(y = output$y,
                 firmid = output$firmid,
                 id = output$id,
                 id_old = output$id_old,
                 firmid_old = output$firmid_old)

DT_controls <- as.data.table(output$controls)
rm(output)
gc()
.kss_stop_if_no_rows(nrow(DT), "leave_out_KSS")

# Drop stayers with a single person year observation
DT[, count := .N, by = id]
DT_controls<-DT_controls[DT$count > 1]
DT <- DT[count > 1]
.kss_stop_if_no_rows(nrow(DT), "leave_out_KSS")


# Resetting ids
DT[, firmid := match(firmid, unique(firmid))]
DT[, id := match(id, unique(id))]
post_pruning_adj <- build_adj(DT$id, DT$firmid)

    # Important Auxiliaries
# Create a flag for the first observation for each worker
DT[, first_obs := .I == 1, by = id]

# Create a lagged firmid, with NA for the first observation per worker
DT[, lagfirmid := shift(firmid), by = id]
DT[first_obs==TRUE, lagfirmid := NA]

# Identify stayers: same firmid as lagfirmid, flag first observation as stayer
DT[, stayer := (firmid == lagfirmid) | first_obs]
DT[is.na(stayer), stayer:=TRUE]

# Sum up stayers by id to determine if all records are stayers
DT[, stayer_sum := sum(stayer, na.rm = TRUE), by = id]
DT[, T := .N, by = id]
DT[, stayer := stayer_sum == T]

# Identify movers: not a stayer in any record
DT[, movers := !stayer]

# Calculate the number of unique movers
id_movers <- DT[movers==TRUE, id]
Nmovers <- uniqueN(id_movers) # Number of unique movers

NT <- nrow(DT)

D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1)
N <- ncol(D)
F <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1)
J <- ncol(F)
var_den <- var(DT$y, na.rm = TRUE)

  leave_one_out_stats <- list(
    mean_outcome = mean(DT$y, na.rm = TRUE),
    variance_outcome = var(DT$y, na.rm = TRUE),
    n_movers = Nmovers,
    n_firms = max(DT$firmid, na.rm = TRUE),
    n_observations = length(DT$y)
  )
  .kss_progress(
    progress,
    sprintf(
      "leave_out_KSS: mean graph degree after pruning %.3f",
      mean(kss_degree_from_adjacency(post_pruning_adj, diag = FALSE))
    )
  )


  # STEP 3: Residualizing and Collapsing
  .kss_progress(progress, "leave_out_KSS: adjusting for controls")

# Assuming D, F, controls, y, J, N, and no_controls are already defined
if (lincom_do == 1) {
  Z <- as.matrix(DT_controls[, (K + 1):ncol(DT_controls), with = FALSE])
  DT_controls <- DT_controls[, 1:K, with = FALSE]
}
lincom_match_id <- NULL

y_original <- DT$y


if (no_controls == 0) {
  # Creating the restriction matrix S
  S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(1, J-1)))

  # Constructing the X matrix
  X <- cbind2(cbind2(D, F %*% S), as.matrix(DT_controls))

  # Calculating xx (X' * X) and xy (X' * y)
  xx <- t(X) %*% X
  xy <- t(X) %*% DT$y



  # Solving for beta
  b <- sanic::solve_cg(a = xx, b = xy, type = "CG", iter = min(1000L, max(dim(xx))), tol = 1e-10, precond = 2, verbose = FALSE)
  b <- as.vector(b)



  # Variance decomposition based on residualized outcome
  y_c_corrected<-DT$y - as.vector(X[, (N + J):ncol(X)] %*% b[(N + J):length(b)])
  DT[, y := y_c_corrected]
}


  control_adjusted_stats <- list(
    mean_outcome = mean(DT$y, na.rm = TRUE),
    variance_outcome = var(DT$y, na.rm = TRUE)
  )

# Collapsing
peso <- rep(1, nrow(DT))
y_py <- DT$y

if (leave_out_level == 'matches') {
  # Generate match_id
  DT[, match_id := as.numeric(factor(paste(id, firmid, sep = "_")))]
  if (lincom_do == 1) {
    lincom_match_id <- DT$match_id
  }

  # Calculate weights based on the length of a given match
  DT[, peso := .N, by = match_id]

  # Collapse data down to match-means
  DT <- DT[, .(
  id = mean(id),
  firmid = mean(firmid),
  id_old = mean(id_old),
  firmid_old = mean(firmid_old),
  y = mean(y),
  peso = .N  # This counts the number of rows in each group
), by = .(match_id)]
peso <- DT$peso
}

  # STEP 4: COMPUTATION OF (Pii, Bii)
  .kss_progress(progress, "leave_out_KSS: computing leverages and bias terms")

# Build Design
NT <- nrow(DT)
D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1)
N <- ncol(D)
FM <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1)
J <- ncol(FM)
X <- cbind2(D, -FM)

# Debugging messages
# Weighting Matrices
X_fe <- cbind2(sparseMatrix(i = integer(0), j = integer(0), x=0 , dims = c(NT, N)) , X[, (N + 1):(N + J)])
X_pe <- cbind2(X[, 1:N], sparseMatrix(i = integer(0), j = integer(0), x=0 , dims = c(NT, J)))
# Replicate rows according to the weights 'peso'
row_indices_fe <- rep(1:NT, times = peso)
row_indices_pe <- rep(1:NT, times = peso)

X_fe <- X_fe[row_indices_fe, , drop = FALSE]
X_pe <- X_pe[row_indices_pe, , drop = FALSE]

# FGLS Transformation
PESO_MAT <- sparseMatrix(i = 1:NT, j = 1:NT, x = sqrt(peso))


y_untransformed <- DT$y

X <- PESO_MAT %*% X

y <- PESO_MAT %*% DT$y

  # Compute xx (X' * X)
  xx <- t(X) %*% X
  if (paral) {
    PBii <- leverages_parallel(
      X_fe = X_fe,
      X_pe = X_pe,
      X = X,
      xx = xx,
      type_algorithm = type_algorithm,
      scale = simulations_JLA,
      progress = progress
    )
  } else {
    PBii <- leverages(
      X_fe = X_fe,
      X_pe = X_pe,
      X = X,
      xx = xx,
      type_algorithm = type_algorithm,
      scale = simulations_JLA,
      progress = progress
    )
  }



# Extract the results
#Pii <- results$Pii
#Mii <- results$Mii
#correction_JLA <- results$correction_JLA
#Bii_fe <- results$Bii_fe
#Bii_cov <- results$Bii_cov
#Bii_pe <- results$Bii_pe

  # STEP 5: ESTIMATION OF VARIANCE COMPONENTS
  # We use the statistical leverages, Pii, and the Bii associated with a given
  # variance component to bias correct these quantities using the KSS approach
  .kss_progress(progress, "leave_out_KSS: estimating biased and bias-corrected components")
# Bind feature matrices and add debugging message
X_fe_lincom <- cbind2(sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(NT, N)), FM)
X_fe <- X_fe_lincom[rep(1:NT, times = peso), , drop = FALSE]


# Bind prediction matrices and add debugging message
X_pe <- cbind2(D, sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(NT, J)))
X_pe <- X_pe[rep(1:NT, times = peso), , drop = FALSE]

# Bind S matrix and add debugging message
S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(1, J-1)))

# Combine D and F matrices, and add debugging message
X_lincom <- cbind2(D, FM %*% S)
y_lincom <- DT$y
X <- PESO_MAT %*% X_lincom

# Compute crossproducts and add debugging messages
xx <- t(X) %*% X
xy <- t(X) %*% y

# Execute the Conjugate Gradient solver and add debugging message
b <- sanic::solve_cg(a = xx, b = as.vector(xy), type = "CG", iter = min(1000L, max(dim(xx))), tol = 1e-10, precond = 2, verbose = FALSE)


# Residual calculation and KSS estimate
eta <- y - X %*% b
eta_h <- eta / PBii$Mii  # Leave one out residual
sigma_i <- (y - mean(y)) * eta_h  # KSS estimate of individual variance
sigma_i <- sigma_i * PBii$correction_JLA # Adjustment for non-linear bias

# Drop the last column from the fixed and random effects design matrices
X_fe_lincom <- X_fe_lincom[, -ncol(X_fe_lincom), drop = FALSE]
X_fe <- X_fe[, -ncol(X_fe), drop = FALSE]
X_pe <- X_pe[, -ncol(X_pe), drop = FALSE]

if (leave_out_level == 'matches') {
  # Calculate the number of times each 'id' is present in the dataset
  DT[, T := .N, by = id]


  # Determine 'stayers' (where T == 1)
  DT[, stayers := T == 1]
  # Calculate 'sigma_stayers' using your 'sigma_for_stayers' function
  # Assuming 'sigma_for_stayers' is a vectorized function that returns a vector of the same length
  sigma_stayers <- sigma_stayers <- sigma_for_stayers(y = y_py,
                                   id = DT$id,
                                   firmid = DT$firmid,
                                   peso = peso,
                                   b = b)
  sigma_i[DT$stayers]<-sigma_stayers[DT$stayers]
  rm(sigma_stayers)
  # Update 'sigma_i' for stayers
}






# Compute the quadratic form for the first set and assign the values to variables
result <- kss_quadratic_form(sigma_i, X_fe, X_fe, b, PBii$Bii_fe)
sigma_2_psi_AKM <- result[[1]]
sigma2_psi <- result[[2]]
# Compute the quadratic form for the second set and assign the values to variables
result <- kss_quadratic_form(sigma_i, X_fe, X_pe, b, PBii$Bii_cov)
sigma_alpha_psi_AKM <- result[[1]]
sigma_psi_alpha <- result[[2]]

# Compute the quadratic form for the third set and assign the values to variables
result <- kss_quadratic_form(sigma_i, X_pe, X_pe, b, PBii$Bii_pe)
sigma_2_alpha_AKM <- result[[1]]
sigma2_alpha <- result[[2]]
rm(result) # Remove the result from memory

  correlation_biased <- sigma_alpha_psi_AKM / (sqrt(sigma_2_psi_AKM) * sqrt(sigma_2_alpha_AKM))
  fraction_variance_explained_biased <- (sigma_2_psi_AKM + 2 * sigma_alpha_psi_AKM + sigma_2_alpha_AKM) / var_den
  correlation_bias_corrected <- sigma_psi_alpha / (sqrt(sigma2_psi) * sqrt(sigma2_alpha))
  fraction_variance_explained_bias_corrected <- (sigma2_psi + 2 * sigma_psi_alpha + sigma2_alpha) / var_den

  biased_estimates <- c(
    variance_firm_effects = sigma_2_psi_AKM,
    covariance_firm_person_effects = sigma_alpha_psi_AKM,
    variance_person_effects = sigma_2_alpha_AKM,
    correlation_firm_person_effects = correlation_biased,
    fraction_variance_explained = fraction_variance_explained_biased
  )
  bias_corrected_estimates <- c(
    variance_firm_effects = sigma2_psi,
    covariance_firm_person_effects = sigma_psi_alpha,
    variance_person_effects = sigma2_alpha,
    correlation_firm_person_effects = correlation_bias_corrected,
    fraction_variance_explained = fraction_variance_explained_bias_corrected
  )


  # STEP 7: CREATING OUTPUT DATA
  .kss_progress(progress, "leave_out_KSS: assembling estimated effects")
OrigD <- data.table(
  id = id,
  firmid = firmid,
  id_original = id_original,
  firmid_original = firmid_original
)
OrigD <- OrigD[, order_id := .I]
setkey(OrigD,id,firmid)

NewD <- data.table(id = rep(DT$id_old,times = peso),
                   pe = as.vector(X_pe %*% b),  # Ensuring the result is a vector
                   firmid = rep(DT$firmid_old, times = peso),
                   fe = as.vector(X_fe %*% b))
NewD <- unique(NewD , by = c("id","firmid"))
setkey(NewD,id,firmid)


OrigD <- NewD[OrigD, nomatch = 0]
setorder(OrigD, order_id)

if (no_controls==0) {
    output_dt <- data.table(
    y_original = y_original,
    y_control_adjusted = y_c_corrected,
    id = OrigD$id_original,
    pe = OrigD$pe,  # Ensuring the result is a vector
    firm_id = OrigD$firmid_original,
    fe = OrigD$fe
  )
} else{
    output_dt <- data.table(
    y_original = y_original,
    id = OrigD$id_original,
    pe = OrigD$pe,  # Ensuring the result is a vector
    firm_id = OrigD$firmid_original,
    fe = OrigD$fe
  )
}


  # STEP 8: LINCOM
  lincom_result <- NULL
  if (lincom_do == 1) {
    if (leave_out_level == "matches") {
      lincom_row_index <- match(lincom_match_id, DT$match_id)
      if (anyNA(lincom_row_index)) {
        stop("Could not expand the match-collapsed lincom inputs back to observation level.", call. = FALSE)
      }
      X_lincom_use <- X_lincom[lincom_row_index, , drop = FALSE]
      y_lincom_use <- y_lincom[lincom_row_index]
      Transform_lincom_use <- X_fe_lincom[lincom_row_index, , drop = FALSE]
      sigma_i_lincom_use <- as.vector(sigma_i)[lincom_row_index]
    } else {
      X_lincom_use <- X_lincom
      y_lincom_use <- y_lincom
      Transform_lincom_use <- X_fe_lincom
      sigma_i_lincom_use <- as.vector(sigma_i)
    }
    .kss_progress(progress, "leave_out_KSS: running lincom projection")
    lincom_result <- lincom_KSS(
      y = y_lincom_use,
      X = X_lincom_use,
      Z = Z,
      Transform = Transform_lincom_use,
      sigma_i = sigma_i_lincom_use,
      labels = labels_lincom
    )
  }

  elapsed_time <- proc.time()[["elapsed"]] - start_time
  result <- structure(
    list(
      call = match.call(),
      method = "leave_out_KSS",
      settings = list(
        leave_out_level = leave_out_level,
        type_algorithm = type_algorithm,
        simulations_JLA = simulations_JLA,
        lincom_do = lincom_do,
        parallel = paral,
        seed = Cd
      ),
      sample_info = list(
        leave_one_out_connected_set = leave_one_out_stats,
        control_adjusted = control_adjusted_stats
      ),
      estimates = list(
        biased = biased_estimates,
        bias_corrected = bias_corrected_estimates,
        table = .kss_estimates_table(biased_estimates, bias_corrected_estimates)
      ),
      effects = output_dt,
      lincom = lincom_result,
      elapsed_seconds = elapsed_time
    ),
    class = "leave_out_kss_result"
  )

  if (!is.null(csv_file)) {
    .kss_progress(progress, sprintf("leave_out_KSS: writing effects to %s", csv_file))
    .kss_write_csv(output_dt, csv_file)
  }
  if (!is.null(txt_file)) {
    .kss_progress(progress, sprintf("leave_out_KSS: writing summary to %s", txt_file))
    .kss_write_text(.kss_format_estimate_lines(result), txt_file)
  }

  result
}
