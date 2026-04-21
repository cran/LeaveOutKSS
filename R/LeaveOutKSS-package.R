#' Leave-Out Variance Component Estimation for Two-Way Fixed Effects Models
#'
#' `LeaveOutKSS` packages an 'R' translation of the original 'MATLAB' package
#' of Kline, Saggio, and Solvsten (2020) for leave-out bias correction and
#' variance decomposition in two-way fixed effects models of the Abowd,
#' Kramarz, and Margolis (1999; AKM) type.
#'
#' The package mirrors the logic of the original script-based implementation in
#' this repository while exposing object-returning, side-effect-free-by-default
#' workflows. Core estimation functions return structured results and write
#' files only when users explicitly provide output paths. The main user-facing
#' workflows are:
#'
#' - [leave_out_KSS()] for leave-out bias-corrected decomposition in two-way
#'   fixed effects models.
#' - [leave_out_KSS_fe()] for the same decomposition when some controls are
#'   supplied as categorical variables to be expanded internally.
#' - [rsquared_comp()] for comparing the fit of a two-way fixed effects model
#'   with a saturated worker-firm specification.
#' - [fast_fe_est()] for recovering fitted values and adjusted outcomes from
#'   one-way or two-way fixed effects models.
#'
#' The implementation follows the structure of the original 'MATLAB' package
#' and the accompanying vignette material:
#'
#' 1. Construct the largest connected set of firms.
#' 2. Prune the sample to a leave-one-worker-out connected set.
#' 3. Partial out controls when requested.
#' 4. Compute statistical leverages exactly or by Johnson-Lindenstrauss
#'    approximation (JLA).
#' 5. Form plug-in and leave-out bias-corrected variance component estimates.
#'
#' A small matched worker-firm panel used by the examples is bundled at
#' `system.file("extdata", "test.csv", package = "LeaveOutKSS")`.
#'
#' @references
#' Abowd, J. M., Kramarz, F., and Margolis, D. N. (1999). High wage workers and
#' high wage firms. *Econometrica*, 67(2), 251-333.
#'
#' Kline, P., Saggio, R., and Solvsten, M. (2020). Leave-out estimation of
#' variance components. *Econometrica*, 88(5), 1859-1898.
#'
#' Johnson, W. B., and Lindenstrauss, J. (1984). Extensions of Lipschitz
#' mappings into a Hilbert space. In *Conference in Modern Analysis and
#' Probability*, 189-206.
#'
#' @keywords internal
#' @import data.table
#' @import Matrix
#' @import igraph
#' @import sanic
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom stats ave rbinom runif var
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console
"_PACKAGE"

utils::globalVariables(
  c(
    "..absorb_col", ".", "count", "firmid_old", "first_obs", "gcs", "i",
    "id_old", "int_id", "interact", "j", "lagfirmid", "match_id", "movers",
    "order_id", "result", "residuals", "stayer", "stayer_sum", "stayers",
    "total", "u", "v", "w", "x", "y_adj", "y_hat"
  )
)
