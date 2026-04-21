#' Internal Helpers for Progress, Formatting, and Optional Export
#'
#' These helpers keep the computational routines side-effect free by default
#' while providing opt-in progress reporting and file export.
#'
#' @name kss_internal_helpers
#' @keywords internal
NULL

#' @keywords internal
.kss_progress <- function(progress, ..., .sep = "") {
  msg <- paste(..., sep = .sep)
  if (isTRUE(progress)) {
    message(msg)
  }
  invisible(msg)
}

#' @keywords internal
.kss_validate_optional_path <- function(path, arg) {
  if (is.null(path)) {
    return(NULL)
  }

  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop(sprintf("`%s` must be NULL or a non-empty character string.", arg), call. = FALSE)
  }

  path
}

#' @keywords internal
.kss_validate_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
  }
  x
}

#' @keywords internal
.kss_validate_scalar_number <- function(x, arg, integerish = FALSE, min_value = -Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("`%s` must be a finite numeric scalar.", arg), call. = FALSE)
  }

  if (integerish && x != as.integer(x)) {
    stop(sprintf("`%s` must be an integer-valued scalar.", arg), call. = FALSE)
  }

  if (x < min_value) {
    stop(
      sprintf("`%s` must be greater than or equal to %s.", arg, min_value),
      call. = FALSE
    )
  }

  x
}

#' @keywords internal
.kss_assert_same_length <- function(..., .arg_names = NULL) {
  args <- list(...)
  if (length(args) <= 1L) {
    return(invisible(NULL))
  }

  arg_names <- .arg_names
  if (is.null(arg_names)) {
    arg_names <- names(args)
  }
  if (is.null(arg_names) || any(!nzchar(arg_names))) {
    arg_names <- paste0("arg", seq_along(args))
  }

  lens <- vapply(args, length, integer(1))
  if (length(unique(lens)) != 1L) {
    stop(
      sprintf(
        "Input lengths must match. Received: %s.",
        paste(sprintf("%s = %s", arg_names, lens), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

#' @keywords internal
.kss_validate_row_aligned_input <- function(x, n, arg) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.vector(x) || is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (is.null(dim(x))) {
    stop(sprintf("`%s` must be a vector, matrix, or data frame.", arg), call. = FALSE)
  }

  if (nrow(x) != n) {
    stop(
      sprintf("`%s` must have %s rows to match the main inputs.", arg, n),
      call. = FALSE
    )
  }

  x
}

#' @keywords internal
.kss_stop_if_missing_or_nonfinite <- function(x, arg) {
  if (anyNA(x)) {
    stop(sprintf("`%s` must not contain missing values.", arg), call. = FALSE)
  }

  if (is.numeric(x) && any(!is.finite(x))) {
    stop(sprintf("`%s` must not contain non-finite values.", arg), call. = FALSE)
  }

  invisible(NULL)
}

#' @keywords internal
.kss_stop_if_no_rows <- function(n, context) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1L) {
    stop(sprintf("%s did not retain any observations.", context), call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
.kss_prepare_output_path <- function(path) {
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  path
}

#' @keywords internal
.kss_estimates_table <- function(biased, bias_corrected) {
  components <- union(names(biased), names(bias_corrected))
  data.frame(
    component = components,
    biased = unname(biased[components]),
    bias_corrected = unname(bias_corrected[components]),
    check.names = FALSE
  )
}

#' @keywords internal
.kss_format_number <- function(x) {
  formatC(x, digits = 6, format = "fg", flag = "#")
}

#' @keywords internal
.kss_format_estimate_lines <- function(result) {
  estimate_tbl <- result$estimates$table
  lines <- c(
    "LeaveOutKSS Decomposition Estimates",
    sprintf("Method: %s", result$method),
    ""
  )

  for (i in seq_len(nrow(estimate_tbl))) {
    lines <- c(
      lines,
      sprintf(
        "%s | biased: %s | bias corrected: %s",
        estimate_tbl$component[[i]],
        .kss_format_number(estimate_tbl$biased[[i]]),
        .kss_format_number(estimate_tbl$bias_corrected[[i]])
      )
    )
  }

  if (!is.null(result$sample_info$leave_one_out_connected_set)) {
    sample_info <- result$sample_info$leave_one_out_connected_set
    lines <- c(
      lines,
      "",
      "Leave-One-Out Connected Set",
      sprintf("Mean outcome: %s", .kss_format_number(sample_info$mean_outcome)),
      sprintf("Variance outcome: %s", .kss_format_number(sample_info$variance_outcome)),
      sprintf("Movers: %s", sample_info$n_movers),
      sprintf("Firms: %s", sample_info$n_firms),
      sprintf("Observations: %s", sample_info$n_observations)
    )
  }

  if (!is.null(result$lincom) && nrow(result$lincom$results) > 0L) {
    lines <- c(lines, "", "Lincom Results")
    for (i in seq_len(nrow(result$lincom$results))) {
      row <- result$lincom$results[i, ]
      lines <- c(
        lines,
        sprintf(
          "%s | estimate: %s | white_se: %s | kss_se: %s | t_statistic: %s",
          row$term,
          .kss_format_number(row$estimate),
          .kss_format_number(row$white_se),
          .kss_format_number(row$kss_se),
          .kss_format_number(row$t_statistic)
        )
      )
    }
  }

  lines
}

#' @keywords internal
.kss_format_rsquared_lines <- function(result) {
  lines <- c("LeaveOutKSS R-Squared Comparison", "")

  for (i in seq_len(nrow(result$summary))) {
    row <- result$summary[i, ]
    lines <- c(
      lines,
      sprintf(
        "%s | r_squared: %s | adjusted_r_squared: %s | n_parameters: %s",
        row$model,
        .kss_format_number(row$r_squared),
        .kss_format_number(row$adjusted_r_squared),
        row$n_parameters
      )
    )
  }

  lines
}

#' @keywords internal
.kss_write_text <- function(lines, path) {
  writeLines(lines, .kss_prepare_output_path(path), useBytes = TRUE)
  invisible(path)
}

#' @keywords internal
.kss_write_csv <- function(x, path) {
  data.table::fwrite(x, .kss_prepare_output_path(path), row.names = FALSE, quote = FALSE)
  invisible(path)
}

#' Print a LeaveOutKSS Decomposition Result
#'
#' @param x A result returned by [leave_out_KSS()] or [leave_out_KSS_fe()].
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @method print leave_out_kss_result
#' @export
print.leave_out_kss_result <- function(x, ...) {
  cat("LeaveOutKSS decomposition result\n")
  print(x$estimates$table, row.names = FALSE)
  invisible(x)
}

#' Summarize a LeaveOutKSS Decomposition Result
#'
#' @param object A result returned by [leave_out_KSS()] or [leave_out_KSS_fe()].
#' @param ... Unused.
#'
#' @return `object`, invisibly.
#' @method summary leave_out_kss_result
#' @export
summary.leave_out_kss_result <- function(object, ...) {
  print(object)
  invisible(object)
}

#' Print a Lincom Result
#'
#' @param x A result returned by [lincom_KSS()].
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @method print lincom_kss_result
#' @export
print.lincom_kss_result <- function(x, ...) {
  cat("LeaveOutKSS lincom result\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' Print an R-Squared Comparison Result
#'
#' @param x A result returned by [rsquared_comp()].
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @method print rsquared_comp_result
#' @export
print.rsquared_comp_result <- function(x, ...) {
  cat("LeaveOutKSS R-squared comparison\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' Print a Fixed Effects Fit Result
#'
#' @param x A result returned by [fast_fe_est()].
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @method print fast_fe_est_result
#' @export
print.fast_fe_est_result <- function(x, ...) {
  cat("LeaveOutKSS fixed-effects fit result\n")
  cat(sprintf("Observations: %s\n", nrow(x$fitted_values)))
  print(utils::head(x$fitted_values))
  invisible(x)
}
