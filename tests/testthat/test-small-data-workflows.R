test_that("small-data leave_out_KSS workflow returns an object and writes files only when requested", {
  skip_on_cran()

  path <- system.file("extdata", "test.csv", package = "LeaveOutKSS")
  dt <- data.table::fread(path, header = FALSE)
  out <- tempfile("leaveoutkss_test_")
  on.exit(unlink(c(paste0(out, ".csv"), paste0(out, ".txt"))), add = TRUE)

  res <- leave_out_KSS(
    y = dt[[4]],
    id = dt[[1]],
    firmid = dt[[2]],
    simulations_JLA = 5,
    paral = FALSE,
    progress = FALSE
  )

  expect_s3_class(res, "leave_out_kss_result")
  expect_true(is.data.frame(res$estimates$table))
  expect_true(is.data.frame(res$effects))
  expect_false(file.exists(paste0(out, ".csv")))
  expect_false(file.exists(paste0(out, ".txt")))

  res_written <- leave_out_KSS(
    y = dt[[4]],
    id = dt[[1]],
    firmid = dt[[2]],
    simulations_JLA = 5,
    paral = FALSE,
    csv_file = paste0(out, ".csv"),
    txt_file = paste0(out, ".txt"),
    progress = FALSE
  )

  expect_s3_class(res_written, "leave_out_kss_result")
  expect_true(file.exists(paste0(out, ".csv")))
  expect_true(file.exists(paste0(out, ".txt")))
})

test_that("small-data leave_out_KSS_fe workflow returns an object and writes files only when requested", {
  skip_on_cran()

  path <- system.file("extdata", "test.csv", package = "LeaveOutKSS")
  dt <- data.table::fread(path, header = FALSE)
  out <- tempfile("leaveoutkss_fe_test_")
  on.exit(unlink(c(paste0(out, ".csv"), paste0(out, ".txt"))), add = TRUE)

  res <- leave_out_KSS_fe(
    y = dt[[4]],
    id = dt[[1]],
    firmid = dt[[2]],
    controls = cbind(year = dt[[3]]),
    absorb_col = 1,
    simulations_JLA = 5,
    paral = FALSE,
    progress = FALSE
  )

  expect_s3_class(res, "leave_out_kss_result")
  expect_true(is.data.frame(res$estimates$table))
  expect_true(is.data.frame(res$effects))
  expect_false(file.exists(paste0(out, ".csv")))
  expect_false(file.exists(paste0(out, ".txt")))

  res_written <- leave_out_KSS_fe(
    y = dt[[4]],
    id = dt[[1]],
    firmid = dt[[2]],
    controls = cbind(year = dt[[3]]),
    absorb_col = 1,
    simulations_JLA = 5,
    paral = FALSE,
    csv_file = paste0(out, ".csv"),
    txt_file = paste0(out, ".txt"),
    progress = FALSE
  )

  expect_s3_class(res_written, "leave_out_kss_result")
  expect_true(file.exists(paste0(out, ".csv")))
  expect_true(file.exists(paste0(out, ".txt")))
})

test_that("small-data rsquared workflow returns an object and writes text only when requested", {
  path <- system.file("extdata", "test.csv", package = "LeaveOutKSS")
  dt <- data.table::fread(path, header = FALSE)
  out <- tempfile("leaveoutkss_r2_test_")
  on.exit(unlink(out), add = TRUE)

  res <- rsquared_comp(
    y = dt[[4]],
    id = dt[[1]],
    firmid = dt[[2]],
    progress = FALSE
  )

  expect_s3_class(res, "rsquared_comp_result")
  expect_true(is.data.frame(res$summary))
  expect_false(file.exists(out))

  res_written <- rsquared_comp(
    y = dt[[4]],
    id = dt[[1]],
    firmid = dt[[2]],
    txt_file = out,
    progress = FALSE
  )

  expect_s3_class(res_written, "rsquared_comp_result")
  expect_true(file.exists(out))
})

test_that("misaligned inputs fail fast with informative errors", {
  path <- system.file("extdata", "test.csv", package = "LeaveOutKSS")
  dt <- data.table::fread(path, header = FALSE)

  expect_error(
    leave_out_KSS(
      y = dt[[4]],
      id = dt[[1]],
      firmid = dt[[2]],
      controls = matrix(1, nrow = nrow(dt) - 1, ncol = 1),
      progress = FALSE
    ),
    "`controls` must have"
  )

  expect_error(
    leave_out_KSS_fe(
      y = dt[[4]],
      id = dt[[1]],
      firmid = dt[[2]],
      controls = cbind(year = dt[[3]]),
      absorb_col = 1.5,
      progress = FALSE
    ),
    "integer-valued"
  )

  expect_error(
    fast_fe_est(
      y = dt[[4]],
      id = dt[[1]][-1],
      progress = FALSE
    ),
    "Input lengths must match"
  )
})
