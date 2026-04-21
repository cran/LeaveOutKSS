## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
path <- system.file("extdata", "test.csv", package = "LeaveOutKSS")
dt <- data.table::fread(path, header = FALSE)
data.table::setorder(dt, V1, V3)
dim(dt)
head(dt)

## ----eval = FALSE-------------------------------------------------------------
# res <- leave_out_KSS(
#   y = dt[[4]],
#   id = dt[[1]],
#   firmid = dt[[2]],
#   leave_out_level = "matches",
#   type_algorithm = "JLA",
#   simulations_JLA = 5,
#   paral = FALSE,
#   progress = FALSE
# )
# 
# print(res)
# res$estimates$table

## ----eval = FALSE-------------------------------------------------------------
# stem <- tempfile("leaveoutkss_")
# 
# leave_out_KSS(
#   y = dt[[4]],
#   id = dt[[1]],
#   firmid = dt[[2]],
#   simulations_JLA = 5,
#   paral = FALSE,
#   csv_file = paste0(stem, ".csv"),
#   txt_file = paste0(stem, ".txt"),
#   progress = FALSE
# )
# 
# unlink(paste0(stem, c(".csv", ".txt")))

## ----eval = FALSE-------------------------------------------------------------
# controls <- model.matrix(~ factor(dt[[3]]) - 1)
# controls <- controls[, -ncol(controls), drop = FALSE]
# 
# leave_out_KSS(
#   y = dt[[4]],
#   id = dt[[1]],
#   firmid = dt[[2]],
#   controls = controls,
#   simulations_JLA = 5,
#   paral = FALSE,
#   progress = FALSE
# )

## ----eval = FALSE-------------------------------------------------------------
# leave_out_KSS_fe(
#   y = dt[[4]],
#   id = dt[[1]],
#   firmid = dt[[2]],
#   controls = cbind(year = dt[[3]]),
#   absorb_col = 1,
#   simulations_JLA = 5,
#   paral = FALSE,
#   progress = FALSE
# )

## ----eval = FALSE-------------------------------------------------------------
# leave_out_KSS(
#   y = dt[[4]],
#   id = dt[[1]],
#   firmid = dt[[2]],
#   leave_out_level = "obs",
#   simulations_JLA = 5,
#   paral = FALSE,
#   progress = FALSE
# )

## ----eval = FALSE-------------------------------------------------------------
# region_dummy <- as.numeric(dt[[3]] <= median(dt[[3]], na.rm = TRUE))
# 
# leave_out_KSS(
#   y = dt[[4]],
#   id = dt[[1]],
#   firmid = dt[[2]],
#   simulations_JLA = 5,
#   paral = FALSE,
#   lincom_do = 1,
#   Z_lincom = region_dummy,
#   labels_lincom = list("Early-Year Indicator"),
#   progress = FALSE
# )

## ----eval = FALSE-------------------------------------------------------------
# rsquared_comp(
#   y = dt[[4]],
#   id = dt[[1]],
#   firmid = dt[[2]],
#   progress = FALSE
# )

