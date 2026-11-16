library(testthat)
library(mcc)

# 1. Define a robust dataset at the top including a grouping factor
mydata <- data.frame(
  id     = c(1, 2, 3, 4, 4, 4, 5, 5),
  time   = c(8, 1, 5, 2, 6, 7, 3, 4),
  status = c(0, 0, 2, 1, 1, 1, 1, 2),
  tstart = c(0, 0, 0, 0, 2, 6, 0, 3),
  group  = c("A", "A", "B", "B", "B", "B", "A", "A")
)

# ------------------------------------------------------------------------------
# TEST 1: Without Left-Truncation
# ------------------------------------------------------------------------------
test_that("mcc calculates accurately WITHOUT left-truncation", {
  result_no_trunc <- mcc(id = mydata$id, time = mydata$time,
                         status = mydata$status, Tstart = 0)

  expected_time <- c(1, 2, 3, 4, 5, 6, 7, 8)
  expected_mcc  <- c(0.00, 0.25, 0.50, 0.50, 0.50, 0.75, 1.00, 1.00)

  res_df <- as.data.frame(result_no_trunc)
  expect_equal(res_df$time, expected_time)
  expect_equal(res_df$MCC, expected_mcc, tolerance = 1e-5)
})

# ------------------------------------------------------------------------------
# TEST 2: With Left-Truncation
# ------------------------------------------------------------------------------
test_that("mcc calculates accurately WITH left-truncation", {
  result_trunc <- mcc(id = mydata$id, time = mydata$time,
                      status = mydata$status, Tstart = mydata$tstart)

  expected_time <- c(1, 2, 3, 4, 5, 6, 7, 8)
  expected_mcc  <- c(0.00, 0.25, 0.50, 0.50, 0.50, 0.75, 0.9166667, 0.9166667)

  res_df <- as.data.frame(result_trunc)
  expect_equal(res_df$time, expected_time)
  expect_equal(res_df$MCC, expected_mcc, tolerance = 1e-5)
})

# ------------------------------------------------------------------------------
# TEST 3: Group Stratification Mechanics
# ------------------------------------------------------------------------------
test_that("mcc stratifies calculations correctly across groups", {
  result_grp <- mcc(id = mydata$id, time = mydata$time,
                    status = mydata$status, Tstart = mydata$tstart,
                    group = mydata$group)

  res_df <- as.data.frame(result_grp)

  # Ensure the group column exists and contains exactly the levels we supplied
  expect_true("group" %in% names(res_df))
  expect_setequal(res_df$group, c("A", "B"))

  # Check group-specific timeline subsets to ensure independent evaluation grids
  df_A <- subset(res_df, group == "A")
  df_B <- subset(res_df, group == "B")

  # Group A event times: 1, 3, 4, 8; Group B event times: 2, 5, 6, 7
  expect_contains(df_A$time, c(1, 3, 4, 8))
  expect_contains(df_B$time, c(2, 5, 6, 7))
})

# ------------------------------------------------------------------------------
# TEST 4: Bootstrap Confidence Interval Components and Wide Summary Formats
# ------------------------------------------------------------------------------
test_that("mcc bootstrap components compile with correct dimensions", {
  # Keep niter tiny (e.g., 5) just to verify structure speeds up automated checks
  result_ci <- mcc(id = mydata$id, time = mydata$time,
                   status = mydata$status, Tstart = mydata$tstart,
                   group = mydata$group, # Explicitly include group to test wide summary conversion
                   ci = TRUE, conf.level = 0.90, niter = 5, seed = 123)

  # 1. Verify object list metadata parameters
  expect_true(result_ci$ci)
  expect_equal(result_ci$conf.level, 0.90)

  # 2. Verify raw data frame extraction retains unpivoted numeric bounds
  res_df <- as.data.frame(result_ci)
  expect_true("lci" %in% names(res_df))
  expect_true("uci" %in% names(res_df))
  expect_true(all(res_df$uci >= res_df$lci, na.rm = TRUE))

  # 3. Verify that summary() cleanly outputs your preferred dynamic wide-format columns
  sum_df <- summary(result_ci, times = c(2, 5))

  # The column headers must match your dynamic naming pattern: MCC_Group (90% CI)
  expected_headers <- c("time", "MCC_A (90% CI)", "MCC_B (90% CI)")
  expect_setequal(colnames(sum_df), expected_headers)

  # Verify that the cell values are formatted cleanly as strings: "MCC (LCI, UCI)"
  expect_type(sum_df$`MCC_A (90% CI)`, "character")
})

# ------------------------------------------------------------------------------
# TEST 5: Input Sanitization and Assertions (Fail-Fast Checks)
# ------------------------------------------------------------------------------
test_that("mcc catches invalid user inputs gracefully", {

  # Status vectors must strictly contain 0, 1, or 2
  expect_error(
    mcc(id = 1:3, time = c(1,2,3), status = c(1, 2, 4)),
    "'status' must only contain values 0 \\(censored\\), 1 \\(event\\), or 2 \\(competing risk\\)\\."
  )

  # Vector lengths must be matching
  expect_error(
    mcc(id = 1:4, time = c(1,2,3), status = c(1,0,1)),
    "Lengths of 'id', 'time', 'status', and 'Tstart' must be identical\\."
  )

  # Missing grouping indices must trigger a hard stop rather than crashing downstream
  bad_group <- c("A", "A", "B", NA, "B", "B", "A", "A")
  expect_error(
    mcc(id = mydata$id, time = mydata$time, status = mydata$status,
        Tstart = mydata$tstart, group = bad_group),
    "'group' contains missing values \\(NA\\)\\."
  )
})

# ------------------------------------------------------------------------------
# TEST 6: S3 Print Method Execution (Console Output Checks)
# ------------------------------------------------------------------------------
test_that("S3 print and summary methods execute and stream console output without errors", {

  # Generate standard unstratified and stratified objects
  result_unstd <- mcc(id = mydata$id, time = mydata$time, status = mydata$status, Tstart = 0)
  result_grp   <- mcc(id = mydata$id, time = mydata$time, status = mydata$status,
                      Tstart = mydata$tstart, group = mydata$group)

  # 1. Test print.mcc() console output for unstratified runs
  expect_output(print(result_unstd), "Mean Cumulative Count \\(MCC\\) Estimates")
  expect_output(print(result_unstd), "Overall Cohort \\(Unstratified Analysis\\)")

  # 2. Test print.mcc() console output for stratified runs
  expect_output(print(result_grp), "Stratified Analysis")
  expect_output(print(result_grp), "--- Group: A ---")
  expect_output(print(result_grp), "Showing the first 3 estimates for each group:")

  # 3. Test print.summary.mcc() console output
  sum_obj <- summary(result_grp, times = c(2, 5))
  expect_output(print(sum_obj), "Mean Cumulative Count \\(MCC\\) Summary")
})

# ------------------------------------------------------------------------------
# TEST 7: Plotting unstratified MCC without bootstrap confidence intervals
# ------------------------------------------------------------------------------
test_that("plotting handles unstratified MCC without confidence intervals", {
  fit <- mcc(id = mydata$id, time = mydata$time,
             status = mydata$status, Tstart = 0)

  expect_message(
    p <- plot_mcc_with_risk_table(fit, risk_times = c(2, 4, 8)),
    "Confidence intervals are unavailable"
  )
  expect_s3_class(p$plot, "patchwork")
  expect_s3_class(p$p_mcc, "ggplot")
  expect_s3_class(p$p_risk, "ggplot")
  expect_equal(p$risk_times, c(2, 4, 8))
  expect_equal(unique(p$plot_data$plot_group), "Unspecified")
  expect_equal(nrow(p$risk_table), 3L)
  expect_equal(p$risk_table$n_risk, c(4L, 4L, 1L))
  expect_silent(ggplot2::ggplot_build(p$p_mcc))
  expect_silent(ggplot2::ggplot_build(p$p_risk))
})

# ------------------------------------------------------------------------------
# TEST 8: Stratified plots and left-truncation-aware risk counts
# ------------------------------------------------------------------------------
test_that("plotting counts distinct subjects in left-truncated risk sets", {
  fit <- mcc(id = mydata$id, time = mydata$time,
             status = mydata$status, Tstart = mydata$tstart,
             group = mydata$group)

  p <- suppressMessages(plot_mcc_with_risk_table(
    fit, risk_times = c(3, 4, 6), conf_int = FALSE,
    control_group = "A", color_values = c(A = "blue", B = "red")
  ))

  expect_s3_class(p$plot, "patchwork")
  expect_setequal(unique(p$plot_data$plot_group), c("A", "B"))
  risk <- p$risk_table
  count_at <- function(g, t) risk$n_risk[risk$group == g & risk$time == t]
  expect_equal(count_at("A", 3), 2L)
  expect_equal(count_at("B", 3), 2L)
  expect_equal(count_at("A", 4), 2L)
  expect_equal(count_at("B", 4), 2L)
  expect_equal(count_at("A", 6), 1L)
  expect_equal(count_at("B", 6), 1L)
  expect_true("A (Control)" %in% as.character(p$plot_data$legend_group))
  expect_silent(ggplot2::ggplot_build(p$p_mcc))
  expect_silent(ggplot2::ggplot_build(p$p_risk))
})

# ------------------------------------------------------------------------------
# TEST 9: Filtering, standard risk sets, and standalone data-frame interface
# ------------------------------------------------------------------------------
test_that("plotting supports group selection and no left truncation", {
  fit <- mcc(id = mydata$id, time = mydata$time,
             status = mydata$status, Tstart = 0,
             group = mydata$group)
  p <- plot_mcc_with_risk_table(
    as.data.frame(fit), data = fit$data, Tstart = NULL,
    group_value = "A", risk_times = c(3, 8), conf_int = FALSE
  )
  expect_equal(unique(p$plot_data$plot_group), "A")
  expect_equal(unique(p$risk_table$group), "A")
  expect_equal(p$risk_table$n_risk, c(2L, 1L))
  expect_error(
    plot_mcc_with_risk_table(as.data.frame(fit),
                             group_value = "A", conf_int = FALSE),
    "data.*must be supplied"
  )
  expect_error(
    plot_mcc_with_risk_table(fit, risk_times = c(NA_real_, 3),
                             conf_int = FALSE),
    "risk_times"
  )
  expect_error(
    plot_mcc_with_risk_table(fit, group_value = "Nonexistent",
                             conf_int = FALSE),
    "None of the requested groups"
  )
})

# ------------------------------------------------------------------------------
# TEST 10: Confidence ribbons and S3 plot method
# ------------------------------------------------------------------------------
test_that("plotting supports confidence intervals and plot.mcc dispatch", {
  fit <- mcc(id = mydata$id, time = mydata$time,
             status = mydata$status, Tstart = 0,
             group = mydata$group)
  # Exercise rendering without running an additional expensive bootstrap.
  fit$estimates$lci <- pmax(0, fit$estimates$MCC - 0.1)
  fit$estimates$uci <- fit$estimates$MCC + 0.1
  fit$ci <- TRUE

  p <- plot_mcc_with_risk_table(fit, conf_int = TRUE,
                                risk_times = c(2, 4, 8))
  expect_s3_class(p$plot, "patchwork")
  expect_true(all(c("lci", "uci") %in% names(p$plot_data)))
  expect_silent(ggplot2::ggplot_build(p$p_mcc))

  # plot.mcc() prints a patchwork object; use a temporary PDF device in tests.
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  expect_s3_class(plot(fit, conf_int = TRUE), "patchwork")
})
