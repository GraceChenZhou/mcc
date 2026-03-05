# 1. Define the input data once at the top so both tests can access it
mydata <- data.frame(
  id = c(1, 2, 3, 4, 4, 4, 5, 5),
  time = c(8, 1, 5, 2, 6, 7, 3, 4),
  status = c(0, 0, 2, 1, 1, 1, 1, 2),
  tstart = c(0, 0, 0, 0, 2, 6, 0, 3)
)

# ------------------------------------------------------------------------------
# TEST 1: Without Left-Truncation
# ------------------------------------------------------------------------------
test_that("mcc calculates accurately WITHOUT left-truncation", {

  # Run package function
  result_no_trunc <- mcc(id = mydata$id, time = mydata$time,
                         status = mydata$status, Tstart = 0)

  # Define the KNOWN ACCURATE results to compare against
  expected_time <- c(1, 2, 3, 4, 5, 6, 7, 8)
  expected_mcc <- c(0.00, 0.25, 0.50, 0.50, 0.50, 0.75, 1.00, 1.00)

  # Run the strict comparisons
  expect_equal(result_no_trunc$time, expected_time)
  expect_equal(result_no_trunc$MCC, expected_mcc, tolerance = 1e-5)

})

# ------------------------------------------------------------------------------
# TEST 2: With Left-Truncation
# ------------------------------------------------------------------------------
test_that("mcc calculates accurately WITH left-truncation", {

  # Run package function
  result_trunc <- mcc(id = mydata$id, time = mydata$time,
                      status = mydata$status, Tstart = mydata$tstart)

  # Define the KNOWN ACCURATE results to compare against
  expected_time <- c(1, 2, 3, 4, 5, 6, 7, 8)
  expected_mcc <- c(0.00, 0.25, 0.50, 0.50, 0.50, 0.75, 0.9166667, 0.9166667)

  # Run the strict comparisons
  expect_equal(result_trunc$time, expected_time)
  expect_equal(result_trunc$MCC, expected_mcc, tolerance = 1e-5)

})
