
library(testthat)

test_that("simulate_subject returns valid recurrent-event data", {
  set.seed(123)

  x <- simulate_subject(
    id = 1,
    group = "A",
    entry_time = 2,
    max_followup = 10,
    event_rate = 0.5
  )

  expect_s3_class(x, "data.frame")
  expect_true(all(c("id", "time", "status", "Tstart", "group") %in% names(x)))
  expect_true(all(x$id == 1))
  expect_true(all(x$Tstart == 2))
  expect_true(all(x$time >= x$Tstart))
  expect_true(all(x$time <= 10))
  expect_true(all(x$status %in% c(0, 1, 2)))

  # Exactly one terminal observation
  expect_equal(sum(x$status %in% c(0, 2)), 1L)
})

test_that("simulate_mcc_data is reproducible", {
  x <- simulate_mcc_data(
    n_per_group = 20,
    seed = 2026
  )

  y <- simulate_mcc_data(
    n_per_group = 20,
    seed = 2026
  )

  expect_equal(x, y)
  expect_equal(length(unique(x$id)), 40L)
  expect_setequal(unique(x$group), c("A", "B"))
})

test_that("simulated subjects have valid terminal records", {
  x <- simulate_mcc_data(
    n_per_group = 30,
    seed = 123
  )

  terminal <- x[x$status %in% c(0, 2), ]

  expect_equal(nrow(terminal), 60L)
  expect_equal(length(unique(terminal$id)), 60L)

  expect_true(all(x$time >= x$Tstart))
  expect_true(all(x$status %in% c(0, 1, 2)))
})

test_that("zero recurrent-event rate produces no recurrent events", {
  x <- simulate_subject(
    id = 1,
    event_rate = 0,
    death_rate = 0.03,
    censor_rate = 0.05
  )

  expect_equal(sum(x$status == 1), 0L)
  expect_equal(nrow(x), 1L)
})

test_that("simulated data can be passed to mcc", {
  x <- simulate_mcc_data(
    n_per_group = 20,
    seed = 123
  )

  fit <- mcc(
    id = x$id,
    time = x$time,
    status = x$status,
    Tstart = x$Tstart,
    group = x$group,
    ci = FALSE
  )

  expect_s3_class(fit, "mcc")
  expect_true(all(c("time", "MCC", "group") %in%
                    names(fit$estimates)))
  expect_setequal(unique(fit$estimates$group), c("A", "B"))
})
