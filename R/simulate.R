#' Simulate One Subject's Recurrent-Event History
#'
#' Simulates recurrent events from a homogeneous Poisson process after study
#' entry, together with an independent competing event and right censoring.
#' The output has the column structure expected by [mcc()]. This is an
#' illustrative data generator, not a model of disease-specific event hazards.
#'
#' @param id A single non-missing subject identifier.
#' @param group A single non-missing group label (default `"A"`).
#' @param entry_time Nonnegative entry time, strictly less than `max_followup`.
#' @param max_followup Positive administrative study end time on the same scale
#'   as `entry_time`.
#' @param event_rate Nonnegative recurrent-event rate per unit time.
#' @param death_rate Nonnegative competing-event hazard per unit time.
#' @param censor_rate Nonnegative independent right-censoring hazard per unit time.
#'
#' @return A data frame with `id`, `time`, `status`, `Tstart`, and `group`.
#'   Each recurrent event has `status = 1`; exactly one final row has
#'   `status = 2` (competing event) or `status = 0` (censoring). `Tstart`
#'   is the subject's fixed study-entry time on every row.
#' @export
#' @examples
#' set.seed(2026)
#' one <- simulate_subject(id = 1, group = "A", entry_time = 0.5)
#' stopifnot(tail(one$status, 1) %in% c(0L, 2L))
simulate_subject <- function(id, group = "A", entry_time = 0,
                             max_followup = 10, event_rate = 0.20,
                             death_rate = 0.03, censor_rate = 0.05) {
  scalar_number <- function(x, name, lower = 0, strict = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        (if (strict) x <= lower else x < lower)) {
      stop("'", name, "' must be a single finite number ",
           if (strict) "greater than " else "at least ", lower, ".",
           call. = FALSE)
    }
  }
  if (length(id) != 1L || is.na(id) || !is.atomic(id)) {
    stop("'id' must be a single non-missing atomic value.", call. = FALSE)
  }
  if (length(group) != 1L || is.na(group) || !is.atomic(group)) {
    stop("'group' must be a single non-missing atomic value.", call. = FALSE)
  }
  scalar_number(entry_time, "entry_time")
  scalar_number(max_followup, "max_followup", strict = TRUE)
  scalar_number(event_rate, "event_rate")
  scalar_number(death_rate, "death_rate")
  scalar_number(censor_rate, "censor_rate")
  if (entry_time >= max_followup) {
    stop("'entry_time' must be less than 'max_followup'.", call. = FALSE)
  }

  # Zero hazard means the corresponding terminal event cannot occur.
  death_time <- if (death_rate == 0) Inf else entry_time + stats::rexp(1L, death_rate)
  censor_time <- if (censor_rate == 0) Inf else entry_time + stats::rexp(1L, censor_rate)
  exit_time <- min(death_time, censor_time, max_followup)
  terminal_status <- if (death_time < min(censor_time, max_followup)) 2L else 0L

  # Conditional on the observed duration, event times are uniform order
  # statistics for a homogeneous Poisson process.
  n_events <- stats::rpois(1L, event_rate * (exit_time - entry_time))
  event_times <- if (n_events) {
    sort(stats::runif(n_events, min = entry_time, max = exit_time))
  } else numeric(0)

  data.frame(
    id = rep(id, n_events + 1L),
    time = c(event_times, exit_time),
    status = c(rep(1L, n_events), terminal_status),
    Tstart = rep(entry_time, n_events + 1L),
    group = rep(as.character(group), n_events + 1L),
    stringsAsFactors = FALSE
  )
}

#' Simulate Recurrent-Event Data for Multiple Groups
#'
#' Generates independent subjects using [simulate_subject()], with group-
#' specific recurrent-event rates and optional delayed entry. Intended for
#' examples, regression tests, and simulation demonstrations.
#'
#' @param n_per_group Positive integer number of subjects in each group.
#' @param event_rates Named numeric vector of nonnegative recurrent-event rates,
#'   e.g. `c(A = 0.22, B = 0.38)`.
#' @param max_followup Positive administrative end time.
#' @param entry_range Length-two numeric vector giving the minimum and maximum
#'   entry time; both must be nonnegative and strictly below `max_followup`.
#' @param death_rate,censor_rate Nonnegative hazards shared across groups.
#' @param seed Optional integer random seed. If supplied, the previous R random
#'   number generator state is restored on exit.
#'
#' @return A data frame with `id`, `time`, `status`, `Tstart`, and `group`.
#'   IDs are unique across groups; each subject has one terminal record.
#' @export
#' @examples
#' dat <- simulate_mcc_data(n_per_group = 10,
#'                          event_rates = c(A = 0.22, B = 0.38), seed = 2026)
#' stopifnot(length(unique(dat$id)) == 20L)
simulate_mcc_data <- function(n_per_group = 150L,
                              event_rates = c(A = 0.22, B = 0.38),
                              max_followup = 10, entry_range = c(0, 2),
                              death_rate = 0.03, censor_rate = 0.05,
                              seed = NULL) {
  if (!is.numeric(n_per_group) || length(n_per_group) != 1L ||
      !is.finite(n_per_group) || n_per_group < 1 ||
      n_per_group != floor(n_per_group)) {
    stop("'n_per_group' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(event_rates) || !length(event_rates) ||
      anyNA(event_rates) || any(!is.finite(event_rates)) || any(event_rates < 0) ||
      is.null(names(event_rates)) || anyNA(names(event_rates)) ||
      any(!nzchar(names(event_rates))) || anyDuplicated(names(event_rates))) {
    stop("'event_rates' must be a non-empty named vector of finite nonnegative rates with unique group names.", call. = FALSE)
  }
  if (!is.numeric(max_followup) || length(max_followup) != 1L ||
      !is.finite(max_followup) || max_followup <= 0) {
    stop("'max_followup' must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(entry_range) || length(entry_range) != 2L ||
      anyNA(entry_range) || any(!is.finite(entry_range)) ||
      entry_range[1] < 0 || entry_range[2] < entry_range[1] ||
      entry_range[2] >= max_followup) {
    stop("'entry_range' must be two ordered nonnegative times below 'max_followup'.", call. = FALSE)
  }
  for (nm in c("death_rate", "censor_rate")) {
    x <- get(nm)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
      stop("'", nm, "' must be a finite nonnegative number.", call. = FALSE)
    }
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
      stop("'seed' must be a nonnegative integer or NULL.", call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  groups <- rep(names(event_rates), each = as.integer(n_per_group))
  entry_times <- stats::runif(length(groups), entry_range[1], entry_range[2])
  rows <- lapply(seq_along(groups), function(i) {
    simulate_subject(
      id = i, group = groups[i], entry_time = entry_times[i],
      max_followup = max_followup,
      event_rate = unname(event_rates[groups[i]]),
      death_rate = death_rate, censor_rate = censor_rate
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
