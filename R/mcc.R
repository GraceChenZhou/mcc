#' Estimate Mean Cumulative Count for Left-Truncated Competing Risks Data
#'
#' @description
#' Estimates the total burden of recurrent events within a population using the Mean
#' Cumulative Count (MCC) method. Specifically designed to handle left-truncated data
#' in the presence of competing risks by applying time-dependent censoring
#' and truncation weights. Optionally stratifies results by a grouping variable.
#'
#' @param id A vector uniquely identifying individual subjects.
#' @param time A numeric vector of event or censoring times.
#' @param status A numeric vector indicating event status (0 = censored, 1 = event of interest, 2 = competing risk).
#' @param Tstart A numeric vector representing the age or time at study entry (left-truncation).
#'   Defaults to 0 if not provided, assuming subjects are followed from time zero.
#' @param group An optional vector of group labels (e.g., treatment arm, sex) with the same length as \code{id}.
#'   MCC is estimated separately within each group. If \code{NULL} (default), MCC is estimated
#'   for the full dataset.
#' @param ci Logical; if \code{TRUE}, calculates bootstrap confidence intervals. Default is \code{FALSE}.
#' @param conf.level A numeric value between 0 and 1 specifying the confidence level for the
#'   bootstrap intervals. Default is 0.95 (95\% CI).
#' @param niter Integer; the number of bootstrap iterations to run if \code{ci = TRUE}. Default is 1000.
#' @param seed Optional integer to set the seed for reproducible bootstrap results.
#'
#' @return An object of class \code{"mcc"}. This object is a named list containing:
#'   \item{estimates}{A data frame of evaluated evaluation time points, estimated MCC, and bootstrap confidence intervals if \code{ci = TRUE}.}
#'   \item{data}{A data frame containing the original sanitized input data vectors used for summary and risk statistics.}
#'   \item{ci}{The logical value indicating whether confidence intervals were calculated.}
#'   \item{conf.level}{The numerical confidence level utilized for interval calculation bounds.}
#'
#' @author Grace Zhou, Department of Biostatistics at St. Jude Children's Research Hospital \email{grace.zhou@@stjude.org}
#' @note Requires \code{crprep2}, a customized internal version of development \code{mstate::crprep}.
#' @export
#'
#' @examples
#' # Create a sample recurrent event dataset
#' mydata <- data.frame(
#'   id     = c(1, 2, 3, 4, 4, 4, 5, 5),
#'   tstop  = c(8, 1, 5, 2, 6, 7, 3, 4),
#'   status = c(0, 0, 2, 1, 1, 1, 1, 2), # 0=censored, 1=event, 2=competing risk
#'   tstart = c(0, 0, 0, 0, 2, 6, 0, 3),
#'   group  = c("A", "A", "B", "B", "B", "B", "A", "A")
#' )
#'
#' # 1. Calculate unstratified MCC
#' result <- mcc(id = mydata$id, time = mydata$tstop,
#'               status = mydata$status, Tstart = mydata$tstart)
#' print(result)
#'
#' # 2. Extract full estimates cleanly as a standard data frame
#' result_df <- as.data.frame(result)
#'
#' # 3. Generate a milestone summary table at specific time points
#' summary(result, times = c(2, 5, 7))
#'
#' # 4. Calculate MCC stratified by treatment group
#' result_grp <- mcc(id = mydata$id, time = mydata$tstop,
#'                   status = mydata$status, Tstart = mydata$tstart,
#'                   group = mydata$group)
#' print(result_grp)
#' summary(result_grp)
#' summary(result_grp, times = c(2, 5, 7))
#' \donttest{
#' # 5. Calculate MCC with bootstrap confidence intervals (run locally)
#' # Note: We wrap this in \donttest to ensure fast package compilation times.
#' result_grp_ci <- mcc(id = mydata$id, time = mydata$tstop,
#'                      status = mydata$status, Tstart = mydata$tstart,
#'                      group = mydata$group, ci = TRUE, niter = 100, seed = 123)
#' print(result_grp_ci)
#' summary(result_grp_ci)
#' }
mcc <- function(id, time, status, Tstart = NULL, group = NULL, ci = FALSE, conf.level = 0.95, niter = 1000, seed = NULL) {

  # Input sanitization and validation
  n_obs <- length(id)

  # Expand Tstart safely if it's NULL or a single placeholder value (like 0)
  if (is.null(Tstart) || (is.numeric(Tstart) && length(Tstart) == 1)) {
    Tstart <- rep(ifelse(is.null(Tstart), 0, Tstart), n_obs)
  }

  if (length(time) != n_obs || length(status) != n_obs || length(Tstart) != n_obs) {
    stop("Lengths of 'id', 'time', 'status', and 'Tstart' must be identical.")
  }

  if (!all(status %in% c(0, 1, 2))) {
    stop("'status' must only contain values 0 (censored), 1 (event), or 2 (competing risk).")
  }

  if (!is.null(group)) {
    if (length(group) != n_obs) {
      stop("'group' must have the same length as 'id'.")
    }
    if (anyNA(group)) {
      stop("'group' contains missing values (NA). Please remove or impute missing group assignments before proceeding.")
    }
  }


  if (is.null(group)) {
    # No grouping: run on full data, group column = NA
    est <- .mcc_single_group(
      id = id, time = time, status = status, Tstart = Tstart,
      grp_label = NA, ci = ci, conf.level = conf.level, niter = niter, seed = seed
    )
  } else {
    grp_levels <- unique(group)
    grp_list   <- vector("list", length(grp_levels))

    for (k in seq_along(grp_levels)) {
      g   <- grp_levels[k]
      idx <- which(group == g)
      message(sprintf("Processing group: %s", g))
      grp_list[[k]] <- .mcc_single_group(
        id = id[idx], time = time[idx], status = status[idx], Tstart = Tstart[idx],
        grp_label = g, ci = ci, conf.level = conf.level, niter = niter, seed = seed
      )
    }

    est <- dplyr::bind_rows(grp_list)
  }

  # Save original input data for risk tables and summary stats
  raw_data <- data.frame(id = id, time = time, status = status, Tstart = Tstart)
  if (!is.null(group)) raw_data$group <- group

  # Bundle into a list
  res <- list(
    estimates  = est,
    data       = raw_data,
    ci         = ci,
    conf.level = conf.level
  )

  class(res) <- "mcc"
  return(res)
}
