#' Mean Cumulative Count for Left-Truncated Competing Risks Data
#'
#' @description
#' The \code{mcc} package provides a robust computational framework to estimate
#' the cumulative burden of recurrent events over long-term cohort follow-up
#' timelines, accounting safely for competing risks and left-truncation.
#'
#' @keywords internal
#'
#' @importFrom dplyr %>% mutate arrange group_by ungroup filter tally
#' @importFrom dplyr pull left_join rename inner_join slice reframe
#' @importFrom dplyr select summarise across starts_with bind_rows
#' @importFrom rlang .data :=
#' @importFrom stats quantile
#' @importFrom zoo na.locf
#' @importFrom tidyr pivot_wider
#' @importFrom survival survfit Surv
#' @importFrom cmprsk cuminc timepoints
#' @importFrom future.apply future_lapply
"_PACKAGE"
