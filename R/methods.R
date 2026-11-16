#' Print Method for MCC Objects
#'
#' @param x An object of class \code{mcc}.
#' @param ... Additional arguments passed to \code{print}.
#' @importFrom utils head
#' @export
print.mcc <- function(x, ...) {
  cat("Mean Cumulative Count (MCC) Estimates\n")
  cat("=======================================\n")

  # Detect grouping assignments properly from the raw data
  has_group <- "group" %in% names(x$data) && !all(is.na(x$data$group))

  # Internal helper to print the cohort breakdown for a given sub-dataset
  print_cohort_breakdown <- function(df, label = NULL) {
    if (!is.null(label)) {
      cat(sprintf("\n--- Group: %s ---\n", label))
    }

    # Summarize history per individual subject in this subset
    subj_summary <- df %>%
      dplyr::group_by(.data$id) %>%
      dplyr::arrange(.data$time) %>%
      dplyr::summarise(
        events = sum(.data$status == 1, na.rm = TRUE),
        terminal_status = dplyr::last(.data$status),
        .groups = "drop"
      )

    # Sample size and cohort breakdown logic
    n_total         <- nrow(subj_summary)
    n_had_events  <- sum(subj_summary$events > 0)

    # Isolate specific sub-cohort states
    n_cens_no_event <- sum(subj_summary$events == 0 & subj_summary$terminal_status == 0)
    n_cr_no_event   <- sum(subj_summary$events == 0 & subj_summary$terminal_status == 2)

    # Terminal status summaries (regardless of prior event history)
    term_censored  <- sum(subj_summary$terminal_status == 0)
    term_event     <- sum(subj_summary$terminal_status == 1)
    term_comp_risk <- sum(subj_summary$terminal_status == 2)

    # Console Output
    cat(sprintf("Total Sample Size (N):                  %d\n", n_total))
    cat(sprintf("  - Subjects with >= 1 recurrent event: %d\n", n_had_events))
    cat(sprintf("  - Censored without experiencing event: %d\n", n_cens_no_event))
    cat(sprintf("  - Competing risk without event:       %d\n", n_cr_no_event))
    cat("\nTerminal Status Summary (Last Observation)\n")
    cat(sprintf("  - Censored (Dropped out/End of study): %d\n", term_censored))
    cat(sprintf("  - Experiencing Event of Interest:     %d\n", term_event))
    cat(sprintf("  - Experiencing Competing Risk:         %d\n", term_comp_risk))
    cat("\nRecurrent Event Distribution Burden:\n")
    cat(sprintf("  - 0 Events:  %d subjects\n", sum(subj_summary$events == 0)))
    cat(sprintf("  - 1 Event:   %d subjects\n", sum(subj_summary$events == 1)))
    cat(sprintf("  - 2 Events:  %d subjects\n", sum(subj_summary$events == 2)))
    cat(sprintf("  - 3+ Events: %d subjects\n", sum(subj_summary$events >= 3)))
  }

  # Execute printing logic based on stratification existence
  if (has_group) {
    cat("Stratified Analysis\n")
    grps <- sort(unique(x$data$group))
    for (g in grps) {
      group_data <- x$data[x$data$group == g, ]
      print_cohort_breakdown(group_data, label = g)
    }
    cat("\n=======================================\n")
  } else {
    cat("Overall Cohort (Unstratified Analysis)\n")
    print_cohort_breakdown(x$data)
    cat("\n=======================================\n")
  }

  if (x$ci) {
    cat(sprintf("Confidence Intervals:  %g%% Bootstrap\n", x$conf.level * 100))
  }

  cat("\nPreview of estimates:\n")

  # Dynamically handles previewing rows conditionally based on grouping rules
  estimates_df <- as.data.frame(x$estimates)
  if (has_group) {
    preview_df <- estimates_df %>%
      dplyr::group_by(.data$group) %>%
      dplyr::slice_head(n = 3) %>%
      dplyr::ungroup()
    cat("Showing the first 3 estimates for each group:\n")
  } else {
    preview_df <- utils::head(estimates_df, n = 6)
  }

  print(as.data.frame(preview_df))

  invisible(x)
}

#' Summary Method for MCC Objects
#'
#' Evaluates mean cumulative count estimates at specified time points.
#' If the analysis is stratified, results are automatically reshaped into a wide
#' format for cross-group comparison.
#'
#' @param object An object of class \code{mcc}.
#' @param times An optional numeric vector of specific time points to evaluate (e.g., c(1, 5, 10)).
#'   If omitted, defaults to 5 evenly spaced across the follow-up window.
#' @param ... Additional arguments passed to \code{summary}.
#'
#' @return A data frame containing mcc estimate at the corresponding time points. Reshaped into a wide format if grouping is present.
#' @export

summary.mcc <- function(object, times = NULL, ...) {
  est <- object$estimates

  if (is.null(times)) {
    max_t <- max(est$time, na.rm = TRUE)
    min_t <- min(est$time, na.rm = TRUE)
    times <- seq(min_t, max_t, length.out = 5)
  }

  # Step-function rollback helper
  get_milestones <- function(sub_df, t_points) {
    sub_df <- sub_df[order(sub_df$time), ]
    res_list <- lapply(t_points, function(t) {
      row_idx <- which(sub_df$time <= t)
      if (length(row_idx) == 0) {
        out <- sub_df[1, ]
        out$time <- t
        out$MCC <- 0
        if ("lci" %in% names(sub_df)) { out$lci <- 0; out$uci <- 0 }
        return(out)
      }
      out <- sub_df[max(row_idx), ]
      out$time <- t
      return(out)
    })
    do.call(rbind, res_list)
  }

  has_group <- "group" %in% names(est) && !all(is.na(est$group))

  if (has_group) {
    grps <- unique(est$group)
    long_milestones <- do.call(rbind, lapply(grps, function(g) {
      get_milestones(est[est$group == g, ], times)
    }))

    # Dynamically pivot to wide format based on whether CIs are included
    if (object$ci) {
      # Construct the custom label based on the object's specific conf.level metadata
      ci_pct <- object$conf.level * 100
      suffix_label <- sprintf(" (%g%% CI)", ci_pct)

      # Combine estimates into "MCC (LCI, UCI)" string format
      milestone_df <- long_milestones %>%
        dplyr::mutate(
          est_str = sprintf("%.3f (%.3f, %.3f)", .data$MCC, .data$lci, .data$uci)
        ) %>%
        dplyr::select("time", "group", "est_str") %>%
        tidyr::pivot_wider(
          names_from = "group",
          values_from = "est_str",
          names_prefix = "MCC_",
          names_glue = "{.value}_{group} " # Dummy spacer used to clean up the suffix paste downstream
        )

      # Apply your custom dynamic header style cleanly across the pivoted columns
      colnames(milestone_df) <- gsub("est_str_", "MCC_", colnames(milestone_df))
      colnames(milestone_df) <- ifelse(
        grepl("^MCC_", colnames(milestone_df)),
        paste0(trimws(colnames(milestone_df)), suffix_label),
        colnames(milestone_df)
      )
    } else {
      # If no CIs, pivot the raw numeric MCC values directly
      milestone_df <- long_milestones %>%
        dplyr::select("time", "group", "MCC") %>%
        tidyr::pivot_wider(
          names_from = "group",
          values_from = "MCC",
          names_prefix = "MCC_"
        )
    }
  } else {
    milestone_df <- get_milestones(est, times)
    if (object$ci) {
      milestone_df <- milestone_df %>%
        dplyr::select("time", "MCC", "lci", "uci")
    } else {
      milestone_df <- milestone_df %>%
        dplyr::select("time", "MCC")
    }
  }

  # Strip internal tibble properties to keep standard df classing
  milestone_df <- as.data.frame(milestone_df)
  class(milestone_df) <- c("summary.mcc", "data.frame")
  return(milestone_df)
}

#' Print Method for summary.mcc Objects
#'
#' @param x An object of class \code{summary.mcc}.
#' @param ... Additional arguments passed to \code{print}.
#' @export
print.summary.mcc <- function(x, ...) {
  cat("Mean Cumulative Count (MCC) Summary\n")
  cat("=============================================\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' Convert MCC Object to Data Frame
#'
#' Extracts the structured data frame of estimates from a calculated \code{mcc} object.
#'
#' @param x An object of class \code{mcc}.
#' @param ... Additional arguments passed to \code{as.data.frame}.
#'
#' @return A data frame containing the time points and mean cumulative count estimates.
#' @export
as.data.frame.mcc <- function(x, ...) {
  return(x$estimates)
}
