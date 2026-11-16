#' Plot Mean Cumulative Count Curves with Number-at-Risk Table
#'
#' @description Plots MCC step curves, optional confidence bands, and subject-level
#'   number-at-risk counts. Left-truncation uses entry < t <= exit.
#' @param mcc_result An `mcc` object or an estimates data frame.
#' @param data Original subject/event data. Optional for an `mcc` object.
#' @param id,time,Tstart Column names for subject ID, exit time, and entry time.
#' @param mcc_group,data_group Group column names in estimates and original data.
#' @param group_value Groups to display; `NULL` displays all groups.
#' @param control_group Optional control group (dashed by default).
#' @param risk_times Numeric times for the number-at-risk table.
#' @param color_values,line_values Optional named vectors indexed by original group names.
#' @param conf_int Show confidence bands if available.
#' @param conf_alpha Confidence band opacity.
#' @param xlab,ylab,title Axis labels and title.
#' @param title_hjust Title alignment.
#' @param legend Whether to display a legend.
#' @param legend.position,legend.title Legend placement and title.
#' @param palette Optional RColorBrewer palette name.
#' @param line_size Step line width.
#' @param ggtheme,tables.theme Themes for curve and risk table.
#' @param risk_table_title,show_risk_table_title Risk table label and visibility.
#' @param risk_table_title_size,risk_table_text_size,risk_table_y_text_size Risk table font sizes.
#' @param title_size,axis_label_size,axis_text_size Plot font sizes.
#' @param legend_text_size,legend_title_size Legend font sizes.
#' @param plot_heights Relative heights of curve and risk table.
#' @return A list with `plot`, `p_mcc`, `p_risk`, `plot_data`, `risk_table`, and `risk_times`.
#' @export


plot_mcc_with_risk_table <- function(mcc_result,
                                     data,
                                     id = "id", #Subject ID column in the original data
                                     time = "time", #Event/censoring/competing risk time column in the original data
                                     Tstart = "Tstart", #Optional left-truncation entry time column; set NULL if not available
                                     mcc_group = "group", #Group variable name in mcc_result; Default: "group"; Used for MCC curve stratification.
                                     data_group = "group", #Group variable name in original dataset; Default: "group"; Used for risk table calculation
                                     group_value = NULL, #Character vector specifying groups to display (Ex. c("A", "B")); if NULL, all groups are plotted
                                     control_group = NULL, #Optional control group; shown as dashed line by default
                                     risk_times = NULL, #User-defined time points for the risk table
                                     color_values = NULL, #Optional named vector for group colors, using original group names, e.g., c(A = "red", B = "blue", C = "black")
                                     line_values = NULL, #Optional named vector for line types, using original group names, e.g., c(A = "solid", B = "solid", C = "dashed")
                                     conf_int = TRUE, #Whether to show confidence interval ribbon
                                     conf_alpha = 0.15, #Transparency of confidence interval ribbon
                                     xlab = "Time", #X-axis label
                                     ylab = "Cumulative burden per individual", #Y-axis label
                                     title = "MCC curves by group", #Plot title
                                     title_hjust = 0.5, #Title alignment: 0 = left, 0.5 = center, 1 = right
                                     legend = TRUE, #Whether to show legend
                                     legend.position = c(0.10, 0.85), #Legend position; can be "right" or coordinates
                                     legend.title = NULL, #Optional legend title
                                     palette = NULL, #Optional RColorBrewer palette, e.g., "Set1", "Dark2"
                                     line_size = 0.5, #MCC curve line width
                                     ggtheme = ggplot2::theme_classic(), #Theme for MCC curve
                                     tables.theme = ggtheme, #Theme for risk table
                                     risk_table_title = "Number at risk", #risk table y-axis title
                                     show_risk_table_title = TRUE, #Whether to show risk table title
                                     risk_table_title_size = 10, #Font size of risk table title
                                     risk_table_text_size = 3, #Font size of risk table numbers
                                     risk_table_y_text_size = 8, #Font size of risk table group labels
                                     title_size = 14, #Font size of plot title
                                     axis_label_size = 12, #Font size of x/y axis labels
                                     axis_text_size = 10, #Font size of x/y axis tick labels
                                     legend_text_size = 10, #Font size of legend labels
                                     legend_title_size = 10, #Font size of legend title
                                     plot_heights = NULL) { #Relative heights of MCC plot and risk table

  if (inherits(mcc_result, "mcc")) {
    if (missing(data) || is.null(data)) data <- mcc_result$data
    mcc_result <- mcc_result$estimates
  }
  if (missing(data) || is.null(data)) stop("'data' must be supplied for data-frame estimates.")
  mcc_result <- as.data.frame(mcc_result)
  data <- as.data.frame(data)
  user_specified_group <- !is.null(group_value)
  if (!mcc_group %in% names(mcc_result)) mcc_result[[mcc_group]] <- NA_character_
  if (!data_group %in% names(data)) data[[data_group]] <- NA_character_
  if (conf_int && !all(c("lci", "uci") %in% names(mcc_result))) {
    message("Confidence intervals are unavailable; plotting MCC estimates without ribbons.")
    conf_int <- FALSE
  }
  # Check required columns in MCC result
  required_mcc_cols <- c("time", "MCC", mcc_group)
  if (!all(required_mcc_cols %in% names(mcc_result))) {
    stop("mcc_result must contain columns: ", paste(required_mcc_cols, collapse = ", "))
  }

  # Check whether left-truncation time is available
  has_Tstart <- !is.null(Tstart) && Tstart %in% names(data)

  # Check required columns in original data
  required_data_cols <- c(id, time, data_group)
  if (has_Tstart) {
    required_data_cols <- c(required_data_cols, Tstart)
  }

  if (!all(required_data_cols %in% names(data))) {
    stop("data must contain columns: ", paste(required_data_cols, collapse = ", "))
  }

  # Preserve overall analyses before filtering requested groups.
  overall <- all(is.na(mcc_result[[mcc_group]])) && all(is.na(data[[data_group]]))
  if (overall) {
    if (user_specified_group) warning("group_value ignored for an unstratified MCC result.")
    group_value <- NULL
  }
  # Check whether requested groups exist
  available_groups <- unique(as.character(data[[data_group]]))

  # If control group is provided, include it in selected groups
  if (!is.null(group_value) && !is.null(control_group)) {
    group_value <- unique(c(group_value, control_group))
  }

  # Check whether requested groups exist
  if (!is.null(group_value)) {
    missing_groups <- group_value[!group_value %in% available_groups]
    # Stop if none of the requested groups exist
    if (all(!group_value %in% available_groups)) {
      stop("None of the requested groups were found in the data.")
    }
    # Warn if only some requested groups do not exist
    if (length(missing_groups) > 0) {
      warning("Requested group(s) not found: ", paste(missing_groups, collapse = ", "))
    }
  }

  # Filter selected groups if group_value is provided
  if (!is.null(group_value)) {
    mcc_result <- mcc_result %>%
      dplyr::filter(.data[[mcc_group]] %in% group_value)

    data <- data %>%
      dplyr::filter(.data[[data_group]] %in% group_value)
  }

  if (nrow(mcc_result) == 0 || nrow(data) == 0) {
    stop("None of the requested groups were found in the data.")
  }

  # Handle missing group information
  # If all group values are NA, plot one overall curve and suppress legend
  missing_group_info <- all(is.na(mcc_result[[mcc_group]])) || all(is.na(data[[data_group]]))

  # Drop records with partially missing group values
  if (!missing_group_info && (anyNA(mcc_result[[mcc_group]]) || anyNA(data[[data_group]]))) {
    message("NA values were found in the group variable. ",
            "Records with missing group values will be dropped before plotting.")
    mcc_result <- mcc_result[!is.na(mcc_result[[mcc_group]]), ]
    data <- data[!is.na(data[[data_group]]), ]
  }

  if (missing_group_info) {
    if (user_specified_group) {
      warning("Group values are missing, so the specified group_value cannot be applied. ")
    } else {
      message("All group values are missing. ","Plotting one unstratified MCC curve.")
    }
    mcc_result[[mcc_group]] <- "Unspecified"
    data[[data_group]] <- "Unspecified"
    control_group <- NULL
  }

  # Show legend only when requested and valid group information exists
  show_legend <- legend && !missing_group_info

  if (nrow(mcc_result) == 0 || nrow(data) == 0) stop("No observations remain for plotting.")
  # Define observed time range from MCC results
  mcc_time_range <- range(mcc_result$time, na.rm = TRUE)

  # Track whether risk_times were specified by the user
  user_specified_risk_times <- !is.null(risk_times)

  # Validate user-specified risk times
  if (user_specified_risk_times) {

    if (!is.numeric(risk_times) ||
        length(risk_times) == 0L ||
        anyNA(risk_times) ||
        any(!is.finite(risk_times))) {

      stop(
        "'risk_times' must be a non-empty numeric vector ",
        "containing only finite, non-missing values.",
        call. = FALSE
      )
    }
  }

  # Default risk table time points
  if (!user_specified_risk_times) {
    # Use pretty() to obtain convenient interior time points
    interior_times <- pretty(mcc_time_range, n = 5)
    # Keep only points strictly inside the observed time range
    interior_times <- interior_times[interior_times > mcc_time_range[1] & interior_times < mcc_time_range[2]]
    # Round default interior points to the nearest integer
    interior_times <- round(interior_times)
    # Keep the exact observed minimum and maximum as endpoints
    risk_times <- c(mcc_time_range[1], interior_times, mcc_time_range[2])
    # Remove possible duplicates after rounding
    risk_times <- sort(unique(risk_times))
    # Default x-axis range
    x_lim <- mcc_time_range
    } else {
    # Preserve exactly the user-specified risk times
    risk_times <- sort(unique(risk_times))
    if (min(risk_times, na.rm = TRUE) < mcc_time_range[1] ||
        max(risk_times, na.rm = TRUE) > mcc_time_range[2]) {
      warning("User-specified risk_times fall outside the observed time range in mcc_result.")
    }
    # Use user-specified range even if outside MCC range
    x_lim <- range(risk_times, na.rm = TRUE)
    }

  # Labels only affect what is printed on the x-axis
  if (user_specified_risk_times) {
    x_labels <- risk_times
  } else {
    x_labels <- round(risk_times)
  }

  if (!is.numeric(risk_times) || !length(risk_times) || any(!is.finite(risk_times))) {
    stop("'risk_times' must contain finite numeric values.")
  }
  # Create long-format risk table
  risk_table <- expand.grid(
    group = unique(as.character(data[[data_group]])),
    time = risk_times,
    stringsAsFactors = FALSE
  )

  # Calculate number at risk
  # If Tstart exists: use left-truncation risk set, Tstart < t <= time
  # If Tstart does not exist: use standard risk set, time >= t
  risk_table$n_risk <- mapply(function(g, t) {
    if (has_Tstart) {
      at_risk <- data[[data_group]] == g & data[[Tstart]] < t & data[[time]] >= t
    } else {
      at_risk <- data[[data_group]] == g & data[[time]] >= t
    }
    length(unique(data[[id]][!is.na(at_risk) & at_risk & !is.na(data[[id]])]))
  }, risk_table$group, risk_table$time)

  # Standardize group variables for plotting
  plot_data <- mcc_result %>%
    dplyr::mutate(plot_group = as.character(.data[[mcc_group]]))

  risk_table <- risk_table %>%
    dplyr::mutate(plot_group = as.character(.data$group))

  # Add control label if control_group is provided
  if (!is.null(control_group)) {
    plot_data <- plot_data %>%
      dplyr::mutate(legend_group = ifelse(
        .data$plot_group == control_group,
        paste0(.data$plot_group, " (Control)"),
        .data$plot_group))

    risk_table <- risk_table %>%
      dplyr::mutate(legend_group = ifelse(
        .data$plot_group == control_group,
        paste0(.data$plot_group, " (Control)"),
        .data$plot_group))
  } else {
    plot_data <- plot_data %>% dplyr::mutate(legend_group = .data$plot_group)
    risk_table <- risk_table %>% dplyr::mutate(legend_group = .data$plot_group)
  }

  # Limit risk-table rows to groups with an estimated MCC curve.
  risk_table <- risk_table[risk_table$plot_group %in% unique(plot_data$plot_group), , drop = FALSE]
  # Keep legend order consistent
  legend_levels <- unique(plot_data$legend_group)
  plot_data$legend_group <- factor(plot_data$legend_group, levels = legend_levels)
  risk_table$legend_group <- factor(risk_table$legend_group, levels = legend_levels)

  # Map original group names to displayed legend labels
  group_label_map <- stats::setNames(
    as.character(unique(plot_data$legend_group)),
    as.character(unique(plot_data$plot_group)))

  # Set default plot heights
  # Increase risk table height when more groups are displayed
  if (is.null(plot_heights)) {
    n_groups <- length(legend_levels)
    plot_heights <- c(4, max(0.8, 0.3 * n_groups))
  }

  # Set line types
  # Default: all groups use solid lines, with the control group dashed.
  # User-specified line types override the defaults for selected groups.

  # Create default line types for all plotted groups
  line_values_final <- stats::setNames(
    rep("solid", length(legend_levels)), legend_levels)

  # Set the control group to dashed by default
  if (!is.null(control_group)) {
    control_label <- paste0(control_group, " (Control)")
    line_values_final[control_label] <- "dashed"
  }

  # Override defaults with user-specified line types
  if (!is.null(line_values)) {
    # Check whether user-specified groups are included in the plotted data
    missing_line_groups <- setdiff(names(line_values), names(group_label_map))

    if (length(missing_line_groups) > 0) {
      warning("Line type value(s) specified for group(s) not found in the plotted data: ",
              paste(missing_line_groups, collapse = ", "))
    }

    # Keep only valid groups
    valid_line_values <- line_values[names(line_values) %in% names(group_label_map)]

    # Convert original group names to displayed group labels
    names(valid_line_values) <- group_label_map[names(valid_line_values)]

    # Override default line types
    line_values_final[names(valid_line_values)] <- valid_line_values
  }

  # Use the final complete line type vector for plotting
  line_values <- line_values_final

  # Set custom colors if provided
  if (!is.null(color_values)) {
    color_values <- color_values[names(color_values) %in% names(group_label_map)]
    names(color_values) <- group_label_map[names(color_values)]
  }

  if (!is.null(color_values)) {
    all_labels <- as.character(unique(plot_data$legend_group))
    absent <- setdiff(all_labels, names(color_values))
    if (length(absent)) stop("'color_values' must specify every plotted group: ", paste(absent, collapse = ", "))
  }
  # MCC curve
  p_mcc <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$time,
      y = .data$MCC,
      color = .data$legend_group,
      fill = .data$legend_group,
      linetype = .data$legend_group)
    )

  # Add confidence interval ribbon if requested
  if (conf_int) {
    p_mcc <- p_mcc +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), alpha = conf_alpha, color = NA)
  }

  p_mcc <- p_mcc +
    ggplot2::geom_step(linewidth = line_size) +
    ggplot2::scale_linetype_manual(values = line_values) +
    ggplot2::scale_x_continuous(breaks = risk_times, labels = x_labels, expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = x_lim) +
    ggplot2::labs(x = xlab, y = ylab, color = legend.title, fill = legend.title, linetype = legend.title, title = title) +
    ggtheme +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, hjust = title_hjust),
      axis.title = ggplot2::element_text(size = axis_label_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      legend.position = if (show_legend) legend.position else "none",
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.title = ggplot2::element_text(size = legend_title_size)
    )

  # Apply custom colors or palette
  if (!is.null(color_values)) {
    p_mcc <- p_mcc +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::scale_fill_manual(values = color_values)
  } else if (!is.null(palette)) {
    p_mcc <- p_mcc +
      ggplot2::scale_color_brewer(palette = palette) +
      ggplot2::scale_fill_brewer(palette = palette)
  }

  # Risk table plot used under MCC curve
  p_risk <- ggplot2::ggplot(
    risk_table,
    ggplot2::aes(x = .data$time, y = .data$group, label = .data$n_risk, color = .data$legend_group)) +
    ggplot2::geom_text(size = risk_table_text_size) +
    ggplot2::scale_x_continuous(breaks = risk_times, expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = x_lim, clip = "off") +
    ggplot2::labs(x = NULL, y = if (show_risk_table_title) risk_table_title else NULL) +
    tables.theme +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = risk_table_title_size),
      axis.text.y = ggplot2::element_text(size = risk_table_y_text_size, margin = ggplot2::margin(r = 10)),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 40))

  # Apply custom colors or palette to risk table
  if (!is.null(color_values)) {
    p_risk <- p_risk + ggplot2::scale_color_manual(values = color_values)
  } else if (!is.null(palette)) {
    p_risk <- p_risk + ggplot2::scale_color_brewer(palette = palette)
  }

  # Combine MCC curve and risk table
  final_plot <- patchwork::wrap_plots(p_mcc, p_risk, ncol = 1, heights = plot_heights)

  # Return plot and underlying data for customization
  return(list(
    plot = final_plot,
    p_mcc = p_mcc,
    p_risk = p_risk,
    plot_data = plot_data,
    risk_table = risk_table,
    risk_times = risk_times
  ))
}




#' Plot an MCC Object
#'
#' @param x An object of class `mcc`.
#' @param ... Additional arguments to [plot_mcc_with_risk_table()].
#' @return The combined patchwork plot, invisibly.
#' @export
plot.mcc <- function(x, ...) {
  out <- plot_mcc_with_risk_table(x, ...)
  print(out$plot)
  invisible(out$plot)
}
