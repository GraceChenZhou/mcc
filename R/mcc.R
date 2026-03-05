#' Sum of Cumulative Incidence
#'
#' Calculates the sum of cumulative incidence for recurrent events in the presence of
#' competing risks and left-truncation. This is an internal helper function.
#'
#' @param id A vector identifying individual subjects.
#' @param time A numeric vector of event or censoring times.
#' @param status A numeric vector indicating event status (1 = event of interest, 0 = censored, 2 = competing risk).
#' @param Tstart A numeric vector representing the age or time at study entry (left-truncation). Defaults to 0 if not provided, assuming subjects are followed from time zero.
#' @return A data frame containing the time points and the estimated mean cumulative count.
#' @author Grace Zhou, Department of Epi and Biostatistics at St. Jude Children's Hospital \email{grace.zhou@stjude.org}
#' @keywords internal
#' @noRd
scumi <- function(id, time, status, Tstart = 0) {

  # Replace magic number with machine precision
  prec <- .Machine$double.eps * 1000

  indata <- data.frame(id = id, time = time, status = status, tstart = Tstart) %>%
    dplyr::mutate(time = ifelse(.data$time == .data$tstart, .data$time + prec, .data$time))

  calc.trunc <- any(indata$tstart != 0)

  indata2 <- indata %>%
    dplyr::arrange(.data$id, .data$time) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::mutate(row_num = dplyr::row_number(), n = dplyr::n()) %>%
    dplyr::ungroup()

  Event <- indata2 %>% dplyr::filter(.data$status == 1) %>% dplyr::tally()
  ftime <- sort(unique(indata2$time))

  # Safer check for zero events
  if (nrow(Event) == 0 || Event$n == 0) {
    MCC.out <- data.frame(time = ftime, MCC = rep(0, length(ftime)))
  } else {
    M <- indata2 %>%
      dplyr::filter(.data$status == 1) %>%
      dplyr::group_by(.data$id) %>%
      dplyr::tally() %>%
      dplyr::pull("n") %>%
      max()

    CumI.list <- vector("list", M)

    for (i in 1:M) {
      in.i.th <- indata2 %>% dplyr::filter(.data$row_num == i)
      rest <- indata2 %>%
        dplyr::filter(.data$n < i & .data$row_num == .data$n) %>%
        dplyr::mutate(status = ifelse(.data$status == 1, 0, .data$status))

      out.i.th <- rbind(in.i.th, rest) %>% dplyr::select(-"row_num", -"n")

      if (calc.trunc) {
        # Note: crprep2 is my customized function (pending Author Dr. Geskus's review and approval)
        count.data <- crprep2(Tstop = "time", status = "status", data = out.i.th,
                              trans = 1, cens = 0, Tstart = "tstart", id = "id")

        if (!"weight.trunc" %in% names(count.data)) {
          count.data$weight.trunc <- 1
        }

        count.data$case_weights <- count.data$weight.cens * count.data$weight.trunc

        fit.i.th <- survival::survfit(survival::Surv(Tstart, Tstop, status == 1) ~ 1,
                                      data = count.data,
                                      weights = count.data$case_weights)

        output <- summary(fit.i.th)
        output2 <- data.frame(time = output$time, CumI = 1 - output$surv)

        CumI.i.th <- data.frame(time = ftime) %>%
          dplyr::left_join(output2, by = 'time') %>%
          dplyr::mutate(CumI = ifelse(dplyr::row_number() == 1 & is.na(.data$CumI), 0, .data$CumI),
                        CumI = zoo::na.locf(.data$CumI, na.rm = FALSE),
                        M = i)
      } else {
        fit.i.th <- cmprsk::cuminc(out.i.th$time, out.i.th$status)
        CumI.i.th <- data.frame(time = ftime,
                                CumI = cmprsk::timepoints(fit.i.th, ftime)$est[1,],
                                M = i)
      }

      CumI.list[[i]] <- CumI.i.th
    }

    MCC.raw <- do.call(rbind, CumI.list)

    if (!calc.trunc) {
      MCC.fill <- MCC.raw %>%
        dplyr::group_by(.data$M) %>%
        dplyr::mutate(CumI = zoo::na.locf(.data$CumI, na.rm = FALSE)) %>%
        dplyr::ungroup()

      MCC.out <- MCC.fill %>%
        dplyr::group_by(.data$time) %>%
        dplyr::summarise(MCC = sum(.data$CumI, na.rm = TRUE)) %>%
        dplyr::ungroup()
    } else {
      MCC.fill <- MCC.raw %>% tidyr::pivot_wider(names_from = "M", values_from = "CumI", names_prefix = 'M')
      # rowSums needs to ignore the 'time' column safely
      MCC <- rowSums(dplyr::select(MCC.fill, -time), na.rm = TRUE)
      MCC.out <- data.frame(time = MCC.fill$time, MCC = MCC)
    }
  }

  MCC.out <- MCC.out %>% dplyr::mutate(time = as.numeric(.data$time))
  return(MCC.out)
}

#' Bootstrap 95% Confidence Interval for Sum of Cumulative Incidence
#'
#' Calculates the bootstrap 95% confidence intervals for the MCC estimates.
#' This is an internal helper function.
#'
#' @param id A vector identifying individual subjects.
#' @param time A numeric vector of event or censoring times.
#' @param status A numeric vector indicating event status (1 = event of interest, 0 = censored, 2 = competing risk).
#' @param Tstart A numeric vector representing the age or time at study entry (left-truncation). Defaults to 0 if not provided, assuming subjects are followed from time zero.
#' @param niter Integer specifying the number of bootstrap iterations. Default is 1000.
#' @param seed Optional integer to set the seed for reproducible bootstrap results.
#'
#' @return A list containing the confidence intervals and the filled MCC matrix.
#' @author Grace Zhou, Department of Epi and Biostatistics at St. Jude Children's Hospital \email{grace.zhou@stjude.org}
#' @keywords internal
#' @noRd
scumi_ci <- function(id, time, status, Tstart, niter, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  MCC.out <- scumi(id, time, status, Tstart)

  MCC.event <- MCC.out %>%
    dplyr::arrange(.data$time) %>%
    dplyr::group_by(.data$MCC) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  indata <- data.frame(id = id, time = time, status = status, tstart = Tstart)
  uid <- unique(id)

  # Pre-allocate list for performance instead of joining in a loop
  boot_list <- vector("list", niter)

  for (boot in 1:niter) {
    if (boot %% 100 == 0) message(sprintf("Bootstrap iteration: %d", boot))

    sid <- sample(uid, replace = TRUE)
    sid2 <- data.frame(new.id = seq_along(sid), id = sid)

    bootdata <- dplyr::inner_join(indata, sid2, by = "id", relationship = "many-to-many") %>%
      dplyr::mutate(id = .data$new.id)

    boot.out <- scumi(id = bootdata$id, time = bootdata$time, status = bootdata$status, Tstart = bootdata$tstart)

    # Store isolated results with a boot_id tag
    boot_list[[boot]] <- boot.out %>%
      dplyr::select("time", "MCC") %>%
      dplyr::mutate(boot_id = paste0("MCC", boot))
  }

  # Bind all iterations and pivot wider simultaneously (O(N) operation)
  all_boots <- dplyr::bind_rows(boot_list)
  boot_matrix <- all_boots %>%
    tidyr::pivot_wider(names_from = "boot_id", values_from = "MCC")

  # Join the consolidated matrix back to the event times
  MCC.event <- MCC.event %>% dplyr::left_join(boot_matrix, by = "time")

  # Handle baseline NAs safely
  MCC.event[1, ][is.na(MCC.event[1, ])] <- 0

  MCC.fill <- MCC.event %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("MCC"), ~zoo::na.locf(.x, na.rm = FALSE)))

  quantiles <- function(x) {
    return(stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
  }

  MCC.boot_cols <- grep("^MCC[0-9]+", names(MCC.fill))
  MCC.result <- t(apply(MCC.fill[, MCC.boot_cols], 1, quantiles))
  colnames(MCC.result) <- c("lci", "uci")

  MCC.result2 <- cbind(MCC.fill[, c("time", "MCC")], MCC.result)

  MCC.output <- MCC.out %>%
    dplyr::left_join(MCC.result2, by = c('time', 'MCC')) %>%
    dplyr::group_by(.data$MCC) %>%
    dplyr::reframe(time = .data$time,
                   lci = zoo::na.locf(.data$lci, na.rm = FALSE),
                   uci = zoo::na.locf(.data$uci, na.rm = FALSE)) %>%
    dplyr::ungroup() %>%
    dplyr::select("time", "MCC", "lci", "uci")

  return(list('MCC.output' = MCC.output, 'MCC.fill' = MCC.fill))
}

#' Mean Cumulative Count for Left-Truncated Competing Risks Data
#'
#' @description
#' Estimates the total burden of recurrent events within a population using the Mean
#' Cumulative Count (MCC) method. Specifically designed to handle left-truncated data
#' in the presence of competing risks by applying time-dependent censoring
#' and truncation weights.
#'
#' @param id A vector identifying individual subjects.
#' @param time A numeric vector of event or censoring times.
#' @param status A numeric vector indicating event status (e.g., 1 = event of interest, 0 = censored, 2 = competing risk).
#' @param Tstart A numeric vector representing the age or time at study entry (left-truncation). Defaults to 0 if not provided, assuming subjects are followed from time zero.
#' @param ci Logical; if \code{TRUE}, calculates 95\% bootstrap confidence intervals. Default is FALSE.
#' @param niter Integer; the number of bootstrap iterations to run if \code{ci = TRUE}. Default is 1000.
#' @param seed Optional integer to set the seed for reproducible bootstrap results.
#'
#' @return A data frame of class \code{"mcc"} containing the time points, estimated mean cumulative count, and optionally the 95\% bootstrap confidence intervals.
#' @author Grace Zhou, Department of Biostatistics at St. Jude Children's Hospital \email{grace.zhou@stjude.org}
#' @importFrom dplyr %>% mutate arrange group_by ungroup filter tally pull left_join rename inner_join slice reframe select summarise across starts_with bind_rows
#' @importFrom rlang .data :=
#' @importFrom stats quantile
#' @importFrom zoo na.locf
#' @importFrom tidyr pivot_wider
#' @importFrom survival survfit Surv
#' @importFrom cmprsk cuminc timepoints
#' @export
#'
#' @examples
#' # Create a sample recurrent event dataset
#' mydata <- data.frame(
#'   id = c(1, 2, 3, 4, 4, 4, 5, 5),
#'   tstop = c(8, 1, 5, 2, 6, 7, 3, 4),
#'   status = c(0, 0, 2, 1, 1, 1, 1, 2),
#'   tstart = c(0, 0, 0, 0, 2, 6, 0, 3)
#' )
#'
#' # Calculate the Mean Cumulative Count without bootstrap CIs
#' result <- mcc(id = mydata$id, time = mydata$tstop,
#'                status = mydata$status, Tstart = mydata$tstart)
#'
#' print(result)
#'
#' # Calculate the Mean Cumulative Count with bootstrap CIs
#' result_ci <- mcc(id = mydata$id, time = mydata$tstop,
#'                status = mydata$status, Tstart = mydata$tstart,
#'                ci = TRUE, niter = 100, seed = 123)
#'
#' print(result_ci)
mcc <- function(id, time, status, Tstart = NULL, ci = FALSE, niter = 1000, seed = NULL) {

  if (is.null(Tstart)) Tstart <- rep(0, length(time))

  if (ci) {
    res_list <- scumi_ci(id = id, time = time, status = status, Tstart = Tstart, niter = niter, seed = seed)
    res <- res_list[['MCC.output']]
  } else {
    res <- scumi(id = id, time = time, status = status, Tstart = Tstart)
  }

  if (is.null(res)) {
    stop("Calculation failed: Internal helper functions returned NULL.")
  }

  class(res) <- c("mcc", "data.frame")
  return(res)
}
