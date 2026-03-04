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
#' @author Grace Zhou, Department of Biostatistics at St. Jude Children's Hospital \email{grace.zhou@@stjude.org}
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

  if (Event$n == 0) {
    MCC.out <- data.frame(time = as.character(ftime), MCC = rep(0, length(ftime)))
  } else {
    M <- indata2 %>%
      dplyr::filter(.data$status == 1) %>%
      dplyr::group_by(.data$id) %>%
      dplyr::tally() %>%
      dplyr::pull("n") %>%
      max()

    CumI.list <- list()

    for (i in 1:M) {
      in.i.th <- indata2 %>% dplyr::filter(.data$row_num == i)
      rest <- indata2 %>%
        dplyr::filter(.data$n < i & .data$row_num == .data$n) %>%
        dplyr::mutate(status = ifelse(.data$status == 1, 0, .data$status))

      out.i.th <- rbind(in.i.th, rest) %>% dplyr::select(-"row_num", -"n")

      if (calc.trunc) {
        count.data <- mstate::crprep(Tstop = "time", status = "status", data = out.i.th,
                                     trans = 1, cens = 0, Tstart = "tstart", id = "id")

        # crprep does not return weight.trunc if a specific subset has no tstart > 0
        if (!"weight.trunc" %in% names(count.data)) {
          count.data$weight.trunc <- 1
        }

        # Handle NA/NaN/Inf in weights — NA weights arise when survival estimate is
        # unavailable (e.g. time is beyond last observed event in KM estimator).
        # Safe fallback: carry forward last known weight, then fill leading NAs with 1.
        count.data <- count.data %>%
          dplyr::arrange(.data$id, .data$Tstop) %>%
          dplyr::group_by(.data$id) %>%
          dplyr::mutate(
            weight.cens  = ifelse(!is.finite(.data$weight.cens),  NA_real_, .data$weight.cens),
            weight.trunc = ifelse(!is.finite(.data$weight.trunc), NA_real_, .data$weight.trunc),
            weight.cens  = zoo::na.locf(.data$weight.cens,  na.rm = FALSE),
            weight.trunc = zoo::na.locf(.data$weight.trunc, na.rm = FALSE),
            # Fill any remaining leading NAs (no prior value to carry forward) with 1
            weight.cens  = ifelse(is.na(.data$weight.cens),  1, .data$weight.cens),
            weight.trunc = ifelse(is.na(.data$weight.trunc), 1, .data$weight.trunc)
          ) %>%
          dplyr::ungroup()

        # Calculate case weights safely inside the dataframe to avoid length mismatches
        count.data$case_weights <- count.data$weight.cens * count.data$weight.trunc

        fit.i.th <- survival::survfit(survival::Surv(Tstart, Tstop, status == 1) ~ 1,
                                      data = count.data,
                                      weight = count.data$case_weights)

        output <- summary(fit.i.th)
        output2 <- data.frame(time = output$time, CumI = 1 - output$surv)

        CumI.i.th <- data.frame(time = ftime) %>%
          dplyr::left_join(output2, by = 'time') %>%
          dplyr::mutate(CumI = ifelse(dplyr::row_number() == 1 & is.na(.data$CumI), 0, .data$CumI),
                        CumI = zoo::na.locf(.data$CumI),
                        M = i)
      } else {
        fit.i.th <- with(out.i.th, cmprsk::cuminc(time, status))
        CumI.i.th <- cmprsk::timepoints(fit.i.th, ftime)$est[1,]
      }

      CumI.list[[i]] <- CumI.i.th
    }

    MCC.raw <- do.call(rbind, CumI.list)

    if (!calc.trunc) {
      MCC.fill <- apply(MCC.raw, 1, zoo::na.locf, na.rm = FALSE)
      MCC <- rowSums(MCC.fill)
      MCC.out <- data.frame(MCC) %>% dplyr::mutate(time = colnames(MCC.raw))
    } else {
      MCC.fill <- MCC.raw %>% tidyr::pivot_wider(names_from = "M", values_from = "CumI", names_prefix = 'M')
      MCC <- rowSums(MCC.fill[, -1])
      MCC.out <- data.frame(MCC) %>% dplyr::mutate(time = MCC.fill$time) %>% dplyr::select("time", "MCC")
    }
  }

  MCC.out <- MCC.out %>% dplyr::mutate_at('time', as.numeric)
  rownames(MCC.out) <- seq(nrow(MCC.out))

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
#'
#' @return A list containing the confidence intervals and the filled MCC matrix.
#' @author Grace Zhou, Department of Biostatistics at St. Jude Children's Hospital \email{grace.zhou@@stjude.org}
#' @keywords internal
#' @noRd
scumi_ci <- function(id, time, status, Tstart, niter) {

  MCC.out <- scumi(id, time, status, Tstart)

  MCC.event <- MCC.out %>%
    dplyr::arrange(.data$time) %>%
    dplyr::group_by(.data$MCC) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  indata <- data.frame(id = id, time = time, status = status, tstart = Tstart)
  uid <- unique(id)

  for (boot in 1:niter) {
    if (boot %% 100 == 0) message(sprintf("Bootstrap iteration: %d", boot))

    sid <- sample(uid, replace = TRUE)
    sid2 <- data.frame(new.id = seq_along(sid), id = sid)

    bootdata <- dplyr::inner_join(indata, sid2, by = "id", relationship = "many-to-many") %>%
      dplyr::mutate(id = .data$new.id)

    boot.out <- scumi(id = bootdata$id, time = bootdata$time, status = bootdata$status, Tstart = bootdata$tstart)

    rename_col <- paste0("MCC", boot)
    boot.out2 <- boot.out %>% dplyr::rename(!!rename_col := "MCC")

    MCC.event <- MCC.event %>% dplyr::left_join(boot.out2, by = "time")
  }

  MCC.event[1, ][is.na(MCC.event[1, ])] <- 0
  MCC.fill <- MCC.event %>% dplyr::mutate_all(zoo::na.locf)

  quantiles <- function(x) {
    return(stats::quantile(x, probs = c(0.025, 0.975)))
  }

  MCC.result <- t(apply(MCC.fill[, -c(1:2)], 1, quantiles))
  MCC.result2 <- cbind(MCC.fill[, c(1:2)], MCC.result)

  MCC.output <- MCC.out %>%
    dplyr::left_join(MCC.result2, by = c('time', 'MCC')) %>%
    dplyr::group_by(.data$MCC) %>%
    dplyr::reframe(time = .data$time,
                   lci = zoo::na.locf(.data$`2.5%`),
                   uci = zoo::na.locf(.data$`97.5%`)) %>%
    dplyr::ungroup() %>%
    dplyr::select("time", "MCC", "lci", "uci")

  return(list('MCC.output' = MCC.output, 'MCC.fill' = MCC.fill))
}

#' Mean Cumulative Count for Left-Truncated Competing Risks Data
#'
#' @description
#' Estimates the total burden of recurrent events within a population using the Mean
#' Cumulative Count (MCC) method. Specifically designed to handle left-truncated data
#' in the presence of competing risks by pulling accurate, time-dependent censoring
#' and truncation weights.
#'
#' @param id A vector identifying individual subjects.
#' @param time A numeric vector of event or censoring times.
#' @param status A numeric vector indicating event status (e.g., 1 = event of interest, 0 = censored, 2 = competing risk).
#' @param Tstart A numeric vector representing the age or time at study entry (left-truncation). Defaults to 0 if not provided, assuming subjects are followed from time zero.
#' @param ci Logical; if \code{TRUE}, calculates 95\% bootstrap confidence intervals. Default is FALSE
#' @param niter Integer; the number of bootstrap iterations to run if \code{ci = TRUE}. Default is 1000.
#'
#' @return A data frame of class \code{"mcc"} containing the time points, estimated mean cumulative count, and optionally the 95\% bootstrap confidence intervals.
#' @author Grace Zhou, Department of Biostatistics at St. Jude Children's Hospital \email{grace.zhou@@stjude.org}
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
#' @export
#'
#' @examples
#' # Create a sample recurrent event dataset
#' mydata <- data.frame(
#'   id = c(1, 2, 3, 4, 4, 4, 5, 5),
#'   time = c(8, 1, 5, 2, 6, 7, 3, 4),
#'   status = c(0, 0, 2, 1, 1, 1, 1, 2),
#'   tstart = c(0, 0, 0, 0, 2, 6, 0, 3)
#' )
#'
#' # Calculate the Mean Cumulative Count without bootstrap CIs
#' result <- mcc(id = mydata$id, time = mydata$time,
#'               status = mydata$status, Tstart = mydata$tstart)
#'
#' print(result)
mcc <- function(id, time, status, Tstart = NULL, ci = FALSE, niter = 1000) {

  # Internally handle the NULL/0 case
  if (is.null(Tstart)) Tstart <- rep(0, length(time))

  if (ci) {
    res <- scumi_ci(id, time, status, Tstart, niter)[['MCC.output']]
  } else {
    res <- scumi(id, time, status, Tstart)
  }

  # Assign custom class for future summary() and plot() S3 methods
  class(res) <- c("mcc", class(res))

  return(res)
}
