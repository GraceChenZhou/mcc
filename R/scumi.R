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
#' @importFrom dplyr %>%
#' @author Grace Zhou, Department of Biostatistics at St. Jude Children's Research Hospital \email{grace.zhou@stjude.org}
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

  ftime <- sort(unique(indata2$time))

  # Single pass: count events per subject (used for both zero-event check and M)
  event_counts <- indata2 %>%
    dplyr::filter(.data$status == 1) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::tally()

  if (nrow(event_counts) == 0) {
    MCC.out <- data.frame(time = ftime, MCC = rep(0, length(ftime)))
  } else {
    M <- max(event_counts$n)

    CumI.list <- vector("list", M)

    for (i in 1:M) {
      in.i.th <- indata2 %>% dplyr::filter(.data$row_num == i)
      rest <- indata2 %>%
        dplyr::filter(.data$n < i & .data$row_num == .data$n) %>%
        dplyr::mutate(status = ifelse(.data$status == 1, 0, .data$status))

      out.i.th <- rbind(in.i.th, rest) %>% dplyr::select(-"row_num", -"n")

      if (calc.trunc) {
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
#' @param conf.level A numeric value between 0 and 1 specifying the confidence level. Default is 0.95 (95\% CI).
#' @param seed Optional integer to set the seed for reproducible bootstrap results.
#'
#' @return A list containing the confidence intervals and the filled MCC matrix.
#' @author Grace Zhou, Department of Biostatistics at St. Jude Children's Research Hospital \email{grace.zhou@stjude.org}
#' @keywords internal
#' @importFrom future.apply future_lapply
#' @noRd
scumi_ci <- function(id, time, status, Tstart, niter, conf.level = 0.95, seed = NULL) {

  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single numeric value strictly between 0 and 1 (e.g. 0.95).")
  }

  alpha <- 1 - conf.level
  probs <- c(alpha / 2, 1 - alpha / 2)

  MCC.out <- scumi(id, time, status, Tstart)

  MCC.event <- MCC.out %>%
    dplyr::arrange(.data$time) %>%
    dplyr::group_by(.data$MCC) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  indata <- data.frame(id = id, time = time, status = status, tstart = Tstart)
  uid <- unique(id)

  # Execute bootstrap in parallel using L'Ecuyer-CMRG RNG streams safely
  boot_list <- future.apply::future_lapply(1:niter, function(boot) {
    sid <- sample(uid, replace = TRUE)
    sid2 <- data.frame(new.id = seq_along(sid), id = sid)

    bootdata <- dplyr::inner_join(indata, sid2, by = "id", relationship = "many-to-many") %>%
      dplyr::mutate(id = .data$new.id)

    boot.out <- scumi(id = bootdata$id, time = bootdata$time, status = bootdata$status, Tstart = bootdata$tstart)

    boot.out %>%
      dplyr::select("time", "MCC") %>%
      dplyr::mutate(boot_id = paste0("MCC", boot))

  }, future.seed = ifelse(is.null(seed), TRUE, seed))


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

  MCC.result <- MCC.fill %>%
    dplyr::select(dplyr::starts_with("MCC")) %>%
    apply(1, function(x) stats::quantile(x, probs = probs, na.rm = TRUE)) %>%
    t()
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

#' Run MCC for a single group subset
#'
#' @param id,time,status,Tstart Subset vectors for one group.
#' @param grp_label Label to attach as the \code{group} column.
#' @param ci,niter,seed Passed through from \code{mcc()}.
#' @return A data frame with a leading \code{group} column.
#' @keywords internal
#' @noRd
.mcc_single_group <- function(id, time, status, Tstart, grp_label, ci, conf.level, niter, seed) {
  if (ci) {
    res_list <- scumi_ci(
      id = id, time = time, status = status, Tstart = Tstart,
      niter = niter, conf.level = conf.level, seed = seed # conf.level passed here
    )
    res <- res_list[["MCC.output"]]
  } else {
    res <- scumi(id = id, time = time, status = status, Tstart = Tstart)
  }
  dplyr::mutate(res, group = grp_label, .before = dplyr::everything())
}
