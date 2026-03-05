devtools::document()
devtools::check()

# 1. Ignore the .github folder (used for your R 4.4.x-4.5.2 testing matrix)
usethis::use_build_ignore(".github")

# 2. Ignore RStudio project file
usethis::use_build_ignore("github_mcc.Rproj")

# 3. Ignore code.r
usethis::use_build_ignore("code.R")

#Test 1: My pacakge vs. Function: Same result====

devtools::install_github("GraceChenZhou/mcc")

library(mcc)

## 1. Create a sample recurrent event dataset
mydata <- data.frame(
  id = c(1, 2, 3, 4, 4, 4, 5, 5),
  tstop = c(8, 1, 5, 2, 6, 7, 3, 4),
  status = c(0, 0, 2, 1, 1, 1, 1, 2), # 1 = event, 0 = censored, 2 = competing risk
  tstart = c(0, 0, 0, 0, 2, 6, 0, 3)  # left-truncation times
)

## 2. Calculate the Mean Cumulative Count
# Set ci = TRUE to calculate 95% bootstrap confidence intervals
mcc(id = mydata$id,
    time = mydata$tstop,
    status = mydata$status,
    Tstart = mydata$tstart)



result_ci <- mcc(id = mydata$id,
                 time = mydata$tstop,
                 status = mydata$status,
                 Tstart = mydata$tstart,
                 ci = TRUE, niter = 100, seed = 123)

print(result_ci)

source("C:/Users/gzhou/OneDrive - St. Jude Children's Research Hospital/4GRACE/MCC/MCC R PACKAGE/MY_MCC_R.R")
library(dplyr)
MY.MCC(id = mydata$id,
       time = mydata$tstop,
       status = mydata$status,
       Tstart = mydata$tstart,
       ci=FALSE)

MY.MCC(id = mydata$id,
       time = mydata$tstop,
       status = mydata$status,
       Tstart = mydata$tstart,
       ci=TRUE, niter = 100)

#Test 2: Correct result? ====

install.packages('mstate') #CRAN V0.3.3
install.packages('zoo')
library(mstate)
library(dplyr)
library(zoo)

data <- data.frame(id=c(1,2,3,4), tstart=c(1,0,2,1), tstop=c(2,3,4,5), status=c(2,0,1,1))

mstate::crprep(Tstop=data$tstop,status=data$status,id=data$id,Tstart=data$tstart,shorten=FALSE) #NOT CORRECT

SCI <- function(id, time, status, Tstart = 0) {
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
    MCC.out <- data.frame(time = as.numeric(ftime), MCC = rep(0, length(ftime)))
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
        w <- if ("weight.trunc" %in% names(count.data))
          count.data$weight.cens * count.data$weight.trunc
        else
          count.data$weight.cens
        fit.i.th <- survival::survfit(
          survival::Surv(Tstart, Tstop, status == 1) ~ 1,
          data = count.data, weight = w)
        output  <- summary(fit.i.th)
        output2 <- data.frame(time = output$time, CumI = 1 - output$surv)
        CumI.i.th <- data.frame(time = ftime) %>%
          dplyr::left_join(output2, by = "time") %>%
          dplyr::mutate(
            CumI = ifelse(dplyr::row_number() == 1 & is.na(.data$CumI), 0, .data$CumI),
            CumI = zoo::na.locf(.data$CumI),
            M    = i)
      } else {
        fit.i.th  <- cmprsk::cuminc(out.i.th$time, out.i.th$status)
        # FIX: consistent data.frame shape for both branches
        CumI.i.th <- data.frame(
          time = ftime,
          CumI = as.numeric(cmprsk::timepoints(fit.i.th, ftime)$est[1, ]),
          M    = i)
      }
      CumI.list[[i]] <- CumI.i.th
    }
    MCC.raw <- do.call(rbind, CumI.list)
    MCC.fill <- MCC.raw %>%
      tidyr::pivot_wider(names_from = "M", values_from = "CumI", names_prefix = "M")
    # FIX: unified path; na.locf applied column-wise before summing
    MCC.fill[, -1] <- lapply(MCC.fill[, -1], zoo::na.locf, na.rm = FALSE)
    MCC <- rowSums(MCC.fill[, -1], na.rm = TRUE)
    MCC.out <- data.frame(time = MCC.fill$time, MCC = MCC)
  }
  MCC.out <- MCC.out %>% dplyr::mutate_at("time", as.numeric)
  rownames(MCC.out) <- seq(nrow(MCC.out))
  return(MCC.out)
}

MY.MCC <- function(id, time, status, Tstart = 0, ci = TRUE, niter = 1000, seed = NULL) {
  if (ci) {
    return(SCI.95CI(id, time, status, Tstart, niter, seed)[["MCC.output"]])
  } else {
    return(SCI(id, time, status, Tstart))
  }
}

MY.MCC(id = data$id,
       time = data$tstop,
       status = data$status,
       Tstart = data$tstart,
       ci=FALSE)

dir <- "C:/Users/gzhou/OneDrive - St. Jude Children's Research Hospital/4GRACE/MCC/MCC R PACKAGE/"
source(paste0(dir, 'mstate_0.3.3/mstate/R/create.wData.omega.r'))
source(paste0(dir, 'mstate_0.3.3/mstate/R/crprep.r'))
crprep(Tstop=data$tstop,status=data$status,id=data$id,Tstart=data$tstart,shorten=FALSE) #CORRECT

SCI_GESKUS <- function(id, time, status, Tstart = 0) {
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
    MCC.out <- data.frame(time = as.numeric(ftime), MCC = rep(0, length(ftime)))
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
        count.data <- crprep(Tstop = "time", status = "status", data = out.i.th,
                                     trans = 1, cens = 0, Tstart = "tstart", id = "id")
        w <- if ("weight.trunc" %in% names(count.data))
          count.data$weight.cens * count.data$weight.trunc
        else
          count.data$weight.cens
        fit.i.th <- survival::survfit(
          survival::Surv(Tstart, Tstop, status == 1) ~ 1,
          data = count.data, weight = w)
        output  <- summary(fit.i.th)
        output2 <- data.frame(time = output$time, CumI = 1 - output$surv)
        CumI.i.th <- data.frame(time = ftime) %>%
          dplyr::left_join(output2, by = "time") %>%
          dplyr::mutate(
            CumI = ifelse(dplyr::row_number() == 1 & is.na(.data$CumI), 0, .data$CumI),
            CumI = zoo::na.locf(.data$CumI),
            M    = i)
      } else {
        fit.i.th  <- cmprsk::cuminc(out.i.th$time, out.i.th$status)
        # FIX: consistent data.frame shape for both branches
        CumI.i.th <- data.frame(
          time = ftime,
          CumI = as.numeric(cmprsk::timepoints(fit.i.th, ftime)$est[1, ]),
          M    = i)
      }
      CumI.list[[i]] <- CumI.i.th
    }
    MCC.raw <- do.call(rbind, CumI.list)
    MCC.fill <- MCC.raw %>%
      tidyr::pivot_wider(names_from = "M", values_from = "CumI", names_prefix = "M")
    # FIX: unified path; na.locf applied column-wise before summing
    MCC.fill[, -1] <- lapply(MCC.fill[, -1], zoo::na.locf, na.rm = FALSE)
    MCC <- rowSums(MCC.fill[, -1], na.rm = TRUE)
    MCC.out <- data.frame(time = MCC.fill$time, MCC = MCC)
  }
  MCC.out <- MCC.out %>% dplyr::mutate_at("time", as.numeric)
  rownames(MCC.out) <- seq(nrow(MCC.out))
  return(MCC.out)
}
MY.MCC.GESKUS <- function(id, time, status, Tstart = 0, ci = TRUE, niter = 1000, seed = NULL) {
  if (ci) {
    return(SCI.95CI(id, time, status, Tstart, niter, seed)[["MCC.output"]])
  } else {
    return(SCI_GESKUS(id, time, status, Tstart))
  }
}

MY.MCC.GESKUS(id = data$id,
       time = data$tstop,
       status = data$status,
       Tstart = data$tstart,
       ci=FALSE)

devtools::install_github("GraceChenZhou/mcc")
library(mcc)
data <- data.frame(id=c(1,2,3,4), tstart=c(1,0,2,1), tstop=c(2,3,4,5), status=c(2,0,1,1))

result <- mcc(id = data$id,
              time = data$stop,
              status = data$status,
              Tstart = data$tstart)

# Test 3: Real data analysis ====

# Test 4: Consistency across Version 4.4.0 - 4.5.2 ====
