#' Internal helper to create framework for dataset
#' @importFrom utils tail
#' @noRd
create.wData.omega <- function(Tstart, Tstop, status, id, stratum, failcode, cens){
  event.times <- sort(unique(Tstop[status==failcode]))
  n <- length(Tstop)
  Nw <- rep(1,n)
  sel.compet <- (status!=cens)&(status!=failcode)&!is.na(status)

  Nw[sel.compet] <-  apply(outer(Tstop[sel.compet],event.times,"<"), 1, sum)+1
  data.weight <- data.frame(id=rep(id,Nw), Tstart=NA, Tstop=NA, status=rep(status,Nw),
                            strata=rep(stratum,Nw))
  data.weight$Tstart <- unlist(lapply(1:n, FUN=function(x,tms,N) {
    if (N[x]==1) {
      Tstart[x]
    } else {
      if (N[x]==2) {
        c(Tstart[x],Tstop[x])
      } else {
        c(Tstart[x], Tstop[x], rev(rev(tms)[2:(N[x]-1)]))
      }
    }
  }, tms=event.times, N=Nw))

  data.weight$Tstop <- unlist(lapply(1:n, FUN=function(x,tms,N) {
    if (N[x]==1) {
      Tstop[x]
    } else {
      c(Tstop[x], tail(tms,N[x]-1))
    }
  }, tms=event.times, N=Nw))

  return(data.weight)
}


#' Internal Function to create weighted data set for competing risks analyses
#'
#' This is a customized, internal version of crprep provided by Dr. Geskus on May 8, 2025
#'
#' @importFrom survival survfit Surv
#' @noRd
crprep2 <- function(Tstop, status, data, trans=1, cens=0, Tstart=0, id, strata,
                    keep, shorten=TRUE, rm.na=TRUE, origin=0,
                    prec.factor=1000, ...) {

  if (!missing(data)) data <- as.data.frame(data)

  if (!(is.numeric(Tstop))) {
    if (!is.character(Tstop) | length(Tstop)!=1) stop("argument \"Tstop\" should be a numeric vector or character string")
    if (missing(data)) stop("missing \"data\" argument not allowed when \"Tstop\" is character")
    tcol <- match(Tstop, names(data))
    if (is.na(tcol)) stop("\"Tstop\" not found in data")
    Tstop <- data[, tcol]
  }
  nn <- length(Tstop)

  if (is.numeric(Tstart)&length(Tstart) == 1) {
    Tstart <- rep(Tstart, nn)
  } else {
    if (!(is.numeric(Tstart))) {
      if (!is.character(Tstart) | length(Tstart)!=1) stop("argument \"Tstart\" should be a numeric vector or character string")
      if (missing(data)) stop("missing \"data\" argument not allowed when \"Tstart\" is character")
      tcol <- match(Tstart, names(data))
      if (is.na(tcol)) stop("\"Tstart\" not found in data")
      Tstart <- data[, tcol]
    }
  }
  if (length(Tstart) != nn) stop("Tstop and Tstart have different lengths")

  calc.trunc <- any(Tstart[!is.na(Tstart)] != 0)

  sel <- !is.na(Tstart) & !is.na(Tstop)
  if (any(Tstart[sel] >= Tstop[sel])) stop("Tstop must be greater than Tstart")

  if (length(status) == 1) {
    if (!is.character(status)) stop("argument \"status\" should be a vector or character string")
    if (missing(data)) stop("missing \"data\" argument not allowed when \"status\" is character")
    tcol <- match(status, names(data))
    if (is.na(tcol)) stop("\"status\" not found in data")
    status <- data[ ,tcol]
  }
  if (length(status) != nn) stop("Tstop and status have different lengths")

  if (missing(strata)) {
    strata.val <- rep(1,nn)
  } else {
    if (is.matrix(strata) | is.data.frame(strata)) stop("only one variable allowed in \"strata\"")
    if (!(is.vector(as.numeric(factor(strata))) & length(strata) > 1)) {
      if (!is.character(strata)) stop("argument \"strata\" should be a character string")
      if (missing(data)) stop("missing \"data\" argument not allowed when \"strata\" is character")
      tcol <- match(strata, names(data))
      if (is.na(tcol)) stop("\"strata\" not found in data")
      strata.name <- strata
      strata.val <- data[ ,tcol]
    } else {
      if (length(strata) != nn) stop("Tstop and strata have different lengths")
      strata.name <- names(strata)
      strata.val <- strata
    }
  }
  strata.num <- as.numeric(factor(strata.val))

  if (missing(id)) {
    id.name <- "id"
    id  <- num.id <- 1:nn
  } else {
    if (is.matrix(id) | is.data.frame(id)) stop("only one variable allowed in \"id\"")
    if (!(is.vector(id) & length(id) > 1)) {
      if (!is.character(id)) stop("argument \"id\" should be a character string")
      if (missing(data)) stop("missing \"data\" argument not allowed when \"id\" is character")
      tcol <- match(id, names(data))
      if (is.na(tcol)) stop("\"id\" not found in data")
      id.name <- id
      num.id <- 1:nn
      id <- data[, tcol]
    } else {
      if (length(id) != nn) stop("Tstop and id have different lengths")
      id.name <- names(id)
      num.id <- 1:nn
    }
  }

  if(rm.na) sel <- sel & !is.na(status)
  Tstart <- Tstart[sel]
  Tstop <- Tstop[sel]
  status <- status[sel]
  strata.val <- strata.val[sel]
  strata.num <- strata.num[sel]
  id <- id[sel]
  num.id <- num.id[sel]
  n <- length(Tstop)

  if(!missing(keep)) {
    if (!(is.matrix(keep) | is.data.frame(keep))) {
      if (is.character(keep)) {
        if (missing(data)) stop("missing \"data\" argument not allowed when \"keep\" is character")
        nkeep <- length(keep)
        kcols <- match(keep, names(data))
        if (any(is.na(kcols))) stop("at least one element of \"keep\" not found in data")
        keep.name <- keep
        keep <- data[, kcols]
      } else {
        nkeep <- 1
        keep.name <- names(keep)
        if (length(keep) != nn) stop("Tstop and keep have different lengths")
      }
    } else {
      nkeep <- ncol(keep)
      if(is.data.frame(keep)){
        keep.name <- names(keep)
      } else {
        keep.name <- colnames(keep)
        if(is.null(keep.name)) keep.name <- paste("V",1:nkeep,sep="")
      }
      if (nrow(keep) != nn) stop("length Tstop and number of rows in keep are different")
      if (nkeep == 1) keep <- keep[, 1]
    }
  }

  Tstart <- Tstart - origin
  Tstop <- Tstop - origin

  prec <- .Machine$double.eps*prec.factor
  Tstop.tie <- ifelse(status==cens, Tstop+prec, Tstop)
  Tstart.tie <- ifelse(Tstart==0, Tstart, Tstart +2* prec)

  surv.cens <- survival::survfit(survival::Surv(Tstart.tie, Tstop.tie, status == cens) ~ strata.num, timefix=FALSE)
  if(calc.trunc) surv.trunc <- survival::survfit(survival::Surv(-Tstop, -Tstart.tie, rep(1, n)) ~ strata.num, timefix=FALSE)

  data.out <- vector("list",length(trans))
  i.list <- 1
  strat <- sort(unique(strata.num),na.last=TRUE)
  len.strat <- length(strat)

  for(failcode in trans) {
    data.weight <- create.wData.omega(Tstart, Tstop, status, num.id, strata.num, failcode, cens)
    if(len.strat==1){
      tmp.time <- data.weight$Tstop
      data.weight$weight.cens <- summary(surv.cens, times=tmp.time)$surv
      if(calc.trunc) data.weight$weight.trunc <- summary(surv.trunc, times=-tmp.time)$surv
    } else {
      data.weight <- split(data.weight,data.weight$strata)
      if(is.na(strat[len.strat])) {
        tmp.sel <- is.na(strata.num)
        data.weight[[len.strat]] <- data.frame(id=num.id[tmp.sel], Tstart=Tstart[tmp.sel], Tstop=Tstop[tmp.sel], status=status[tmp.sel], strata=NA,  weight.cens=NA)
        if(calc.trunc) data.weight[[len.strat]]$weight.trunc <- NA
      }
      for(tmp.strat in 1:(len.strat-is.na(strat[len.strat]))){
        tmp.sel <- !is.na(strata.num) & strata.num==tmp.strat
        tmp.time <- data.weight[[tmp.strat]]$Tstop
        data.weight[[tmp.strat]]$weight.cens <- summary(surv.cens[tmp.strat], times=tmp.time)$surv
        if(calc.trunc) data.weight[[tmp.strat]]$weight.trunc <- summary(surv.trunc[tmp.strat], times=-tmp.time)$surv
      }
      data.weight <- do.call("rbind", data.weight)
    }

    data.weight <- data.weight[order(data.weight$id,data.weight$Tstop), ]

    data.weight$weight.cens <- unlist(tapply(data.weight$weight.cens, data.weight$id,
                                             FUN = function(x) {
                                               if (length(x) == 1 & !is.na(x[1])) {
                                                 return(1)
                                               } else if (length(x) != 1 & !is.na(x[1]) & x[1] != 0) {
                                                 return(x / x[1])
                                               } else {
                                                 return(rep(NA, length(x)))
                                               }
                                             })) #Grace Zhou debug to handle x[1]=0
    if(calc.trunc) {
      data.weight$weight.trunc <- unlist(tapply(data.weight$weight.trunc, data.weight$id,
                                                FUN = function(x) {
                                                  if (length(x) == 1 & !is.na(x[1])) {
                                                    return(1)
                                                  } else if (length(x) != 1 & !is.na(x[1]) & x[1] != 0) {
                                                    return(x / x[1])
                                                  } else {
                                                    return(rep(NA, length(x)))
                                                  }
                                                })) #Grace Zhou debug to handle x[1]=0
    }

    tbl <- table(data.weight$id)

    if(!missing(keep)) {
      if (is.null(keep.name)) {
        m <- match.call(expand.dots = FALSE)
        m <- m[match("keep", names(m))]
        if(!is.null(m)) {
          keep.name <- as.character(m[1])
          keep.name.split <- strsplit(keep.name, '')[[1]]
          tag <- which(keep.name.split == '$')
          if(length(tag) != 0) {
            keep.name <- substring(keep.name, tag[length(tag)]+1)
          } else {
            tag <- which(keep.name.split == '"')
            if(length(tag) != 0) {
              keep.name <- substring(keep.name, tag[1]+1, tag[2]-1)
            }
          }
        }
      }

      if (nkeep > 0) {
        if (nkeep == 1) {
          keep <- keep[sel]
          ddcovs <- rep(keep, tbl)
          ddcovs <- as.data.frame(ddcovs)
          names(ddcovs) <- as.character(keep.name)
        } else {
          keep <- keep[sel, ]
          ddcovs <- lapply(1:nkeep, function(i) rep(keep[, i], tbl))
          ddcovs <- as.data.frame(ddcovs)
          names(ddcovs) <- keep.name
        }
        data.weight <- cbind(data.weight, ddcovs)
      }
    }

    if (shorten) {
      if(calc.trunc) {
        keep.rows <- with(data.weight, c(diff(id)!=0 | diff(weight.cens)!=0 | diff(weight.trunc)!=0, TRUE))
      } else {
        keep.rows <- with(data.weight, c(diff(id)!=0 | diff(weight.cens)!=0, TRUE))
      }
      keep.rows[unlist(mapply(seq,1,tbl))==1] <- TRUE

      # GZ debug: Handles NA values generated by diff(NA) to prevent row collapse
      keep.rows[which(is.na(keep.rows))] <- TRUE

      keep.start <- data.weight$Tstart[unlist(tapply(keep.rows, data.weight$id, FUN=function(x) if(length(x)==1) x else c(TRUE,x[-length(x)])))]
      data.weight <- data.weight[keep.rows,]
      data.weight$Tstart <- keep.start
    }

    tbl <- table(data.weight$id)
    data.weight$count <- unlist(mapply(seq,1,tbl))
    data.weight$failcode <- failcode
    data.weight$id <- rep(id,tbl)

    data.out[[i.list]] <- data.weight
    i.list <- i.list+1
  }

  out <- do.call("rbind", data.out)

  if(!missing(strata)) {
    if (is.null(strata.name)) {
      m <- match.call(expand.dots = FALSE)
      m <- m[match("strata", names(m))]
      if(!is.null(m)) {
        strata.name <- as.character(m[1])
        strata.name.split <- strsplit(strata.name, '')[[1]]
        tag <- which(strata.name.split == '$')
        if(length(tag) != 0) {
          strata.name <- substring(strata.name, tag[length(tag)]+1)
        } else {
          tag <- which(strata.name.split == '"')
          if(length(tag) != 0) {
            strata.name <- substring(strata.name, tag[1]+1, tag[2]-1)
          }
        }
      }
    }
    if(is.factor(strata.val)) {
      out$strata <- factor(out$strata, labels=levels(strata.val))
    } else {
      out$strata <- levels(factor(strata.val))[out$strata]
      if(is.numeric(strata.val)) out$strata <- as.numeric(out$strata)
    }
    tmp.sel <- match("strata", names(out))
    names(out)[tmp.sel] <- strata.name
  } else {
    out$strata <- NULL
  }

  if (is.null(id.name)) {
    m <- match.call(expand.dots = FALSE)
    m <- m[match("id", names(m))]
    if(!is.null(m)) {
      id.name <- as.character(m[1])
      id.name.split <- strsplit(id.name, '')[[1]]
      tag <- which(id.name.split == '$')
      if(length(tag) != 0) {
        id.name <- substring(id.name, tag[length(tag)]+1)
      } else {
        tag <- which(id.name.split == '"')
        if(length(tag) != 0) {
          id.name <- substring(id.name, tag[1]+1, tag[2]-1)
        }
      }
    }
  }

  row.names(out) <- as.character(1:nrow(out))
  names(out)[1] <- id.name

  return(out)
}
