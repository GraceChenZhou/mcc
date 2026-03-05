# mcc: Mean Cumulative Count for Left-Truncated Competing Risks Data

[![R-CMD-check](https://github.com/GraceChenZhou/mcc/actions/workflows/r.yml/badge.svg)](https://github.com/GraceChenZhou/mcc/actions/workflows/r.yml)

The **mcc** R package provides a robust and efficient tool to estimate the total burden of recurrent events within a population using the Mean Cumulative Count (MCC) method. 

It is specifically designed to handle complex survival data, accurately applying time-dependent censoring and truncation weights for left-truncated data in the presence of competing risks. The methodology implements an optimized Geskus (2011) weighting approach, making it highly suitable for complex, long-term cohort studies.

# Installation

You can install the development version of the package directly from GitHub:

```r
# Install devtools if you haven't already
# install.packages("devtools")

devtools::install_github("GraceChenZhou/mcc")
```

# Dependencies

To maintain mathematical accuracy when computing weights for left‑truncated data, this package requires R version 4.4.0 or later and survival version 3.7-0 or later. 

# Basic Usage

The primary function is `mcc()`. Here is a quick, reproducible example showing how to calculate the Mean Cumulative Count for a sample cohort:

```r
library(mcc)

## 1. Create a sample recurrent event dataset
mydata <- data.frame(
  id = c(1, 2, 3, 4, 4, 4, 5, 5),
  time = c(8, 1, 5, 2, 6, 7, 3, 4),
  status = c(0, 0, 2, 1, 1, 1, 1, 2), # 1 = event, 0 = censored, 2 = competing risk
  tstart = c(0, 0, 0, 0, 2, 6, 0, 3)  # left-truncation times
)

## 2. Calculate the Mean Cumulative Count 
# Set ci = TRUE to calculate 95% bootstrap confidence intervals
result <- mcc(id = mydata$id, 
              time = mydata$time, 
              status = mydata$status, 
              Tstart = mydata$tstart)
                 
print(result)
```

