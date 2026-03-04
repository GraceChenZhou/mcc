devtools::document()
devtools::check()
# 1. Ignore the .github folder (used for your R 4.4.x-4.5.2 testing matrix)
usethis::use_build_ignore(".github")

# 2. Ignore your RStudio project file
usethis::use_build_ignore("github_mcc.Rproj")

usethis::use_build_ignore("code.R")

#Test

devtools::install_github("GraceChenZhou/mcc")

devtools::install_github("GraceChenZhou/mcc", force=TRUE)

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

result_ci <- mcc(id = mydata$id,
                 time = mydata$time,
                 status = mydata$status,
                 Tstart = mydata$tstart,
                 ci = TRUE, niter = 500)

print(result_ci)

source("C:/Users/gzhou/OneDrive - St. Jude Children's Research Hospital/4GRACE/MCC/MCC R PACKAGE/MY_MCC_R.R")
library(dplyr)
MY.MCC(id = mydata$id,
       time = mydata$time,
       status = mydata$status,
       Tstart = mydata$tstart,
       ci=FALSE)

MY.MCC(id = mydata$id,
       time = mydata$time,
       status = mydata$status,
       Tstart = mydata$tstart,
       ci=TRUE, niter = 500)
