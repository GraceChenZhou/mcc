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

mydata <- data.frame(
  id = c(1, 2, 3, 4, 4, 4, 5, 5),
  tstop = c(8, 1, 5, 2, 6, 7, 3, 4),
  status = c(0, 0, 2, 1, 1, 1, 1, 2), # 1 = event, 0 = censored, 2 = competing risk
  tstart = c(0, 0, 0, 0, 2, 6, 0, 3)  # left-truncation times
)

mcc(id = mydata$id,
    time = mydata$tstop,
    status = mydata$status,
    Tstart = 0)

mcc(id = mydata$id,
    time = mydata$tstop,
    status = mydata$status)

my_package_rst <- mcc(id = mydata$id,
    time = mydata$tstop,
    status = mydata$status,
    Tstart = mydata$tstart)

my_package_rst_ci <- mcc(id = mydata$id,
   time = mydata$tstop,
   status = mydata$status,
   Tstart = mydata$tstart,
   ci = TRUE, niter = 100, seed = 123)

dir <- "C:/Users/gzhou/OneDrive - St. Jude Children's Research Hospital/4GRACE/MCC/MCC R PACKAGE/"
source(paste0(dir, 'CRPREP_20250508.R')) # Load MY.MCC.GESKUS function

geskus_rst <- MY.MCC.GESKUS(id = mydata$id,
              time = mydata$tstop,
              status = mydata$status,
              Tstart = mydata$tstart,
              ci=FALSE)

geskus_rst_ci <- MY.MCC.GESKUS(id = mydata$id,
              time = mydata$tstop,
              status = mydata$status,
              Tstart = mydata$tstart,
              ci=TRUE, niter=100, seed=123)

all.equal(my_package_rst, geskus_rst) #Class mismatch is acceptable
all.equal(my_package_rst_ci, geskus_rst_ci) #Class mismatch is acceptable

#Test 2: Correct result? ====

install.packages('mstate') #CRAN V0.3.3
install.packages('zoo')
library(mstate)
library(dplyr)
library(zoo)

rm(list = ls())

data <- data.frame(id=c(1,2,3,4), tstart=c(1,0,2,1), tstop=c(2,3,4,5), status=c(2,0,1,1))

mstate::crprep(Tstop=data$tstop,status=data$status,id=data$id,Tstart=data$tstart,shorten=FALSE) #NOT CORRECT

dir <- "C:/Users/gzhou/OneDrive - St. Jude Children's Research Hospital/4GRACE/MCC/MCC R PACKAGE/"
source(paste0(dir, 'CRPREP_20250508.R')) # Load MY.MCC.GESKUS function

crprep(Tstop=data$tstop,status=data$status,id=data$id,Tstart=data$tstart,shorten=FALSE) #CORRECT

source('R/crprep2.R')

crprep2(Tstop=data$tstop,status=data$status,id=data$id,Tstart=data$tstart,shorten=FALSE) #CORRECT

mcc_rst <- mcc(id = data$id,
            time = data$tstop,
            status = data$status,
            Tstart = data$tstart,
            ci=FALSE)

mcc_rst

# Test 3: Crash data  ====

data <- data.frame(id=c(1,2,3,4), tstart=c(1,0,2,1), tstop=c(1.001,1,2.001,5), status=c(2,0,1,1))

dir <- "C:/Users/gzhou/OneDrive - St. Jude Children's Research Hospital/4GRACE/MCC/MCC R PACKAGE/"
source(paste0(dir, 'CRPREP_20250508.R')) # Load MY.MCC.GESKUS function

crprep(Tstop=data$tstop,status=data$status,id=data$id,Tstart=data$tstart, shorten=FALSE) #weight.cens NaN
crprep(Tstop=data$tstop,status=data$status,id=data$id,Tstart=data$tstart, shorten=TRUE) #crash

MY.MCC.GESKUS(id = data$id,
              time = data$tstop,
              status = data$status,
              Tstart = data$tstart,
              ci=FALSE) #crash

crprep2(Tstop=data$tstop,status=data$status,id=data$id,Tstart=data$tstart,shorten=FALSE) #weight.cens NA
crprep2(Tstop=data$tstop,status=data$status,id=data$id,Tstart=data$tstart,shorten=TRUE) #weight.cens NA

mcc_rst <- mcc(id = data$id,
               time = data$tstop,
               status = data$status,
               Tstart = data$tstart,
               ci=FALSE) #not crashed

mcc_rst

mcc(id = data$id,
    time = data$tstop,
    status = data$status,
    Tstart = data$tstart,
    ci=TRUE, niter=100, seed=123)

# Test 4: Consistency across Version 4.4.0 - 4.5.2 ====

## Step 1: Initialize the Testing Infrastructure
usethis::use_testthat()

## Create a Test File for MCC
usethis::use_test("mcc")

# Test 5: Survival package ====

# Define the historical versions to check
versions_to_test <- c("3.5-8", "3.7-0")

  # 1. Create a temporary folder to act as a sandbox library
  my_temp_lib <- file.path(tempdir(), "old_survival")
  dir.create(my_temp_lib, showWarnings = FALSE)

  # 2. Tell remotes to install the old version ONLY into this new folder
  # (I turned quiet = FALSE so you can see if it throws any C++ compilation errors)
  remotes::install_version("survival", version = "3.7-0", lib = my_temp_lib, quiet = FALSE)

  # 3. Unload your current modern survival package
  if ("package:survival" %in% search()) detach("package:survival", unload=TRUE)

  try(unloadNamespace("mstate"), silent = TRUE)
  try(unloadNamespace("survival"), silent = TRUE)
  .libPaths(c(my_temp_lib, .libPaths()))
  # 4. Load the historical version explicitly from your sandbox
  library(survival, lib.loc = my_temp_lib)

  # 5. Verify it worked
  packageVersion("survival")

  # Run your exact test
  tmp <- as.data.frame(matrix(c(1,2,6,1,3,1,3,2,4,4,6,1,5,5,7,1),byrow=TRUE,nrow=4))
  names(tmp) <- c("id","tstart","tstop","status")
  fit <- survfit(Surv(-tmp$tstop, -tmp$tstart, c(1,1,1,1)) ~ 1)
  res <- summary(fit, times = c(-5, -6, -3, -6, -7, -6, -7))

  print(res)

# Remember to reinstall the latest version when you are done!
# install.packages("survival")
