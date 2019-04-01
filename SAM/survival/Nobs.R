library(TMB)
load("Nobs.RData")
matplot(log(Nobs$Nobs), main="logN")

compile("survival.cpp")
dyn.load(dynlib("survival"))

# Prepare for TMB

param <- list(logN = matrix(0 ,nrow = nrow(Nobs$Nobs), ncol = ncol(Nobs$Nobs)),
              log_sigma_Nobs = 0, 
              log_sigma_logN = 0,
              log_sigma_logR = 0)

data <- list(Nobs = Nobs$Nobs,
             F = Nobs$F,
             M = Nobs$M)

obj <- MakeADFun(data, param, random = "logN", DLL = "survival")
opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- sdreport(obj)

logN_rep <- summary(rep, "random")

matplot(as.list(rep ,"Est")$logN, type="l", add=TRUE)
