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


#### Estimate SSB #### 


load("Nobs2.RData")

compile("survival_ssb.cpp")
dyn.load(dynlib("survival_ssb"))


data <- list(Nobs = Nobs$Nobs, 
             M = Nobs$M, 
             F = Nobs$F,
             SW = Nobs$SW,
             MO = Nobs$MO,
             PF = Nobs$PF,
             PM = Nobs$PM,
             method = 1)

param <- list(logN = matrix(0 ,nrow = nrow(Nobs$Nobs), ncol = ncol(Nobs$Nobs)),
              log_sigma_Nobs = 0, 
              log_sigma_logN = 0,
              log_sigma_logR = 0,
              ricker = if(data$method == 1){c(1,1)} else{numeric(0)}, # Ricker is sensitive to inital conditions. Fail to converge when it starts in c(0,0)
              bh = if(data$method == 2){numeric(2)} else{numeric(0)})


obj <- MakeADFun(data, param, random = "logN", DLL = "survival_ssb")
opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- sdreport(obj)

# Get estimated logN and ssb 
matplot(as.list(rep, "Est")$logN, type = "l", add = TRUE)

# Plot ssb 
srep <- summary(rep)
ssb <- srep[rownames(srep) == "ssb", ]

ts.plot(ssb[, 1], col = "red")
lines(ssb[, 1] + ssb[, 2], col = "red", lty = 2)
lines(ssb[, 1] - ssb[, 2], col = "red", lty = 2)


# Random walk = 1299
# Ricker nll = 1318
# Beverton-Holt = 1293
