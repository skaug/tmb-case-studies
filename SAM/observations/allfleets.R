library(TMB)
load("allfleets.RData")
compile("allfleets.cpp")
dyn.load(dynlib("allfleets"))

allfleets$keyQ <- rbind(rep(NA, 9), 
                        c(NA, 0, 1, 2, 3, 4, 5, 6, NA),
                        c(7, 8, 9, 10, 11, 12, NA, NA, NA))

# keyQ starts from 0 and adds 1 for each age class that is observed. This continued across surveys. We have one standard deviation per age group pr fleet

allfleets$keySd <- rbind(rep(0, 9), 
                            c(NA, 1, 1, 1, 1, 1, 1, 1, NA),
                            c(2, 2, 2, 2, 2, 2, NA, NA, NA))
# keyFleet is an indicator of which fleet we are observing - proportionality coefficient

param <- list(logQ = numeric(max(allfleets$keyQ, na.rm = T) + 1), 
              log_sigma_s = numeric(max(allfleets$keySd, na.rm = T) + 1),
              missing = numeric(sum(is.na(allfleets$obs))))


obj <- MakeADFun(allfleets, param, random = "missing", DLL = "allfleets")
opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- summary(sdreport(obj), select = "random")
est <- rep[, 1]
est <- obj$report()$logPred

par(mfrow=c(1,1))
for(f in 1:3){
  idx<-which(allfleets$aux[,2]==f)
  matplot(xtabs(log(allfleets$obs[idx])~allfleets$aux[idx,1]+allfleets$aux[idx,3], na.action = na.pass), ylab="Log Obs")
  matplot(xtabs(est[idx]~allfleets$aux[idx,1]+allfleets$aux[idx,3]), type="l", add=TRUE)
}
