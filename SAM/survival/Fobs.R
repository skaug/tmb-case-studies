library(TMB)
compile("Fobs.cpp")
dyn.load(dynlib("Fobs"))


load("Fobs.Rdata")
matplot(Fobs$year, log(Fobs$Fobs), xlab="Year", ylab="logF", pch=colnames(Fobs$Fobs))

data <- list(Fobs = Fobs$Fobs,
             cormode = 1)

param <- list(logF = matrix(0, nrow = nrow(Fobs$Fobs), ncol = ncol(Fobs$Fobs)),
              log_sigma_logF = rep(0, ncol(Fobs$Fobs)),
              trans_rho = if(data$cormode < 2){numeric(0)} else{0.1},
              log_sigma_Fobs = 0)

obj <- MakeADFun(data, param, random = "logF", DLL = "Fobs")
opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- sdreport(obj)

matplot(Fobs$year, as.list(rep, "Est")$logF, type="l", add=TRUE)
