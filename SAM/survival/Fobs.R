library(TMB)
compile("Fobs.cpp")
dyn.load(dynlib("Fobs"))


load("Fobs.Rdata")
matplot(Fobs$year, log(Fobs$Fobs), xlab="Year", ylab="logF", pch=colnames(Fobs$Fobs))

data <- list(Fobs = Fobs$Fobs,
             cormode = 2)

param <- list(logF = matrix(0, nrow = nrow(Fobs$Fobs), ncol = ncol(Fobs$Fobs)),
              log_sigma_logF = rep(0, ncol(Fobs$Fobs)),
              trans_rho = if(data$cormode < 2){numeric(0)} else{10},
              log_sigma_Fobs = 0)

# Important: trans_rho has to be a parameter vector of size 1 to be able to use numeric()
map <- list(trans_rho = factor(NA))
obj <- MakeADFun(data, param, random = "logF", DLL = "Fobs", map = map)#, map = list(log_sigma_logF = factor(rep(1, ncol(Fobs$Fobs)))))
opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- sdreport(obj)

matplot(Fobs$year, as.list(rep, "Est")$logF, type="l", add=TRUE)


Sigma <- rep$value[names(rep$value) == "Sigma"]
