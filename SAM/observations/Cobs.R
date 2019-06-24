library(TMB)
load("Cobs.RData")

compile("Cobs.cpp")
dyn.load(dynlib("Cobs"))


data <- list(Cobs = Cobs$Cobs, 
             N = Cobs$N,
             F = Cobs$F,
             M = Cobs$M, 
             aux = Cobs$aux,
             minAge = Cobs$minAge,
             minYear = Cobs$minYear)

param <- list(log_sigma = 0)

obj <- MakeADFun(data, param, DLL = "Cobs")
opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- sdreport(obj)
srep <- summary(rep)
logPred <- srep[rownames(srep) == "logPred", 1]

matplot(rownames(Cobs$N), xtabs(log(Cobs$Cobs)~Cobs$aux[, 1] + Cobs$aux[, 3]), ylab="Log C", xlab="Year")
matplot(rownames(Cobs$N), xtabs(logPred~Cobs$aux[,1]+Cobs$aux[,3]), type="l", add=TRUE)


# Plot year 4

Cobs_true <-  xtabs(log(Cobs$Cobs)~Cobs$aux[, 1] + Cobs$aux[, 3])
Cobs_estimated <-  xtabs(logPred~Cobs$aux[,1]+Cobs$aux[,3])

plot(rownames(Cobs_true), Cobs_true[, 4], type = "l")
lines(rownames(Cobs_estimated), Cobs_estimated[, 4], type = "l", col = "red")
