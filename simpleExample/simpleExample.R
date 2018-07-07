library(TMB)

compile("simpleExample.cpp")
dyn.load(dynlib("simpleExample"))

TN = read.table("TeethNitrogen.txt",header = T)

data = list()
data$NL = TN$X15N[TN$Tooth=="Moby"]
data$x = TN$Age[TN$Tooth=="Moby"]

parameters = list(
  beta0 = 0,
  beta1 = 0,
  logSigma = 0
)

obj = MakeADFun(data,parameters,DLL = "simpleExample")
opt = nlminb(obj$par,obj$fn, obj$gr)
rep = sdreport(obj)

plot(data$x,data$NL)
abline(a = rep$par.fixed[1], b = rep$par.fixed[2])

