library(TMB)
compile("example_mean.cpp")
dyn.load(dynlib("example_mean"))

n <- 50
mu <- 10 

x <- rnorm(n, mu, 1)

param <- list(mu = 5)

data <- list(x = x, 
             n = n)

obj <- MakeADFun(data = data, parameters = param, DLL = "example_mean")
opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- sdreport(obj)

summary(rep)
