# How to turn parameters on-and-off using the "map" argument
library(TMB)

compile("map_example.cpp")
dyn.load(dynlib("map_example"))

# Simulate artificial data from random effects model
set.seed(234)
n = 100         # Number of individuals
n_i = 4
n_j = 10
i = sample(1:n_i,size=n,replace=TRUE)   # Levels of factor 1
j = sample(1:n_j,size=n,replace=TRUE)   # Levels of factor 2
beta = rnorm(n_i,0,1)                     # Coeficients of factor 1
u = rnorm(n_j,0,1)                       # Random effects vector associated with 2

y = beta[i] + u[j] + rnorm(n,0,.2)         # Data

data <- list(y = y,
             i = i-1,                   # 0-indexing in C++
             j = j-1)

parameters <- list(beta = rep(0,n_i),
                   u = rep(0,n_j),
                   sigma = 1,
                   tau = 1)

# Fit full model
obj <- MakeADFun(data, parameters, random="u", DLL="map_example")
opt <- nlminb(obj$par, obj$fn, obj$gr)

# Make all beta's equal
map=list(beta=as.factor(c(1,1,1,1)))
obj <- MakeADFun(data, parameters, random="u", DLL="map_example",map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

# Group beta's in two groups
map=list(beta=as.factor(c(1,1,2,2)))
obj <- MakeADFun(data, parameters, random="u", DLL="map_example",map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

# Do not estimate beta[1] (fix at initial value 0)
map=list(beta=as.factor(c(NA,2,3,4)))
obj <- MakeADFun(data, parameters, random="u", DLL="map_example",map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

# Do not estimate any beta's (fix at initial value 0)
map=list(beta=as.factor(c(NA,NA,NA,NA)))
obj <- MakeADFun(data, parameters, random="u", DLL="map_example",map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

# Remove random effect from model
map=list(u=as.factor(rep(NA,length(u))),tau=as.factor(NA))
obj <- MakeADFun(data, parameters, random="u", DLL="map_example",map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

# Bounds on parameters (not "u")
U = c(rep(10,n_i),10,10)
names(U) = c(rep("beta",n_i),"sigma","tau")
L = c(rep(-10,4),1e-10,1e-10)
names(L) = names(U)
map=list()
obj <- MakeADFun(data, parameters, random="u", DLL="map_example",map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr,lower=L,upper=U)

# Model without random effect but with parameter bounds
map=list(u=as.factor(rep(NA,length(u))),tau=as.factor(NA))
obj <- MakeADFun(data, parameters, random="u", DLL="map_example",map=map)
print(names(U))  # To see which component of L and U should be removed
opt <- nlminb(obj$par, obj$fn, obj$gr,lower=L[-6],upper=U[-6])
