# Code from Thygesen et al., 2017 is used.
library(TMB)
compile("OSA.cpp")
dyn.load(dynlib("OSA"))

set.seed(123)
## Simulate data with these parameters
mu = 0.75
sigma = 1
s = 1
huge = 1e3
nT = 100
X = c(0,cumsum(rnorm(nT-1,mean=mu,sd=sigma)))
Y = X + rnorm(nT,sd=s)

data = list(y=Y,huge=huge)
parameters = list(x=X,mu=0,logsigma=log(sigma),logs=log(s))
obj = MakeADFun(data,parameters,random=c("x"),DLL="OSA")
opt = nlminb(obj$par, obj$fn,obj$gr)


### One-step predictions
predict  <- oneStepPredict(obj,observation.name="y",method="fullGaussian",data.term.indicator = "keep")
plot(predict$residual,xlab="Time", main = "OSA residuals", ylab = "Residual")
qqnorm(predict$residual)
qqline(predict$residual, col = "steelblue", lwd = 2)
