load("input.RData") #Data and code downloaded from https://github.com/fishfollower/samex
library(TMB)
compile("sam.cpp")
dyn.load(dynlib("sam"))

#Fit model
obj <- MakeADFun(data,parameters,random=c("logN","logF"),DLL="sam")
opt<-nlminb(obj$par,obj$fn,obj$gr,control=list(trace=1,eval.max=1200,iter.max=900))
sdrep<-sdreport(obj)


#Plot results
pl <- as.list(sdrep,"Est")
plsd <- as.list(sdrep,"Std")
plotit <-function (what, x=data$years, ylab=what, xlab="Years", trans=function(x)x ,...){
  idx<-names(sdrep$value)==what
  y<-sdrep$value[idx]
  ci<-y+sdrep$sd[idx]%o%c(-2,2)
  plot(x,trans(y), xlab=xlab, ylab=ylab, type="l", lwd=3, ylim=range(c(trans(ci),0)), las=1,...)
  polygon(c(x,rev(x)), y = c(trans(ci[,1]),rev(trans(ci[,2]))), border = gray(.5,alpha=.5), col = gray(.5,alpha=.5))
  grid(col="black")
}

par(mfrow=c(2,2))
options(scipen=-3)
plotit("logssb", ylab="SSB", trans=exp)
plotit("logfbar", ylab="Fbar", trans=exp)
plotit("logR", ylab="R", trans=exp)
plotit("logCatch", ylab="Catch", trans=exp, x=data$years[-length(data$years)])

