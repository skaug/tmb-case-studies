library(TMB)
compile("debug_tutorial.cpp",flags="-O0 -g",DLLFLAGS="",libtmb=FALSE)
# Note to Windows user: flags="-O1 -g" must be used for larger programs,
#                       but does not gives worse debugger functionality

dyn.load(dynlib("debug_tutorial"))

dat <- list(X=matrix(1:6,nrow=3,ncol=2),y=1:5)
pars <- list(a=0)
obj <- MakeADFun(data=dat, parameters=pars, DLL="debug_tutorial")
