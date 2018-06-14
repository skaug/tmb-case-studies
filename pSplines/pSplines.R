library(TMB)
library(mgcv) #Use gam
library(Matrix) #Use sparse matrices

# Tutorial: sparse matrixes in R and TMB
M = matrix(1:4,2,2)   # ordinary 2x2 matrix
M_block_diag = .bdiag(list(M,M)) # Block diagonal (sparse) matrix
data.class(M_block_diag)     # Check data.class
print(M_block_diag)    # dots means 0 value

#Load the data and compile c++ code------
Vegetation <- read.table(file = "Vegetation.txt", header = TRUE, dec = ".")
Vegetation = Vegetation[!is.na(Vegetation$Richness),]

compile("pSplines.cpp")
dyn.load(dynlib("pSplines"))
#----------------------------------------

#Set up spline structure by using mgcv---
gam_setup = gam(Richness ~ s(ROCK, bs = "cs") +
      s(LITTER, bs = "cs") + s(BARESOIL, bs = "cs") +
      s(FallPrec, bs = "cs") + s(SprTmax, bs = "cs"),
    data = Vegetation,fit=FALSE)

#Extrtact penelization matrices
S_ROCK = gam_setup$smooth[[1]]$S[[1]]
S_LITTER = gam_setup$smooth[[2]]$S[[1]]
S_BARESOIL = gam_setup$smooth[[3]]$S[[1]]
S_FallPrec = gam_setup$smooth[[4]]$S[[1]]
S_SprTmax = gam_setup$smooth[[5]]$S[[1]]

S_list = list(S_ROCK,S_LITTER,S_BARESOIL,S_FallPrec,S_SprTmax)
S_combined = .bdiag(S_list)         # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S
#----------------------------------------


#For report, used for constructing plots----
ROCK=seq(min(Vegetation$ROCK),max(Vegetation$ROCK),by = 1)
LITTER=seq(min(Vegetation$LITTER),max(Vegetation$LITTER),by = 1)
BARESOIL=seq(min(Vegetation$BARESOIL),max(Vegetation$BARESOIL),by = 1)
FallPrec=seq(min(Vegetation$FallPrec),max(Vegetation$FallPrec),by = 0.2)
SprTmax=seq(min(Vegetation$SprTmax),max(Vegetation$SprTmax),by = 0.2)

rockReport = PredictMat(gam_setup$smooth[[1]],data = data.frame(ROCK))
litterReport = PredictMat(gam_setup$smooth[[2]],data = data.frame(LITTER))
soilReport = PredictMat(gam_setup$smooth[[3]],data = data.frame(BARESOIL))
fallReport = PredictMat(gam_setup$smooth[[4]],data = data.frame(FallPrec))
sprReport = PredictMat(gam_setup$smooth[[5]],data = data.frame(SprTmax))

designMatrixForReport = list(rockReport,litterReport,soilReport,fallReport,sprReport)
#-------------------------------------------

#Define data object which is given to TMB---
data = list(Y = Vegetation$Richness, # Response
            X = gam_setup$X[,-1],  # Design matrix, without intercept
            S = S_combined,      # Combined penalty matrix
            Sdims = Sdims,
            designMatrixForReport = .bdiag(designMatrixForReport))
#-------------------------------------------

#Define parameter object given to TMB-------
par = list(
  beta0 = 0,  # Intercept
  beta = rep(0,sum(Sdims)),  # Spline coefficients
  log_lambda = rep(rep(0,length(Sdims))), #Log spline penalization coefficients
  log_sigma = 0
)
#-------------------------------------------

#Fit model----------------------------------
obj = MakeADFun(data = data, parameters = par,random="beta",DLL = "pSplines")
opt = nlminb(obj$par,obj$fn,obj$gr)
rep = sdreport(obj)
#-------------------------------------------


#Plot results
muSpline = rep$value[names(rep$value)=="splineForReport"]
sdSpline<-rep$sd[names(rep$value)=="splineForReport"]

par(mfrow=c(2,3))
start = 1
stop = start + length(ROCK) -1
plot(ROCK, muSpline[start:stop], lty=1,type = 'l',ylim = c(-6,5),ylab = "f(rock)",main = "Spline for ROCK")
lines(ROCK, muSpline[start:stop]- 1.96*sdSpline[start:stop], lty=2)
lines(ROCK, muSpline[start:stop]+ 1.96*sdSpline[start:stop], lty=2)
abline(h = 0)

start = stop +1
stop = start+ length(LITTER)-1
plot(LITTER, muSpline[start:stop], lty=1,type = 'l',ylim = c(-6,5),ylab = "f(litter)",main = "Spline for LITTER")
lines(LITTER, muSpline[start:stop]- 1.96*sdSpline[start:stop], lty=2)
lines(LITTER, muSpline[start:stop]+ 1.96*sdSpline[start:stop], lty=2)
abline(h = 0)

start = stop +1
stop = start+ length(BARESOIL)-1
plot(BARESOIL, muSpline[start:stop], lty=1,type = 'l',ylim = c(-6,5),ylab = "f(soil)",main = "Spline for BARESOIL")
lines(BARESOIL, muSpline[start:stop]- 1.96*sdSpline[start:stop], lty=2)
lines(BARESOIL, muSpline[start:stop]+ 1.96*sdSpline[start:stop], lty=2)
abline(h = 0)

start = stop +1
stop = start+ length(FallPrec)-1
plot(FallPrec, muSpline[start:stop], lty=1,type = 'l',ylim = c(-6,5),ylab = "f(fallPrec)",main = "Spline for FallPrec")
lines(FallPrec, muSpline[start:stop]- 1.96*sdSpline[start:stop], lty=2)
lines(FallPrec, muSpline[start:stop]+ 1.96*sdSpline[start:stop], lty=2)
abline(h = 0)

start = stop +1
stop = start+ length(SprTmax)-1
plot(SprTmax, muSpline[start:stop], lty=1,type = 'l',ylim = c(-6,5),ylab = "f(sprTMax)",main = "Spline for SprTmax")
lines(SprTmax, muSpline[start:stop]- 1.96*sdSpline[start:stop], lty=2)
lines(SprTmax, muSpline[start:stop]+ 1.96*sdSpline[start:stop], lty=2)
abline(h = 0)
