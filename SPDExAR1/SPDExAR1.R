library(INLA)
library(TMB)
library(fields) #use image.plot()
library(animation) #use saveGIF()

#Comile and load c++ code-------
compile("SPDExAR1.cpp")
dyn.load(dynlib("SPDExAR1"))
#--------------------------------

#Read data------------
borders <-read.table("Piemonte_borders.csv",header=TRUE,sep=",")
Piemonte_data <-read.table("Piemonte_data_byday.csv",header=TRUE,sep=",")
coordinates <-read.table("coordinates.csv",header=TRUE,sep=",")
rownames(coordinates) = coordinates[,"Station.ID"]
#--------------------------------

#Define some variables used------
n_stations <- length(coordinates$Station.ID)
n_data <- length(Piemonte_data$Station.ID)
n_days <- as.integer(n_data/n_stations)
#--------------------------------

##Standardize covariates and calculate log of PM10--------
mean_covariates = apply(Piemonte_data[,3:10],2,mean)
sd_covariates = apply(Piemonte_data[,3:10],2,sd)
Piemonte_data[,3:10] =
  scale(Piemonte_data[,3:10],
        mean_covariates, sd_covariates)
Piemonte_data$logPM10 = log(Piemonte_data$PM10)
Piemonte_data$logPM10[which(is.na(Piemonte_data$logPM10))]=-999
#---------------------------------------------------------

#Calcualtes the days from first observation---------------
Piemonte_data$day = rep(0,n_data)
for(i in 1:n_days){
  for(j in 1:n_stations){
    Piemonte_data$day[(i-1)*n_stations + j] = i-1
  }
}
#---------------------------------------------------------

#Group the days into intervalls used as time steps in the AR1 process--------
dtLength = 21#Length of time intervalls
Piemonte_data$dt = floor(Piemonte_data$day/dtLength)
lengthDt = sum(Piemonte_data$dt==0)#Number of observations within each time intervall
maxDt = max(Piemonte_data$dt)+1
#---------------------------------------------------------

#Defines the GMRF and represenation of the sparce precision matrix------
locStations = cbind(coordinates$UTMX,
            coordinates$UTMY)
mesh =
  inla.mesh.2d(loc=locStations,
               loc.domain=borders,
               offset=c(10, 140),
               max.edge=c(40, 80),
               min.angle=c(26, 21),
               cutoff=0
  )

spde = inla.spde2.matern(mesh=mesh, alpha=2)
spdeMatrices = spde$param.inla[c("M0","M1","M2")]
A = inla.spde.make.A(mesh,locStations)
aLoc = rep((1:dim(locStations)[1]),dtLength)-1
#---------------------------------------------------------

#Plot the triangulation-----------------------------------
plot(mesh)
lines(borders, lwd=3)
points(coordinates$UTMX, coordinates$UTMY, pch=20, cex=2, col=2)
#---------------------------------------------------------

#Define the design matrix for the fixed effects-------
X <- model.matrix( ~ 1 + Piemonte_data$A + Piemonte_data$UTMX + Piemonte_data$UTMY +
                     Piemonte_data$WS + Piemonte_data$TEMP + Piemonte_data$HMIX +
                     Piemonte_data$PREC + Piemonte_data$EMI , data = Piemonte_data)
#-----------------------------------------------------

#Define data and parameter object--------------------
data <- list(logPM10 = Piemonte_data$logPM10,
             n_data = n_data,
             lengthDt = lengthDt,
             maxDt = maxDt,
             A = A,
             aLoc = aLoc,
             X = as.matrix(X),
             spdeMatrices = spdeMatrices
)

parameters <- list(beta      = c(3,rep(0,8)),
                   log_tau   = 4,
                   log_kappa = -4,
                   rhoTan = 0,
                   x  = array(0,dim = c(mesh$n,maxDt)),
                   logSigmaE = -1)

#-----------------------------------------------------

#Fit model with SEPARABLE formulation-----------------
startTime <- Sys.time()
data$flag = 1
obj <- TMB::MakeADFun(data,parameters,random = c("x"),DLL = "SPDExAR1")
obj <- normalize(obj, flag="flag")
opt<-stats::nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=1000, iter.max=1000))
rep<- sdreport(obj)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)
#-----------------------------------------------------


#Extract the range---------------
rangeIndex = which(row.names(summary(rep,"report"))=="range")
range = summary(rep,"report")[rangeIndex,]
#--------------------------------

#Construct .git file with spatio-termporal illustration------
saveGIF({
for(dtPlot in 1:(maxDt-1))
  {
    proj = inla.mesh.projector(mesh)
    gammaindeks = which(names(rep$par.random)=="gamma")
    gammaindeks = which(names(rep$par.random)=="xVec")
    XYUTM = mesh$loc[,1:2]
    indexOfGammaToPlot = ((dtPlot)*mesh$n+1):((dtPlot+1)*mesh$n)
    gamma = rep$par.random[indexOfGammaToPlot]/exp(rep$par.fixed[which(names(rep$par.fixed)=="log_tau")])

    latentFieldML = gamma
    image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldML),col =  colorRampPalette(c("white","yellow", "red"))(12),
               xlab = 'Easting', ylab = 'Northing',
               zlim = c(-2,2),
               main = paste("MAP estimate of the spatio-temporal field day ", (dtPlot-1)*dtLength +1 ," to ",dtPlot*dtLength,sep = ""),
               cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)
    contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldML),nlevels = 6 ,add = T,labcex  = 1,cex = 1)
    lines(borders, lwd=3)
  }
}, movie.name = "SPDExAR1.gif")
#-----------------------------------------------------------
