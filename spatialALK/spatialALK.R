library(TMB)
library(INLA) #Use SPDE functionality
library(rgdal) #use readOGR
library(mgcv) #Use for spline functionality
library(maptools) #use pointLabel()
library(fields) #use image.plot()
library(rworldxtra) #Use dataset for ploting

compile("spatialALK.cpp")
dyn.load(dynlib("spatialALK"))

#Read data---------------------------------------------
load("cod2015Q1.RData")
plusGroup = 6;
ca_hh$Age[ca_hh$Age>plusGroup] = plusGroup

polygons <- readOGR("Roundfish_shapefiles")
data(countriesHigh)
map <- countriesHigh
#------------------------------------------------------


#Extrtact stuff needed for splines---------------------
#Use functionality in mgcv
gam_setup = gam(Age ~ s(LngtCm, bs = "cs"),
                data = ca_hh,fit=FALSE)

#Penelization matrix
S = gam_setup$smooth[[1]]$S[[1]]
S_list = list(rep(list(S),plusGroup-1))[[1]]
S_combined = .bdiag(S_list)         # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S

#Design matrix
X_list = list(rep(list(gam_setup$X[,-1]),plusGroup-1))[[1]]

#For report
LngtCm=seq(min(ca_hh$LngtCm),100,by = 1)
lengthReport = PredictMat(gam_setup$smooth[[1]],data = data.frame(LngtCm))
lengthReport_list = list(rep(list(lengthReport),plusGroup-1))[[1]]
#----------------------------------------------------------------------------


#Define the SPDE-------------------------------------------------------------
#Extract the locations and convert to UTM-coordinates
n = dim(ca_hh)[1]
XY = data.frame(ca_hh$lon, ca_hh$lat)
names(XY) = c("X","Y")

UTMcenter = 31
latLon = SpatialPoints(XY, proj4string = CRS(as.character(NA)))
proj4string(latLon) ="+proj=longlat"
XYUTM= spTransform(latLon,paste("+proj=utm +",UTMcenter,"V +units=km",sep = ""))

cutoff = 50
maxEdge = c(100,140)
points = matrix(0,n,2)
points[,1] = XYUTM$X
points[,2] = XYUTM$Y
boundary = INLA::inla.nonconvex.hull(points)
boundary2 = INLA::inla.nonconvex.hull(points,convex = -0.45)
mesh = INLA::inla.mesh.2d(loc = points, boundary = list(boundary,boundary2),max.edge = maxEdge,cutoff = cutoff)
A = inla.spde.make.A(mesh,points)
spde = INLA::inla.spde2.matern(mesh = mesh , alpha = 2)#Defines the parts needed in the precision matrix
#----------------------------------------------------------------------------


#Define parameters--------------------------------
data = list(
  age = ca_hh$Age,
  antObs = length(ca_hh$Age),
  length = ca_hh$LngtCm,
  Sdims = Sdims,
  S = S_combined,
  X = .bdiag(X_list), # Design matrix, without intercept
  designMatrixForReport = .bdiag(lengthReport_list),
  A =  A,
  spde = spde$param.inla[c("M0","M1","M2")],
  plusGroup = plusGroup
)


parameters <- list(
  beta0 = rep(0,plusGroup-1),
  betaLength = rep(0,dim(S)[1]*(plusGroup-1)),
  log_lambda = rep(0,plusGroup-1),
  x1 = rep(0,mesh$n),
  x2 = rep(0,mesh$n),
  x3 = rep(0,mesh$n),
  x4 = rep(0,mesh$n),
  x5 = rep(0,mesh$n),
  logKappa = -4,
  logTau = 4
)
#-----------------------------------------------------------------------

#Estimating the model and extract results-------------
timeBefore=Sys.time();
data$flag = 1
obj <- TMB::MakeADFun(data,parameters,
                      random = c("x1","x2","x3","x4","x5", "betaLength"),
                      DLL = "spatialALK")
obj <- normalize(obj, flag="flag")
opt<-stats::nlminb(obj$par,obj$fn,obj$gr,
                   control = list(iter.max = 300,
                                  eval.max = 300))
rep<-TMB::sdreport(obj, getJointPrecision=TRUE)
timeUsed = difftime(Sys.time(),timeBefore,units="mins")
print(paste("Time used to fit model and report: ", round(timeUsed,2),".",sep=""))
#---------------------------------------------------


#Plot length contribution to the linear predictor---------
index = which(names(rep$value)=="repLength")
muSpline1 = exp(rep$par.fixed[1] + rep$value[index[1:dim(lengthReport)[1]]])
muSpline2 = exp(rep$par.fixed[2] + rep$value[index[(dim(lengthReport)[1]+1):(2*dim(lengthReport)[1])]])
muSpline3 = exp(rep$par.fixed[3] + rep$value[index[(2*dim(lengthReport)[1]+1):(3*dim(lengthReport)[1])]])
muSpline4 = exp(rep$par.fixed[4] + rep$value[index[(3*dim(lengthReport)[1]+1):(4*dim(lengthReport)[1])]])
muSpline5 = exp(rep$par.fixed[5] + rep$value[index[(4*dim(lengthReport)[1]+1):(5*dim(lengthReport)[1])]])
sum = muSpline1 + muSpline2 + muSpline3 + muSpline4 + muSpline5

prob1 = muSpline1/(1 + sum)
prob2 = muSpline2/(1 + sum)
prob3 = muSpline3/(1 + sum)
prob4 = muSpline4/(1 + sum)
prob5 = muSpline5/(1 + sum)
prob6 = 1/(1+sum)


plot(LngtCm,prob1 ,
     ylim = c(0,1),
     xlab = "Length", ylab = "Prob", main = "P(age| length) in 2015 Q1 for cod",
     lwd=3, col="red",type = 'l',
     cex.lab=1.5, cex.main = 1.8,cex.axis = 1.2)
lines(LngtCm,prob2 ,
      ylim = c(0,1),
      lwd=3, col="blue",type = 'l')
lines(LngtCm,prob3 ,
      ylim = c(0,1),
      lwd=3, col="black",type = 'l')
lines(LngtCm,prob4 ,
      ylim = c(0,1),
      lwd=3, col="green",type = 'l')
lines(LngtCm,prob5 ,
      ylim = c(0,1),
      lwd=3, col="yellow",type = 'l')
lines(LngtCm,prob6 ,
      ylim = c(0,1),
      lwd=3, col="brown",type = 'l')
abline(h=0)
legend(x=75,y=0.99,legend = c("1 year","2 year","3 year","4 year","5 year",">5 year")
   ,col=c("red","blue","black","green","yellow","brown"),lty = 1,lwd = 3)
#---------------------------------------------------


#Plot spatial and length contribution to linear predictor------------

#Convert the mesh to xy-coordinates for plotting-----------------------
tmp = data.frame(mesh$loc[,1:2])
names(tmp) = c("X","Y")
tmp = SpatialPoints(tmp, proj4string = CRS(as.character(NA)))
proj4string(tmp) =paste("+proj=utm +",UTMcenter,"U +units=km",sep = "")
xyTmp =  spTransform(tmp,"+proj=longlat")
meshXY = mesh
meshXY$loc[,1] = xyTmp$X
meshXY$loc[,2] = xyTmp$Y
projXY = inla.mesh.projector(meshXY)
#-----------------------------------------------------------------------

x1 = which(names(rep$par.random)=="x1")
x2 = which(names(rep$par.random)=="x2")
x3 = which(names(rep$par.random)=="x3")
x4 = which(names(rep$par.random)=="x4")
x5 = which(names(rep$par.random)=="x5")

proj = inla.mesh.projector(mesh)
field1 = rep$par.random[x1]/exp(rep$par.fixed[which(names(rep$par.fixed)=="logTau")])
field2 = rep$par.random[x2]/exp(rep$par.fixed[which(names(rep$par.fixed)=="logTau")])
field3 = rep$par.random[x3]/exp(rep$par.fixed[which(names(rep$par.fixed)=="logTau")])
field4 = rep$par.random[x4]/exp(rep$par.fixed[which(names(rep$par.fixed)=="logTau")])
field5 = rep$par.random[x5]/exp(rep$par.fixed[which(names(rep$par.fixed)=="logTau")])

length = 55
mu1 = exp(rep$par.fixed[1]+ rep$value[index[length-min(ca_hh$LngtCm)+1]] +field1)
mu2 = exp(rep$par.fixed[2]+ rep$value[index[(dim(lengthReport)[1])+length-min(ca_hh$LngtCm)+1]] +field2)
mu3 = exp(rep$par.fixed[3]+ rep$value[index[(2*dim(lengthReport)[1])+length-min(ca_hh$LngtCm)+1]] +field3)
mu4 = exp(rep$par.fixed[4]+ rep$value[index[(3*dim(lengthReport)[1])+length-min(ca_hh$LngtCm)+1]]+ field4)
mu5 = exp(rep$par.fixed[5]+ rep$value[index[(4*dim(lengthReport)[1])+length-min(ca_hh$LngtCm)+1]]+field5)
muSum = mu1 + mu2 + mu3 + mu4 + mu5


probField1 = mu1/(1 + muSum)
probField2 = mu2/(1 + muSum)
probField3 = mu3/(1 + muSum)
probField4 = mu4/(1 + muSum)
probField5 = mu5/(1 + muSum)
probField6 = 1/(1+muSum)

probField1[probField1<0.01] = 0
probField2[probField2<0.01] = 0
probField3[probField3<0.01] = 0
probField4[probField4<0.01] = 0
probField5[probField5<0.01] = 0
probField6[probField6<0.01] = 0


xlim = c(-5,16)
ylim = c(49,62.5)

breaks   = c(0.05,0.1,0.15,0.20,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.999999999999)

xVidde = 5.5;yVidde = 4.5;mainVidde = 3.5
par(mfrow = c(2,3),mar=c(xVidde,yVidde,mainVidde,3.2))
image.plot(projXY$x,projXY$y, inla.mesh.project(projXY, probField1),col =  colorRampPalette(c("white","yellow", "red"))(12),
           main = paste("Prob of age 1 given ", length, " cm",sep = ""), xlab = 'Degrees east', ylab = 'Degrees north',
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,
           xlim = xlim,
           ylim = ylim,
           breaks   =breaks)

contour(projXY$x, projXY$y,inla.mesh.project(projXY, probField1) ,add = T,labcex  = 1,cex = 1,
        breaks   = breaks,
        zlim = c(0.005,0.9))
plot(map, col="grey", add=T)
if (!is.null(polygons)){
  plot(polygons, add=T,lwd = 3)
}
pointLabel(coordinates(polygons),labels=polygons$AreaName)


image.plot(projXY$x, projXY$y, inla.mesh.project(projXY, probField2),col =  colorRampPalette(c("white","yellow", "red"))(12),
           main = paste("Prob of age 2 given ", length, " cm",sep = ""), xlab = 'Degrees east', ylab = 'Degrees north',
           #     main = paste("Probability of age 2 given ", length-1, " cm",sep = ""), xlab = 'Degrees east', ylab = 'Degrees north',
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,
           xlim = xlim,
           ylim = ylim,
           breaks   = breaks)
contour(projXY$x, projXY$y,inla.mesh.project(projXY, probField2) ,add = T,labcex  = 1,cex = 1,
        levels = breaks)
plot(map, col="grey", add=T)
if (!is.null(polygons)){
  plot(polygons, add=T,lwd = 3)
}
pointLabel(coordinates(polygons),labels=polygons$AreaName)


image.plot(projXY$x, projXY$y, inla.mesh.project(projXY, probField3),col =  colorRampPalette(c("white","yellow", "red"))(12),
           main = paste("Prob of age 3 given ", length, " cm",sep = ""), xlab = 'Degrees east', ylab = 'Degrees north',
           #         main = paste("Probability of age 3 given ", length-1, " cm",sep = ""), xlab = 'Degrees east', ylab = 'Degrees north',
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,
           xlim = xlim,
           ylim = ylim,
           breaks   = breaks)
contour(projXY$x, projXY$y,inla.mesh.project(projXY, probField3) ,add = T,labcex  = 1,cex = 1,
        levels = breaks)
plot(map, col="grey", add=T)
if (!is.null(polygons)){
  plot(polygons, add=T,lwd = 3)
}
pointLabel(coordinates(polygons),labels=polygons$AreaName)


image.plot(projXY$x, projXY$y, inla.mesh.project(projXY, probField4),col =  colorRampPalette(c("white","yellow", "red"))(12),
           main = paste("Prob of age 4 given ", length, " cm",sep = ""), xlab = 'Degrees east', ylab = 'Degrees north',
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,
           xlim = xlim,
           ylim = ylim,
           breaks   = breaks)
contour(projXY$x, projXY$y,inla.mesh.project(projXY, probField4) ,add = T,labcex  = 1,cex = 1,
        levels = breaks)
plot(map, col="grey", add=T)
if (!is.null(polygons)){
  plot(polygons, add=T,lwd = 3)
}
pointLabel(coordinates(polygons),labels=polygons$AreaName)



image.plot(projXY$x, projXY$y, inla.mesh.project(projXY, probField5),col =  colorRampPalette(c("white","yellow", "red"))(12),
           main = paste("Prob of age 5 given ", length, " cm",sep = ""), xlab = 'Degrees east', ylab = 'Degrees north',
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,
           xlim = xlim,
           ylim = ylim,
           breaks   = breaks)
contour(projXY$x, projXY$y,inla.mesh.project(projXY, probField5) ,add = T,labcex  = 1,cex = 1,
        levels = breaks)
plot(map, col="grey", add=T)
if (!is.null(polygons)){
  plot(polygons, add=T,lwd = 3)
}
pointLabel(coordinates(polygons),labels=polygons$AreaName)



image.plot(projXY$x, projXY$y, inla.mesh.project(projXY, probField6),col =  colorRampPalette(c("white","yellow", "red"))(12),
           main = paste("Prob of age >5 given ", length, " cm",sep = ""), xlab = 'Degrees east', ylab = 'Degrees north',
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,
           xlim = xlim,
           ylim = c(49.1,62.5),
           breaks   = breaks)
contour(projXY$x, projXY$y,inla.mesh.project(projXY, probField6) ,add = T,labcex  = 1,cex = 1,
        levels = breaks)
plot(map, col="grey", add=T)
if (!is.null(polygons)){
  plot(polygons, add=T,lwd = 3)
}
pointLabel(coordinates(polygons),labels=polygons$AreaName)
#---------------------------------------------------

