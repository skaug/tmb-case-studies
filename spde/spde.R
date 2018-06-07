library(TMB)
library(INLA)
library(fields)

#Compile and load c++ code-------
compile("spde.cpp")
dyn.load(dynlib("spde"))
#--------------------------------

#Read data-----------------------
map = read.table("Leuk.map")
data(Leuk,package = "INLA")
#--------------------------------

#Define mesh and components representing the  precision matrix----
loc = cbind(Leuk$xcoord, Leuk$ycoord)
boundary = INLA::inla.nonconvex.hull(loc)
boundary2 = INLA::inla.nonconvex.hull(loc,convex = -0.35)
mesh = INLA::inla.mesh.2d(
  loc=loc,
  boundary = list(boundary,boundary2),
  max.edge=c(0.05, 0.2),
  cutoff=0.05
)
A = inla.spde.make.A(mesh,loc)
spde = inla.spde2.matern(mesh, alpha=2)
spdeMatrices = spde$param.inla[c("M0","M1","M2")]
#-------------------------------------------------------------------

#Define the design matrix for the fixed effects-------
X <- model.matrix( ~ 1 + sex + age + wbc + tpi, data = Leuk)
#-----------------------------------------------------

#Define the data and parameters given to TMB----------
data <- list(time       = Leuk$time,
             notcens    = Leuk$cens,
             meshidxloc = mesh$idx$loc - 1,
             A = A,
             X          = as.matrix(X),
             spdeMatrices = spdeMatrices
             )

parameters <- list(beta      = c(0.0,0,0,0,0),
                   log_tau   = 0,
                   log_kappa = 0,
                   log_omega = -1,
                   x         = rep(0.0, nrow(data$spdeMatrices$M0)) )
#-----------------------------------------------------

#Estimating the model and extract results-------------
#data$flag = 1
startTime <- Sys.time()
obj <- MakeADFun(data, parameters, random="x", DLL="spde")
#obj <- normalize(obj, flag="flag")
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)
#-----------------------------------------------------


rangeIndex = which(row.names(summary(rep,"report"))=="range")
fieldIndex = which(row.names(summary(rep,"report"))=="x")
range = summary(rep,"report")[rangeIndex,]

proj = inla.mesh.projector(mesh)
latentFieldMAP = rep$par.random[names(rep$par.random)=="x"]/exp(rep$par.fixed[which(names(rep$par.fixed)=="log_tau")])
image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),col =  colorRampPalette(c("white","yellow", "red"))(12),
           xlab = 'Easting', ylab = 'Northing',
           main = "MAP estimate of the spatial latent field",
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)
contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldMAP) ,add = T,labcex  = 1,cex = 1)
points(map,type = 'l',lwd = 2)

latentFieldSD = sqrt(rep$diag.cov.random)/exp(rep$par.fixed[which(names(rep$par.fixed)=="log_tau")])
image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldSD),col =  colorRampPalette(c("white","yellow", "red"))(12),
           xlab = 'Easting', ylab = 'Northing',
           main = "Standard deviation of the estimated spatial latent field",
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)
contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldSD) ,add = T,labcex  = 1,cex = 1)
points(map,type = 'l',lwd = 2)

library(SparseM)
rep <- sdreport(obj,getJointPrecision=T)
nameIndex = which(colnames(rep$jointPrecision)=="x")
Q = rep$jointPrecision[nameIndex,nameIndex]
Q[Q!=0] = 1
image(Q, main = "Sparsness structure of the precision matrix of the GMRF")

