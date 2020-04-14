#Currently needs a version of TMB which includes the barrier procedure
#devtools::install_github("OlavNikolaiBreivik/adcomp/TMB", ref = "barrierAndAD")

library(INLA)
library(fields)
library(TMB)

#Compile and load c-code
compile("spdeBarrier.cpp")
dyn.load("spdeBarrier")

#Download and read data
download.file(url = "https://haakonbakka.bitbucket.io/data/WebSiteData-Archipelago.RData", destfile = "WebSiteData-Archipelago.RData")
load(file = "WebSiteData-Archipelago.RData")

#Construtc mesh
mesh = inla.mesh.2d(boundary = poly.water,
                    loc=cbind(df$locx, df$locy),
                    max.edge = c(1,5)*0.6,
                    cutoff = 0.06,
                    offset = c(0.6, 4.6))
#Plot mesh
png(file = "mapAndMesh.png",height = 900,width = 900)
plot(mesh) 
points(df$locx, df$locy, col="red",cex = sqrt(df$y.smelt)+1)
dev.off()

#Construct spatial interpolation matrix
A = inla.spde.make.A(mesh, loc = cbind(df$locx, df$locy))


#Construct barrier.triangles
tl = length(mesh$graph$tv[,1])
posTri = matrix(0, tl, 2)
for(t in 1:tl){
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)]
}
posTri = SpatialPoints(posTri)
normal = unlist(over(poly.water, posTri, returnList=T))
barrier.triangles = setdiff(1:tl, normal)


#Construct spde object and structure needed in the spde-procedure
spde = inla.spde2.matern(mesh, alpha=2)
spdeMatrices = spde$param.inla[c("M0","M1","M2")] #Matrices needed in Standard spde-procedure
fem = INLA:::inla.barrier.fem(mesh = mesh, barrier.triangles = barrier.triangles)
spdeMatricesBarrier = list(C0 = fem$C[[1]],C1 =fem$C[[2]] ,D0 = fem$D[[1]],D1 = fem$D[[2]],I = fem$I )



#Construct design matrix for fixed effects
#X <- model.matrix( ~ 1 + dptLUKE + dptavg15km + dist30m + joetdsumsq + lined15km + swmlog10 + temjul15,
#                   data = df) #Use covariates
X <- model.matrix( ~ 1,
                   data = df) #No covariates



#Define data and parameters-------
data = list(y = df$y.smelt,
            A = A,
            spdeMatrices = spdeMatrices,
            spdeMatricesBarrier = spdeMatricesBarrier,
            barrier = 1,
            c = c(1,0.2),
            X = as.matrix(X)
)

par = list(beta = rep(0,dim(X)[2]),
           log_tau =1,
           log_kappa = -3,
           log_sigmaIID = -3,
           x = rep(0,mesh$n),
           xiid = rep(0,length(df$y.smelt)))
#--------------------------------

#Estimating the model and extract results-------------
startTime <- Sys.time()
obj <- MakeADFun(data, par, random=c("x", "xiid"), DLL="spdeBarrier")
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport(obj)
endTime = Sys.time()
timeUsed = endTime - startTime
print(timeUsed)
#------------------------

#Extract range
rangeIndex = which(row.names(summary(rep,"report"))=="range")
range = summary(rep,"report")[rangeIndex,]

#Plot spatial effect
png("spatialBarrierTMB.png")
proj = inla.mesh.projector(mesh)
latentFieldMAP = rep$par.random[names(rep$par.random)=="x"]/exp(rep$par.fixed[which(names(rep$par.fixed)=="log_tau")])
image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),col =  colorRampPalette(c("white","yellow", "red"))(12),
           xlab = 'Easting', ylab = 'Northing',zlim = c(-4,4),
           main = "MAP estimate of spatial latent field",
           cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)
plot(inla.barrier.polygon(mesh, barrier.triangles), add=T)
contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldMAP) ,add = T,labcex  = 1,cex = 1)
dev.off()








