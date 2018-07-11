## Script to demonstrate testing the Laplace approximation of a TMB model
## using MCMC. Verisons are run with and without the LA turned on for the
## random effects, and differences in posteriors for fixed effects are
## compared.

## Working directory needs to be set to the folder of this file before
## proceeding.

## Note: these models include priors, and use exponentiation in the
## template with a Jacobian adjustment instead of external bounds of
## (0,Inf) for hypervariance parameters. Thus, the SD parameters are really
## in log space. See model files for more details.

library(TMB)
library(tmbstan)
cores <- parallel::detectCores()-1
options(mc.cores = cores) ## parallel TMB runs with tmbstan

##### Setup inputs to model
set.seed(321)
data <- readRDS('bayesian.RDS')
## Use random inits for each core
inits <- lapply(1:cores, function(x)
  list(yearInterceptSD = runif(1, .05, .15),
       plantInterceptSD = runif(1, .05, .15),
       plantSlopeSD = runif(1, .05, .15),
       intercept = rnorm(data$Nstage, 0, 1),
       slope = rnorm(1, 0, 1),
       yearInterceptEffect_raw= rnorm(data$Nyear, 0, 1),
       plantInterceptEffect_raw= rnorm(data$Nplant, 0, 1),
       plantSlopeEffect_raw= rnorm(data$Nplant, 0, 1)))

## Build and run MLE version of model (but has priors still)
compile('bayesian.cpp')
dyn.load('bayesian')
obj <- MakeADFun(data=data, parameters=inits[[1]],
            random=c('yearInterceptEffect_raw',
                     'plantInterceptEffect_raw',
                     'plantSlopeEffect_raw'), DLL='bayesian')

## Run with thinning since these models mix at different rates (the full
## model mixes slower). We want to compare samples between chains with
## roughly equal ESS otherwise the results are skewed.
iter <- 1000
wm <- 1000
## ___!!! Warning, this make take an hour or two to run !!!!___

## LA turned off (i.e., full MCMC integration)
mcmc <- tmbstan(obj, iter=iter*5+wm, warmup=wm, thin=5,
                init=inits, chains=cores, seed=1,
                control=list(adapt_delta=0.8))
## LA turned on
mcmc.la <- tmbstan(obj, iter=iter+wm, warmup=wm, thin=1,
                   chains=cores, laplace=TRUE, seed=1,
                   init=inits, control=list(adapt_delta=0.8))

## If you wanted to compare efficiency you would do minESS/time between the
## different versions. Here we just look at potential issues with the LA.
x1 <- as.matrix(mcmc)
x2 <- as.matrix(mcmc.la)
stopifnot(nrow(x1)==nrow(x2))
## Get the fixed effects
post <- data.frame(rbind(x1[,dimnames(x2)[[2]]], x2[,dimnames(x2)[[2]]]))
post$model <- as.factor(rep(c("Full MCMC", "LA"), each=nrow(x1)))
## Randomize order to prevent overplotting
set.seed(2342)
ind <- sample(1:nrow(post), size=nrow(post))
post <- post[ind,]
png("pairs_bayesian_LA.png", width=7, height=5, units='in', res=500)
pairs(post[,1:5], col=post$model, pch=16, cex=.5)
dev.off()
## library(ggplot2)
## ggplot(post, aes(plantSlopeSD)) + geom_histogram() + facet_wrap('model')
