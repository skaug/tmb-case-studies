load("Robs.RData")

library(TMB)
library(tidyverse)
compile("Robs.cpp")
dyn.load(dynlib("Robs"))

plot(Robs$year, log(Robs$Robs), type = "l")


compile("recruitment.cpp")
dyn.load(dynlib("recruitment"))

# Prepare for TMB
param <- list(log_sigma_r = -0.5)
data <- list(Robs = Robs$Robs)


# Make objective function
obj <- MakeADFun(data, param, DLL = "recruitment")
opt <- nlminb(obj$par, obj$fn, obj$gr)

# Get standard deviations
rep <- sdreport(obj)

# Value of our transformed parameter
rep$value


compile("recruitment_state_space.cpp")
dyn.load(dynlib("recruitment_state_space"))

param <- list(logR = rep(0, length(Robs$Robs)),
              log_sigma_Robs = -0.5,
              log_sigma_logR = -0.5)

data <- list(Robs = Robs$Robs)

obj <- MakeADFun(data, param, random = "logR", DLL = "recruitment_state_space")
opt <- nlminb(obj$par, obj$fn, obj$gr)

# Get standard deviation of parameters and latent process
rep <- sdreport(obj)

# Get summary of fixed effects 
summary(rep, "fixed", p.value = T)
logR_rep <- summary(rep, "random")

# Plot observations and estimated latent process

dat <- tibble(year = Robs$year, log_Robs = log(Robs$Robs), logR = logR_rep[, 1], logR_sd = logR_rep[, 2])

ggplot(dat) + geom_point(aes(year, log_Robs)) + 
  geom_line(aes(year, logR)) + 
  geom_ribbon(aes(year, ymin = log_Robs - 2 * logR_sd, ymax = log_Robs + 2 * logR_sd), fill = "#CC6666", alpha = 0.25)
  ggtitle("Estimated Fish recruitment ")

  
  
# RW

Robs$mode <- 0
par <- list()
par$logsdo <- 0
par$logsdp <- 0
par$logR <- rep(0,length(Robs$Robs))
par$rickerpar <- if(Robs$mode==1){numeric(2)}else{numeric(0)}
par$bhpar <- if(Robs$mode==2){numeric(2)}else{numeric(0)}

obj <- MakeADFun(Robs, par, random="logR", DLL="Robs")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
logR<-as.list(sdr, "Est")$logR
sdlogR<-as.list(sdr, "Std")$logR
lines(Robs$year, logR, col="black")
lines(Robs$year, logR-2*sdlogR, col="black", lty="dashed")
lines(Robs$year, logR+2*sdlogR, col="black", lty="dashed")

# Ricker

Robs$mode <- 1
par <- list()
par$logsdo <- 0
par$logsdp <- 0
par$logR <- rep(0,length(Robs$Robs))
par$rickerpar <- if(Robs$mode==1){numeric(2)}else{numeric(0)}
par$bhpar <- if(Robs$mode==2){numeric(2)}else{numeric(0)}

obj <- MakeADFun(Robs, par, random="logR", DLL="Robs")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
logR<-as.list(sdr, "Est")$logR
lines(Robs$year, logR, col="red")

# Beverton-Holt

Robs$mode <- 2
par <- list()
par$logsdo <- 0
par$logsdp <- 0
par$logR <- rep(0,length(Robs$Robs))
par$rickerpar <- if(Robs$mode==1){numeric(2)}else{numeric(0)}
par$bhpar <- if(Robs$mode==2){numeric(2)}else{numeric(0)}

obj <- MakeADFun(Robs, par, random="logR", DLL="Robs")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
logR<-as.list(sdr, "Est")$logR
lines(Robs$year, logR, col="blue")
