
# Load data
load("Robs.RData")

library(TMB)

# Compile TMB code for for random walk
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
summary(rep)


# Compile TMB code for random walk state space model
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


# Plot estimated recruitment process
plot(Robs$year, log(Robs$Robs))
lines(Robs$year, logR_rep[, 1], col="black")
lines(Robs$year, logR_rep[, 1] - 2 * logR_rep[, 2], col="black", lty="dashed")
lines(Robs$year, logR_rep[, 1] + 2 * logR_rep[, 2], col="black", lty="dashed")



# Rickers method

compile("recruitment_state_space_SSB.cpp")
dyn.load(dynlib("recruitment_state_space_SSB"))


data <- list(Robs = Robs$Robs,
             SSB = Robs$ssb,
             method = 1)

param <- list(logR = rep(0, length(Robs$Robs)),
              log_sigma_Robs = 0,
              log_sigma_logR = 0, 
              ricker = if(data$method == 1){numeric(2)} else{numeric(0)},
              bh = if(data$method == 2){numeric(2)} else{numeric(0)})


  
obj <- MakeADFun(data, param, random = "logR", DLL = "recruitment_state_space_SSB")
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
logR_rep <- summary(rep, "random")

lines(Robs$year, logR_rep[, 1], col="red")



  

# Beverton-Holt


data <- list(Robs = Robs$Robs,
             SSB = Robs$ssb,
             method = 2)

param <- list(logR = rep(0, length(Robs$Robs)),
              log_sigma_Robs = 0,
              log_sigma_logR = 0, 
              ricker = if(data$method == 1){numeric(2)} else{numeric(0)},
              bh = if(data$method == 2){numeric(2)} else{numeric(0)})



obj <- MakeADFun(data, param, random = "logR", DLL = "recruitment_state_space_SSB")
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
logR_rep <- summary(rep, "random")

lines(Robs$year, logR_rep[, 1], col = "blue")
legend("bottomright", c("Random walk", "Ricker", "Berverton-Hold"), lty = 1, col = c("black", "red", "blue"))

