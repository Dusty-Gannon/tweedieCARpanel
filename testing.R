# Testing

library(TMB)
library(glmmTMB)
library(tidyverse)
library(tweedie)
devtools::load_all()

### Simulating some data
S <- 100
# Neighbors matrix
B <- matrix(data = 0, nrow = S, ncol = S)
for(i in 1:S){
  for(j in 1:S){
    if(abs(i - j) == 1){
      B[i, j] <- 1
    }
  }
}

# row-standardization
D <- diag(rowSums(B))
rho <- 0.7
sigma2 <- 0.3

Sigma = sigma2 * solve(D - rho * B)

s <- mvtnorm::rmvnorm(1, sigma = Sigma) %>% as.double()

# Construct random effects model
dat <- tibble(
  grp = rep(1:S, each = 5),
  x = rnorm(5 * S, sd = 2),
  y = rtweedie(5 * S, xi = 1.5, mu = exp(-0.1 + .5 * x + s[grp]), phi = 1)
)

# see if we can fit this model

test <- fit_panel_CAR_tw(
  form = y ~ x,
  speff = "grp",
  B = B,
  data = dat
)











### Testing with Kamana's data

load("TMB/testing.RDS")

# tweedie_model <- glmmTMB(
#   update(mod_form, . ~ . + (1 | Vendor.County_f)),
#   family = tweedie,
#   data = data_analysis
# )
#
# # get random effects values
# u <- ranef(tweedie_model)$cond$Vendor.County_f[,1]
# lsig_s <- tweedie_model$fit$par["theta"]



datlist <- list(
  y = data_analysis$county_LIF_M,
  X = model.matrix(mod_form, data = data_analysis),
  group = data_analysis$Vendor.County_f,
  B =
)

parlist <- list(
  beta = rep(0, ncol(datlist$X)),
  log_phi = 0,
  logit_xi = 0,
  log_sigma_s = 0,
  s = rep(0, length(u))
)
obj <- MakeADFun(datlist, parlist, DLL = "tweedieCAR", random = "s")
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
summary(rep)



####
library(TMB)
compile("TMB/linreg.cpp")
dyn.load(dynlib("TMB/linreg"))
set.seed(123)
data <- list(Y = rnorm(10) + 1:10, x=1:10)
parameters <- list(a=0, b=0, logSigma=0)
obj <- MakeADFun(data, parameters, DLL="linreg")
opt <- do.call("optim", obj)
test <- sdreport(obj)



