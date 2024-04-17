# Testing

library(tidyverse)
library(tweedie)
library(tweedieCARpanel)
library(DHARMa)

### Simulating some data
S <- 100

# AR1 Neighbors matrix
B <- matrix(data = 0, nrow = S, ncol = S)
for(i in 1:S){
  for(j in 1:S){
    if(abs(i - j) == 1){
      B[i, j] <- 1
    }
  }
}

# row-standardization for simulation
D <- diag(rowSums(B))
rho <- 0.7
sigma2 <- 0.3

Sigma = sigma2 * solve(D - rho * B)

s <- mvtnorm::rmvnorm(1, sigma = Sigma) %>% as.double()

# Construct data
dat <- tibble(
  grp = rep(1:S, each = 5),
  x = rnorm(5 * S, sd = 2),
  y = rtweedie(5 * S, xi = 1.5, mu = exp(-0.1 + .5 * x + s[grp]), phi = 1)
)

# fit the model and simulate for checking with DHARMa

fit <- fit_panel_CAR_tw(
  form = y ~ x,
  speff = "grp",
  B = B,
  data = dat, sim = T
)


# now use simulations with DHARMa to check residuals
qresids <- createDHARMa(
  simulatedResponse = fit$y_rep,
  observedResponse = dat$y,
  integerResponse = T
)
plot(qresids)




