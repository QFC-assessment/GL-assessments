# ---------------------------------------------------------
## Equilibrium Fmsy calculations 
# Walters and Martell 2004 - box 3.1 equations
# conversion between finite and instantaneous fishing mortality
# U = 1 - exp(-F)
# F = -log(1 - U)
# ---------------------------------------------------------

rm(list = ls())

# ---------------------------------------------------------
## Leading Parameters ####
amax <- 20 # maximum age
ages <- seq(1, amax) # age vector
nages <- length(ages) # number of ages
M <- 0.2 # natural mortality

# length at age - von Bertalanffy growth
linf <- 630 # von Bertalanffy growth parameter - asymptotic length
vbk <- 0.3 # von Bertalanffy growth parameter - growth rate
to <- -0.01 # age at time = 0
la <- linf * (1 - exp(-vbk * (ages - to)))

# calculate proportions mature
mat50 <- 378.6 # length or age
sigma_mat <- 0.215
mat <- 1 / (1 + exp(-sigma_mat * (la - mat50))) # for length
# mat <- 1 / (1 + exp(-sig_mat * (ages - mat50))) # for ages

# length-weight regression
alw <- 4.17e-6
alwb <- 3.11
wa <- alw * la^alwb

# vulnerability at age
vul50 <- exp(6.178401) # length or age
sigma_vul <- exp(-3.190373)
vul <- 1 / (1 + exp(-sigma_vul * (la - vul50))) # for length
# vul <- 1 / (1 + exp(-sigma_vul * (ages - vul50))) # for ages


# ---------------------------------------------------------
## At age calculations ####
# survivorship in unfished conditions
S <- exp(-M)
lo <- numeric(nages)
lo[1] <- 1
for (a in 2:nages) {
  lo[a] <- lo[a - 1] * S
}
lo[nages] <- lo[nages] / (1 - S) # plus group
# spawner biomass per recruit in unfished condition
sbro <- sum(lo * wa * mat)

# ---------------------------------------------------------
## Recruitment parameters ####
Ro <- 1 # recruitment scaler
cr <- 7 # compensation ratio
# alpha parameter for stock-recruitment relationship (in log space) - same for Ricker or Beverton-Holt
ln_ar <- log(cr / sbro)

# Beverton-Holt beta parameter
bh <- (exp(ln_ar) * sbro - 1) / (Ro * sbro)

# plot stock-recruitment curves
ssb <- 0:(sbro) * 1.1
bh_recruits <- ssb * ((exp(ln_ar) / (1 + bh * ssb)))
plot(ssb, bh_recruits) # does it look correct?


# ---------------------------------------------------------
## Fmsy calculations - for loop ####
# per recruit calculations
F_vec <- seq(0, 2, length.out = 1000)
Req <- numeric(length(F_vec))
Yeq <- numeric(length(F_vec))
Yield <- 0
Fmsy <- 0

# loop through low to high fishing mortality values
for (i in 1:length(F_vec)) {
  YPR <- sbrf <- lf <- 0
  lf[1] <- Ro # equivalent to Ro = 1
  for (a in 1:nages) {
    # spawner biomass per recruit in fished condition
    sbrf <- sbrf + lf * mat[a] * wa[a]
    # yield per recruit
    YPR <- YPR + lf * wa[a] * vul[a]
    # survivorship in the fished condition
    # lf <- lf * exp(-M) * (1 - U_vec[i] * vul[a])
    lf <- lf * exp(-M) * (1 - (1 - exp(-F_vec[i])) * vul[a])
  }
  # equilibrium recruitment
  Req[i] <- (exp(ln_ar) * sbrf - 1.0) / (bh * sbrf)
  # equilibrium yield
  # Yeq[i] <- YPR * Req[i] * U_vec[i]
  Yeq[i] <- YPR * Req[i] * (1 - exp(-F_vec[i]))

  # find F where yield is maximized
  if (Yeq[i] > Yield) {
    Yield <- Yeq[i]
    Fmsy <- F_vec[i]
  }
  if (Yeq[i] < 0) {
    Yeq[i] <- 0
  }
}

# plot fishing mortality vs yield
plot(F_vec, Yeq)
abline(v = Fmsy, lty = 2)

# Fmsy
print(Fmsy)


# ---------------------------------------------------------
## Fmsy calculations - uniroot ####
# use the "tape" function from RTMB, which finds the derivative (jacobian - first-order partial)
library(RTMB)
# "data" list
data <- list()
data$amax <- 20 # maximum age
data$ages <- seq(1, amax) # age vector
data$nages <- length(ages) # number of ages
data$M <- 0.2 # natural mortality

# length at age - von Bertalanffy growth
linf <- 630 # von Bertalanffy growth parameter - asymptotic length
vbk <- 0.3 # von Bertalanffy growth parameter - growth rate
to <- -0.01 # age at time = 0
data$la <- linf * (1 - exp(-vbk * (ages - to)))

# calculate proportions mature
mat50 <- 378.6 # length or age
sigma_mat <- 0.215
data$mat <- 1 / (1 + exp(-sigma_mat * (la - mat50))) # for length
# mat <- 1 / (1 + exp(-sig_mat * (ages - mat50))) # for ages

# length-weight regression
alw <- 4.17e-6
alwb <- 3.11
data$wa <- alw * la^alwb

# vulnerability at age
vul50 <- exp(6.178401) # length or age
sigma_vul <- exp(-3.190373)
data$vul <- 1 / (1 + exp(-sigma_vul * (la - vul50))) # for length
# vul <- 1 / (1 + exp(-sigma_vul * (ages - vul50))) # for ages

get_eq_yield <- function(F, cr) {
  getall(data)

  ## Recruitment parameters
  Ro <- 1 # recruitment scaler
  # alpha parameter for stock-recruitment relationship (in log space) - same for Ricker or Beverton-Holt
  ln_ar <- log(cr / sbro)

  # Beverton-Holt beta parameter
  bh <- (exp(ln_ar) * sbro - 1) / (Ro * sbro)

  ## Fmsy calculations
  Req <- Yeq <- Yield <- Fmsy <- YPR <- sbrf <- lf <- 0
  lf[1] <- Ro # equivalent to Ro = 1
  for (a in 1:nages) {
    # spawner biomass per recruit in fished condition
    sbrf <- sbrf + lf * mat[a] * wa[a]
    # yield per recruit
    YPR <- YPR + lf * wa[a] * vul[a]
    # survivorship in the fished condition
    lf <- lf * exp(-M) * (1 - (1 - exp(-F)) * vul[a])
  }
  # equilibrium recruitment
  Req <- (exp(ln_ar) * sbrf - 1.0) / (bh * sbrf)
  # equilibrium yield
  Yeq <- YPR * Req * (1 - exp(-F))

  return(Yeq)
}

# Create a tape for the function once
# find Fmsy based on a known compensation ratio
cr <- 7
tape <- MakeTape(function(F) get_fmsy(F, cr), numeric(1))
jacfun <- tape$jacfun()
# uniroot - find Fmsy
root <- uniroot(function(x) jacfun(x), c(1e-3, 2))$root
# Fmsy
print(root)

# compare Fmsy values
print(Fmsy)
print(root)
