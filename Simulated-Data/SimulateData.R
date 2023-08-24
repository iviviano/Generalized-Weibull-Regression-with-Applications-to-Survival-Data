library(flexsurv)

###############################################################################
#Censoring Mechanism
###############################################################################
#data is the data to be censored
#frac is the fraction of data to be censored
#returns a list with the censored data and a vector indicating which were censored
#Algorithm: 
#randomly choose the specified fraction of data to be censored
#for each censored measurement, the reported value is chosen uniformly at 
#random from the interval from 0 to the actual end time
censor <- function(data, frac) {
  N <- length(data)
  trun = numeric(N)
  
  for (i in sample(1:N, frac * N, replace = FALSE)) {
    trun[i] <- 1
    
    data[i] <- runif(1, 0, data[i])
  }
  list(data = data, trun = trun)
}

set.seed(1)

beta0 <- 2.5
betas <- c(.3, -.6)
N <- 100
x1 <- rnorm(N, 0, .5); x2 <- rbinom(N, size = 2, prob = .2)
X <- cbind((x1 - mean(x1)) * .5 / sd(x1), x2 - mean(x2))
P <- ncol(X)

set.seed(runif(1))

sigma <- 2; lambda <- 2

DFR <- numeric(N)

for (i in 1:N)
  DFR[i] <- rgengamma.orig(1, shape = lambda, scale = sigma, k = (beta0 - X %*% betas)[i])


data <- censor(DFR, .25)

censoring <- numeric(N)
cens <- data$trun

T <- DFR
t <- data$data



samples <- list()

###DFR: 
#Gen Weib: alpha >= 1, lambda <= 0
#Gen gam: sigma > 1, lambda >= 1 / sigma

sigma <- 2; lambda <- 2

DFR <- numeric(N)

for (i in 1:N)
  DFR[i] <- rgengamma.orig(1, shape = lambda, scale = sigma, k = (beta0 - X %*% betas)[i])

samples$DFR1 <- DFR


sigma <- 2; lambda = .95

for (i in 1:N)
  DFR[i] <- rgengamma.orig(1, shape = lambda, scale = sigma, k = (beta0 - X %*% betas)[i])

samples$DFR2 <- DFR

###Bathtub: 
#Gen Weib: alpha > 1, lambda > 0
#Gen gam: lambda > max(sigma, 1 / sigma)

sigma <- 1; lambda <- 2

bathtub <- numeric(N)

for (i in 1:N)
  bathtub[i] <- rgengamma.orig(1, shape = lambda, scale = sigma, k = (beta0 - X %*% betas)[i])

samples$bath1 <- bathtub


sigma <- 2; lambda <- 3

for (i in 1:N)
  bathtub[i] <- rgengamma.orig(1, shape = lambda, scale = sigma, k = (beta0 - X %*% betas)[i])

samples$bath2 <- bathtub

###Unimodal: 
#Gen Weib: alpha < 1, lambda < 0
#Gen Gam: lambda < min(sigma, 1 / sigma)

sigma <- 1; lambda <- .5

uni <- numeric(N)

for (i in 1:N)
  uni[i] <- rgengamma.orig(1, shape = lambda, scale = sigma, k = (beta0 - X %*% betas)[i])

samples$uni1 <- uni


sigma <- 2; lambda <- .25

for (i in 1:N)
  uni[i] <- rgengamma.orig(1, shape = lambda, scale = sigma, k = (beta0 - X %*% betas)[i])

samples$uni <- uni

###IFR: 
#Gen Weib: alpha <= 1, lambda >= 0
#Gen gam: 0 < sigma < 1, lambda <= 1 / sigma

sigma <- .5; lambda <- .5

IFR <- numeric(N)

for (i in 1:N)
  IFR[i] <- rgengamma.orig(1, shape = lambda, scale = sigma, k = (beta0 - X %*% betas)[i])

samples$IFR1 <- IFR


sigma <- .25; lambda <- 3
for (i in 1:N)
  IFR[i] <- rgengamma.orig(1, shape = lambda, scale = sigma, k = (beta0 - X %*% betas)[i])

samples$IFR2 <- IFR






sim_data <- as.data.frame(X)

names <- c("cont", "cat")
censoring_fracs <- c(.1, .3, .5, .7)

for (i in 1:length(samples)) {
  col <- 20 * samples[[i]]
  for (frac in censoring_fracs) {
    censored <- censor(col, frac = frac)
    sim_data <- cbind(sim_data, censored$data, censored$trun, col)
    names <- append(names, c(sprintf("%s (%d%%)", names(samples)[[i]], 100 * frac), "status", "true"))
  }
}

names(sim_data) <- names

write.csv(sim_data, file = "sim_data.csv")
