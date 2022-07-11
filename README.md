# rmSCGLR
Supervised Component-based Generalised Linear Regression for mixture model on the responses

## Installation

``` r
# Install development version from GitHub
remotes::install_github("julien-gibaud/rmSCGLR")
```
## Example

``` r
library(rmSCGLR)

# load sample data
data <- genus
data <- as.data.frame(apply(data, 2, as.numeric ))

# get variable names from dataset
n <- names(data)
ny <- n[grep("^gen",n)]
nx <- n[-grep("^gen",n)]
na <- c("geology")
nx <- nx[!nx%in%c("geology", "surface")]

# build multivariate formula
form <- multivariateFormula(Y = ny, X = nx, A = na)

# define family
fam <- rep("poisson", length(ny))

# run function
H <- c(2,2)
met <- methodSR.RMSCGLR(l=4, s=0.1, t=0.5)
res <- ResponseMixtureSCGLR(formula=form, data=data, H=H,
                           family=fam, method=met, offset=data$surface)
                           
# print results
res$U
res$comp
res$cluster

# plot the results
plot_RMSCGLR(x=res, thresold=0.5, group=1, plan=c(1,2))
```

## Simulations from the paper ([Gibaud et al. 2022](#ref-gibaud22))
### Simulation 1
```r
library(rmSCGLR)
library(mvtnorm)

#**********#
# settings #
#**********#
N <- 100     # observations
K <- 100     # responses
rho <- 0.9   # correlation between the first and second latent variables
variance_bundle <- 0.1 # variance within the bundles

#*****************************#
# create the latent variables #
#*****************************#
Gamma <- matrix(rho, 2, 2); for(i in 1:2) Gamma[i, i] <- 1
PSI.sim <- rmvnorm(N, mean = rep(0, 2), sigma = Gamma)
psi1.sim <- PSI.sim[,1]   # first latent variable
psi2.sim <- PSI.sim[,2]   # second latent variable
psi3.sim <- rnorm(n = N)  # third latent variable
psi4.sim <- rnorm(n = N)  # fourth latent variable

#**********************************#
# create the regression parameters #
#**********************************#
# regression parameters for latent variables 1 and 2 
gamma.sim <- runif(n = K, min = -4, max = 4)
# regression parameters for latent variables 3 and 4
gamma.sim1 <- runif(n = K, min = -2, max = 2)

#*********************************#
# simulate the response variables #
#*********************************#
Y <- cbind()
# first group
for(i in 1:20){ # gaussian
  eta <- psi1.sim*gamma.sim[i]+psi3.sim*gamma.sim1[i]
  mu <- eta
  Y <- cbind(Y, rnorm(n=N, mean=mu, sd=sqrt(1)))
}
for(i in 21:70){ # poisson
  eta <- psi1.sim*gamma.sim[i]*0.25+psi3.sim*gamma.sim1[i]*0.25
  mu <- exp(eta)
  Y <- cbind(Y, rpois(n=N, lambda=mu))
}
# second group
for(i in 71:80){ # gaussian
  eta <- psi2.sim*gamma.sim[i]+psi4.sim*gamma.sim1[i]
  mu <- eta
  Y <- cbind(Y, rnorm(n=N, mean=mu, sd=sqrt(1)))
}
for(i in 81:100){ # bernoulli
  eta <- psi2.sim*gamma.sim[i]+psi4.sim*gamma.sim1[i]
  mu <- exp(eta)/(1+exp(eta))
  Y <- cbind(Y, rbinom(n=N, size=1, prob=mu))
}

#************************************#
# simulate the explanatory variables #
#************************************#
X <- cbind()
for(i in 1:20)  X <- cbind(X, psi1.sim + rnorm(n = N, sd = sqrt(variance_bundle)))
for(i in 1:20)  X <- cbind(X, psi2.sim + rnorm(n = N, sd = sqrt(variance_bundle)))
for(i in 1:10)  X <- cbind(X, psi3.sim + rnorm(n = N, sd = sqrt(variance_bundle)))
for(i in 1:10)  X <- cbind(X, psi4.sim + rnorm(n = N, sd = sqrt(variance_bundle)))
for(i in 1:40)  X <- cbind(X, rnorm(n = N, sd = 1))

#**************#
# run function #
#**************#
# build data
data <- data.frame(cbind(Y, X))
# build multivariate formula
ny <- paste("y", 1:ncol(Y), sep = "")
nx <- paste("x", 1:ncol(X), sep = "")
colnames(data) <- c(ny, nx)
form <- multivariateFormula(Y = ny, X = nx)
# define family
fam <- c(rep("gaussian", 20), rep("poisson", 50),
         rep("gaussian", 10), rep("bernoulli", 20))
# define method
met <- methodSR.RMSCGLR(l=4, s=0.1, t=0.4)
# run
res.scglrRM <- ResponseMixtureSCGLR(formula=form, data=data, H=c(2,2),
                                    family=fam, method=met)

#**************#
# show results #
#**************#
# the correlation scatterplots
plot1 <- plot_RMSCGLR(x=res.scglrRM, thresold = 0.5, group=1, plan = c(1,2))
plot2 <- plot_RMSCGLR(x=res.scglrRM, thresold = 0.5, group=2, plan = c(1,2))
# the prior probability of belonging to the group
res.scglrRM$pg
# the posterior group membership probabilities
res.scglrRM$postProbs
# the cluster to which each response is allocated
res.scglrRM$cluster
# the supervised components
res.scglrRM$comp
```
### Simulation 2

```r
library(rmSCGLR)
library(mvtnorm)

#**********#
# settings #
#**********#
N <- 100     # observations
K <- 100     # responses
rho <- 0.5   # correlation between the first, third and fifth latent variables
variance_bundle <- 0.1 # variance within the bundles

#*****************************#
# create the latent variables #
#*****************************#
Gamma <- matrix(rho, 3, 3); for(i in 1:3) Gamma[i, i] <- 1
PSI.sim <- rmvnorm(N, mean = rep(0, 3), sigma = Gamma)
psi1.sim <- PSI.sim[,1]  # first latent variable
psi3.sim <- PSI.sim[,2]  # third latent variable
psi5.sim <- PSI.sim[,3]  # fifth latent variable
psi2.sim <- rnorm(n = N) # second latent variable
psi4.sim <- rnorm(n = N) # fourth latent varieble

#**********************************#
# create the regression parameters #
#**********************************#
# regression parameters for latent variables 1, 2 and 3 
gamma.sim <- sample(c(runif(n = K/2, min = 2, max = 4),
                      runif(n = K/2, min = -4, max = -2)))
# regression parameters for latent variables 4 and 5
gamma.sim1 <- sample(c(runif(n = K/2, min = 1, max = 2), 
                       runif(n = K/2, min = -2, max = -1)))

#*********************************#
# simulate the response variables #
#*********************************#
Y <- cbind()
# first group
for(i in 1:20){ # gaussian
  eta <- psi1.sim*gamma.sim[i]+psi4.sim*gamma.sim1[i]
  mu <- eta
  Y <- cbind(Y, rnorm(n=N, mean=mu, sd=sqrt(1)))
}
# second group
for(i in 21:70){ # poisson
  eta <- psi2.sim*gamma.sim[i]*0.25+psi5.sim*gamma.sim1[i]*0.25
  mu <- exp(eta)
  Y <- cbind(Y, rpois(n=N, lambda=mu))
}
# third group
for(i in 71:100){ # bernoulli
  eta <- psi3.sim*gamma.sim[i]
  mu <- exp(eta)/(1+exp(eta))
  Y <- cbind(Y, rbinom(n=N, size=1, prob=mu))
}

#************************************#
# simulate the explanatory variables #
#************************************#
X <- cbind()
for(i in 1:50)  X <- cbind(X, psi1.sim + rnorm(n = N, sd = sqrt(variance_bundle)))
for(i in 1:40)  X <- cbind(X, psi2.sim + rnorm(n = N, sd = sqrt(variance_bundle)))
for(i in 1:30)  X <- cbind(X, psi3.sim + rnorm(n = N, sd = sqrt(variance_bundle)))
for(i in 1:20)  X <- cbind(X, psi4.sim + rnorm(n = N, sd = sqrt(variance_bundle)))
for(i in 1:10)  X <- cbind(X, psi5.sim + rnorm(n = N, sd = sqrt(variance_bundle)))
for(i in 1:50)  X <- cbind(X, rnorm(n = N, sd = 1))

#**************#
# run function #
#**************#
# build data
data <- data.frame(cbind(Y, X))
# build multivariate formula
ny <- paste("y", 1:ncol(Y), sep = "")
nx <- paste("x", 1:ncol(X), sep = "")
colnames(data) <- c(ny, nx)
form <- multivariateFormula(Y = ny, X = nx)
# define family
fam <- c(rep("gaussian", 20), rep("poisson", 50),
         rep("bernoulli", 30))
# define method
met <- methodSR.RMSCGLR(l=4, s=0.1, t=0.6)
# run
res.scglrRM <- ResponseMixtureSCGLR(formula=form, data=data, H=c(2,2,1),
                                    family=fam, method=met)

#**************#
# show results #
#**************#
# the correlation scatterplots
plot1 <- plot_RMSCGLR(x=res.scglrRM, thresold = 0.5, group=1, plan = c(1,2))
plot2 <- plot_RMSCGLR(x=res.scglrRM, thresold = 0.5, group=2, plan = c(1,2))
# the prior probability of belonging to the group
res.scglrRM$pg
# the posterior group membership probabilities
res.scglrRM$postProbs
# the cluster to which each response is allocated
res.scglrRM$cluster
# the supervised components
res.scglrRM$comp
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

If you use this package, please cite

<div id="ref-gibaud22" class="csl-entry"> 

Gibaud J., Bry X., Trottier C., Mortier F. and Réjou-Méchain M., (2022), ``Response mixture models based on supervised components: Clustering floristic taxa'', *Statistical Modelling*

</div>
  
</div>
