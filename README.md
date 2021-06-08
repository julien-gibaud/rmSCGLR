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
