#' @title Function that fits the response mixture scglr model
#'
#' @import ClustOfVar
#' @import pls
#' @importFrom stats rnorm dnorm dpois dbinom
#' @param formula an object of class \code{MultivariateFormula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data a data frame to be modeled.
#' @param H a vector of interger representing the number of components per group
#' @param family a vector of character of the same length as the number of dependent variables:
#' "bernoulli", "binomial", "poisson" or "gaussian" is allowed.
#' @param size describes the number of trials for the binomial dependent variables.
#' @param offset used for the poisson dependent variables.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param crit a list of two elements : maxit and tol, describing respectively the maximum number of iterations and
#' the tolerance convergence criterion for the Expectation Maximization algorithm. Default is set to 50 and 10e-6 respectively.
#' @param method structural relevance criterion. Object of class "method.RMSCGLR"
#'
#' @return \item{U}{the set of loading vectors}
#' @return \item{comp}{the set of components}
#' @return \item{eta}{the fitted linear predictor for each group}
#' @return \item{coef}{contains regression parameters : intercept, gamma and delta}
#' @return \item{pg}{the prior probability of belonging to the group}
#' @return \item{cluster}{a vector of integers indicating the cluster to which each response is allocated}
#' @return \item{postProbs}{the posterior group membership probabilities}
#' @return \item{Y}{the response matrix}
#' @return \item{X}{the standardized explanatory matrix}
#' @return \item{A}{the supplementary explanatory variables}
#' @export
#'
#' @examples \dontrun{
#'
#' library(rmSCGLR)
#'
#' # load sample data
#' data <- genus
#' data <- as.data.frame(apply(data, 2, as.numeric ))
#'
#' # get variable names from dataset
#' n <- names(data)
#' ny <- n[grep("^gen",n)]
#' nx <- n[-grep("^gen",n)]
#' na <- c("geology")
#' nx <- nx[!nx%in%c("geology", "surface")]
#'
#' # build multivariate formula
#' form <- multivariateFormula(Y = ny, X = nx, A = na)
#'
#' # define family
#' fam <- rep("poisson", length(ny))
#'
#' # run function
#' H <- c(2,2)
#' met <- methodSR.RMSCGLR(l=4, s=0.1, t=0.4)
#' res <- ResponseMixtureSCGLR(formula=form, data=data, H=H,
#'                            family=fam, method=met, offset = data$surface)
#' }
#'

ResponseMixtureSCGLR <- function(formula, data, H=c(2,2), family, size = NULL,
                                 offset = NULL, subset = NULL,
                                 crit = list(), method = methodSR.RMSCGLR()){

  if(!inherits(formula,"MultivariateFormula"))
    formula <- multivariateFormula(formula,data=data)

  if(length(formula$X)>1)
    stop("SCGLR Response Mixture deals with only one theme of covariates")

  additional <- formula$additional

  # check data
  if(!inherits(data, "data.frame"))
    stop("data must be compatible with a data.frame!")
  data_size <- nrow(data)

  # check crit
  crit <- do.call("critConvergence", crit)

  # extract left-hand side (Y)
  # Y is a formula of the form ~...
  if(length(formula)[1] != 1)
    stop("Left hand side part of formula (Y) must have ONE part!")
  terms_Y <- stats::terms(formula, lhs=1, rhs=0)
  Y_vars <- all.vars(terms_Y)

  # check family
  if(!is.character(family)||any(!(family%in%c("gaussian","poisson","bernoulli","binomial"))))
    stop("Expecting character vector of gaussian, poisson, bernoulli or binomial for family")
  if(!(length(family)%in%c(1,length(Y_vars))))
    stop("Length of family must be equal to one or number of Y variables")
  if(length(family)==1)
    family <- rep(family,length(Y_vars))

  # check and preprocess size parameter
  binomial_count <- sum(family=="binomial")
  if(!is.null(size)&&(binomial_count>0)) {
    if(is.vector(size)) {
      if(sum(family=="binomial")>1)
        message("Assuming that each binomial variable has same size!")
      size <- matrix(size, data_size, binomial_count)
    } else {
      size <- as.matrix(size)
    }
  }

  # check and preprocess offset parameter
  poisson_count <- sum(family=="poisson")
  if(!is.null(offset)&&(poisson_count>0)) {
    if(is.vector(offset)) {
      offset <- matrix(offset, data_size, poisson_count)
    } else {
      offset <- as.matrix(offset)
    }
  }

  # check H (number of components to keep per group)
  H <- as.integer(H)
  if(any(H<=0))
    stop("H must be a vector of positive integers")
  if((length(H) == 1) & (method$t > 0))
    stop("If H is an integer, t must be equal 0")

  # extract additional variables (A)
  # A is a formula of the form ~...
  if(additional) {
    terms_A <- stats::terms(formula, lhs=0, rhs=length(formula)[[2]])
  } else {
    terms_A <- NULL
  }
  # extract vars
  data_vars <- names(data)
  terms_X <- stats::terms(formula, lhs=0, rhs=length(formula)[[1]])
  X_vars <- all.vars(terms_X)#unique(unlist(lapply(theme_X, all.vars)))
  A_vars <- all.vars(terms_A)
  YXA_vars <- unique(c(Y_vars, X_vars, A_vars))

  # check if all variables can be found within data
  missing_vars <- YXA_vars[!YXA_vars %in% data_vars]
  if(length(missing_vars))
    stop("Some variable(s) where not found in data! '", paste(missing_vars, collapse="', '"),"'")

  # check that Y and X+A variable do not overlap
  error_vars <- Y_vars[Y_vars %in% c(X_vars, A_vars)]
  if(length(error_vars))
    stop("LHS and RHS variables must be different! '", paste(error_vars, collapse="', '"),"'")

  # check that A variables do not overlap with X
  error_vars <- A_vars[A_vars %in% X_vars]
  if(length(error_vars))
    stop("Additional variables must be different of theme variables! '", paste(error_vars, collapse="', '"),"'")

  # build data
  data <- data[,YXA_vars]

  # subset it if needed
  if(!is.null(subset)) {
    data <- data[subset,]
    offset <- offset[subset,]
    size <- size[subset,]
    weights <- weights[subset,]
    data_size <- nrow(data)
  }

  #***************#
  #initialisation #----
  #***************#
  Y <- as.matrix(data[,Y_vars])
  X <- as.matrix(data[,X_vars])
  n <- nrow(X)
  X <- apply(X, 2, FUN =  function(x) return(wtScale(x=x, w=1/n)) )
  if(!is.null(terms_A)) A <- as.matrix(data[,A_vars]) else A <- NULL

  K <- length(Y_vars)
  G <- length(H)
  if(G > K)
    stop("The number of groups is greater than the number of responses")

  p <- ncol(X)
  n <- nrow(X)
  if(!is.null(terms_A)) r <- ncol(A)

  # initialisation de Z et de W, de alpha, de u et de a
  if (!is.null(offset))
    loffset <- log(offset)

  ### Initialization, Z working variables
  mu0 <- apply(Y, 2, function(x) mean(x))
  mu0 <- matrix(mu0, n, K, byrow = TRUE)
  Z1 <- Y

  if ("bernoulli" %in% family) {
    tmu0 <- mu0[, family == "bernoulli"]
    Z1[, family == "bernoulli"] <- log(tmu0/(1 - tmu0)) +
      (Y[, family == "bernoulli"] - tmu0)/(tmu0 * (1 - tmu0))
  }
  if ("binomial" %in% family) {
    tmu0 <- mu0[, family == "binomial"]
    Z1[, family == "binomial"] <- log(tmu0/(1 - tmu0)) +
      (Y[, family == "binomial"] - tmu0)/(tmu0 * (1 - tmu0))
  }
  if ("poisson" %in% family) {
    tmu0 <- mu0[, family == "poisson"]
    if (is.null(offset)) {
      Z1[, family == "poisson"] <- log(tmu0) + (Y[, family == "poisson"] - tmu0)/tmu0
    } else {
      Z1[, family == "poisson"] <- log(tmu0) - loffset + (Y[, family == "poisson"] - tmu0)/tmu0
    }
  }
  if ("gaussian" %in% family) {
    Z1[, family == "gaussian"] <- Y[, family == "gaussian"]
  }
  Z <- list()
  for(g in 1:G) Z[[g]] <- Z1
  rm(Z1)


  # initialisation des groupes
  quanti <- apply(Y, 2, FUN = function(x) x <- x+rnorm(n=length(x), sd=0.01))
  tree <- ClustOfVar::hclustvar(X.quanti = quanti)
  clusters <- ClustOfVar::cutreevar(obj = tree, k = G)$cluster
  clusters <- clusters[Y_vars]
  rm(tree)

  #initialisation des u et f
  U <- list()
  FF <- list()
  for(g in 1:G){
    ind <- which(clusters == g)
    U[[g]] <- cbind(pls::plsr(formula = Z[[g]][,ind]~X, ncomp = H[g])$loading.weights)
    FF[[g]] <- cbind(X%*%U[[g]])
  }
  # initialisation des eta[[g]]

  eta.old <- list()
  for(g in 1:G){
    if (is.null(terms_A)) {
      reg <- cbind(1, FF[[g]])
      sol <- solve(crossprod(reg), crossprod(reg, Z[[g]]))
      eta.old[[g]] <- apply(sol, 2, function(x) x[1] + FF[[g]] %*% x[2:(H[g]+1)])
    } else {#browser()
      reg <- cbind(1, A, FF[[g]])
      sol <- solve(crossprod(reg), crossprod(reg, Z[[g]]))
      eta.old[[g]] <- apply(sol, 2, function(x) x[1] + A %*% x[2:(r + 1)] + FF[[g]] %*% x[(r+2):(r+H[g]+1)])
    }
  }
  # Update initialization of Z and initialization of W Z <- eta[[g]]
  muinf <- 1e-05
  W <- list()
  for(g in 1:G){
    W[[g]] <- matrix(1, n, K)

    if ("bernoulli" %in% family) {
      etainf <- log(muinf/(1 - muinf))
      indinf <- 1 * (eta.old[[g]][, family == "bernoulli"] < etainf)
      eta.old[[g]][, family == "bernoulli"] <- eta.old[[g]][, family == "bernoulli"] *
        (1 - indinf) + etainf * indinf
      indsup <- 1 * (eta.old[[g]][, family == "bernoulli"] > -etainf)
      eta.old[[g]][, family == "bernoulli"] <- eta.old[[g]][, family == "bernoulli"] * (1 - indsup) -
        etainf * indsup
      mu <- exp(eta.old[[g]][, family == "bernoulli"])/(1 + exp(eta.old[[g]][, family == "bernoulli"]))
      Z[[g]][, family == "bernoulli"] <- eta.old[[g]][, family == "bernoulli"] +
        (Y[, family == "bernoulli"] - mu)/(mu * (1 - mu))
      W[[g]][, family == "bernoulli"] <- mu * (1 - mu)
    }
    if ("binomial" %in% family) {
      etainf <- log(muinf/(1 - muinf))
      indinf <- 1 * (eta.old[[g]][, family == "binomial"] < etainf)
      eta.old[[g]][, family == "binomial"] <- eta.old[[g]][, family == "binomial"] * (1 - indinf) +
        etainf * indinf
      indsup <- 1 * (eta.old[[g]][, family == "binomial"] > -etainf)
      eta.old[[g]][, family == "binomial"] <- eta.old[[g]][, family == "binomial"] * (1 - indsup) -
        etainf * indsup
      mu <- exp(eta.old[[g]][, family == "binomial"])/(1 + exp(eta.old[[g]][, family == "binomial"]))
      Z[[g]][, family == "binomial"] <- eta.old[[g]][, family == "binomial"] +
        (Y[, family == "binomial"] - mu)/(mu * (1 - mu))
      W[[g]][, family == "binomial"] <- mu * (1 - mu) * size
    }
    if ("poisson" %in% family) {
      etainf <- log(muinf)
      indinf <- 1 * (eta.old[[g]][, family == "poisson"] < etainf)
      eta.old[[g]][, family == "poisson"] <- eta.old[[g]][, family == "poisson"] * (1 - indinf) +
        etainf * indinf
      if (is.null(offset)) {
        mu <- exp(eta.old[[g]][, family == "poisson"])
      } else {
        mu <- exp(eta.old[[g]][, family == "poisson"] + loffset)
      }
      Z[[g]][, family == "poisson"] <- eta.old[[g]][, family == "poisson"] +
        (Y[, family == "poisson"] - mu)/mu
      W[[g]][, family == "poisson"] <- mu
    }
    if ("gaussian" %in% family) {
      Z[[g]][, family == "gaussian"] <- Y[, family == "gaussian"]
    }
    W[[g]] <- apply(W[[g]], 2, function(x) x/sum(x))
  }

  # initialisation des probas a priori
  pg <- rep(0, G)
  for(g in 1:G){
    ind <- which(clusters == g)
    pg[g] <- length(ind)/K
  }

  f_old <- FF
  pg.old <- pg

  tol1 <- rep(Inf, G)
  tol2 <- rep(Inf, G)
  tol3 <- rep(Inf, G)

  i <- 1
  while((any(c(tol1,tol2,tol3)>method$epsilon)) && (i<method$maxiter))
  {
    #****#
    # EM #----
    #****#
    res.EM <- EM.RMSCGLR(Y=Y,A=A,Z=Z,W=W,eta=eta.old,pg=pg,H=H,FF=FF,
                         family=family,size=size,offset=offset,crit=crit)
    pg.new <- res.EM$pg
    postProbs <- res.EM$postProbs
    eta.new <- res.EM$eta
    W <- res.EM$W
    Z <- res.EM$Z
    #******#
    # PING #----
    #******#
    for(g in 1:G){
      for(h in 1:H[g]){
        if(h==1){
          U[[g]][,h] <- ping.RMSCGLR(Z=Z[[g]],X=X,AX=A,W=W[[g]],F=NULL,u=U[[g]][,h],E=FF[-g],
                                     method=method,poids=postProbs[,g])
        } else {
          u <- U[[g]][,h]
          F <- cbind(FF[[g]][,1:(h-1)])
          C <- crossprod(X, F)
          u <- u - C %*% solve(crossprod(C), crossprod(C,u))
          u <- c(u/sqrt(sum(u^2)))
          U[[g]][,h] <- ping.RMSCGLR(Z=Z[[g]],X=X,AX=A,W=W[[g]],F=F,u=u,E=FF[-g],
                                     method=method,poids=postProbs[,g])
        }
        FF[[g]][,h] <- cbind(X%*%(U[[g]][,h]))
      }
      f_old[[g]] <- apply(f_old[[g]],2,function(x) x/(sqrt(c(crossprod(x)/n))))
      f_new <- apply(FF[[g]],2,function(x) x/(sqrt(c(crossprod(x)/n))))
      tol1[g] <- abs(sum(1 - diag(crossprod(f_old[[g]], f_new)/n)^2))
      f_old[[g]] <- FF[[g]]

      tol3[g] <- mean((eta.old[[g]]-eta.new[[g]])^2)
      eta.old[[g]] <- eta.new[[g]]
    }
    tol2 <- abs(pg.old-pg.new)
    pg.old <- pg.new
    i <- i+1
  }
  #****************#
  # Classification #----
  #****************#
  res.EM <- EM.RMSCGLR(Y=Y,A=A,Z=Z,W=W,eta=eta.new,pg=pg,H=H,FF=FF,
                       family=family,size=size,offset=offset,crit=crit)
  postProbs <- res.EM$postProbs
  cluster <- as.data.frame(rbind(apply(postProbs, 1, which.max)))
  colnames(cluster) <- colnames(Y)
  row.names(cluster) <- ""

  #*****#
  # out #----
  #*****#
  UU <- cbind()
  comp <- cbind()
  names.comp <- c()
  names.u <- c()
  coef <- res.EM$sol
  for(g in 1:G){
    UU <- cbind(UU, U[[g]])
    comp <- cbind(comp, FF[[g]])
    names.u <- c(names.u, paste("G",g,"_U",1:H[g], sep = ""))
    names.comp <- c(names.comp, paste("G",g,"_SC",1:H[g], sep = ""))
    names.coef <- c("intercept", A_vars, paste("comp", 1:H[g], sep = ""))
    row.names(coef[[g]]) <- names.coef
  }
  colnames(UU) <- names.u
  colnames(comp) <- names.comp
  postProbs <- res.EM$postProbs
  colnames(postProbs) <- paste("G", 1:G, sep = "")
  row.names(postProbs) <- colnames(Y)
  pg <- as.data.frame(rbind(res.EM$pg))
  colnames(pg) <- paste("G", 1:G, sep = "")
  row.names(pg) <- ""

  out <- list(U = UU,
              comp = comp,
              eta = res.EM$eta,
              coef = coef,
              pg = pg,
              cluster = cluster,
              postProbs = postProbs,
              Y = Y,
              X = X,
              A = A)

  class(out) <- "rmSCGLR"

  return(out)
}


EM.RMSCGLR <- function(Y, A, Z, W, eta, pg, H, FF, family,
                       size=NULL, offset=NULL, crit)
{
  G <- length(pg)
  K <- ncol(Y)
  n <- nrow(Y)
  if(!is.null(A)) r <- ncol(cbind(A)) else r <- 0
  if (!is.null(offset)) loffset <- log(offset)
  muinf <- 1e-05
  step <- 1
  pg.old <- pg
  pg.dif <- pg+1
  eta.old <- eta
  eta.new <- eta
  tol <- rep(Inf, G)
  sol <- list()

  while(step < 10){
    #**********#
    # step E   #----
    #**********#
    postProbs <- matrix(data = 0, nrow = K, ncol = G)
    for(g in 1:G){
      if ("gaussian" %in% family){
        postProbs[family == "gaussian",g] <- log(pg[g])+
          apply(dnorm(Y[,family == "gaussian"], eta.old[[g]][,family == "gaussian"], log = TRUE),2,sum)
      }
      if ("poisson" %in% family) {
        if (is.null(offset)) {
          mu <- exp(eta.old[[g]][, family == "poisson"])
        } else {
          mu <- exp(eta.old[[g]][, family == "poisson"] + loffset)
        }
        postProbs[family == "poisson",g] <- log(pg[g])+
          apply(dpois(Y[,family == "poisson"], mu, log = TRUE),2,sum)
      }
      if ("bernoulli" %in% family) {
        mu <- exp(eta.old[[g]][, family == "bernoulli"])/(1 + exp(eta.old[[g]][, family == "bernoulli"]))
        postProbs[family == "bernoulli",g] <- log(pg[g])+
          apply(dbinom(Y[,family == "bernoulli"], 1, mu, log = TRUE),2,sum)
      }
      if ("binomial" %in% family) {
        mu <- exp(eta.old[[g]][, family == "binomial"])/(1 + exp(eta.old[[g]][, family == "binomial"]))
        postProbs[family == "binomial",g] <- log(pg[g])+log(factorial(size))-
          apply(dbinom(Y[,family == "binomial"], size, mu, log = TRUE),2,sum)
      }
    }
    postProbs <- exp(postProbs - apply(postProbs, 1, sumlogs))+1e-100
    if(step < 5){
      alpha <- (1-0.8*G)/(0.8*(2-G)-1)
      postProbs <- (2*alpha*postProbs-alpha+1)/(2*alpha-alpha*G+G)
    }

    #**********#
    # step M   #----
    #**********#
    pg.new <- colMeans(postProbs)
    for(g in 1:G){
      sol[[g]] <- matrix(1, nrow = 1+r+H[g], ncol = K)
      if (is.null(A)) {
        reg <- cbind(1, FF[[g]])
        for (j in seq(K)) {
          sol[[g]][,j] <- solve(crossprod(reg, W[[g]][,j]*reg), crossprod(reg, W[[g]][,j]*Z[[g]][,j]))
          eta.new[[g]][,j] <- sol[[g]][1,j] + FF[[g]] %*% sol[[g]][2:(H[g]+1),j]
        }
      } else {
        reg <- cbind(1, A, FF[[g]])
        for (j in seq(K)) {
          sol[[g]][,j] <- solve(crossprod(reg, W[[g]][,j]*reg), crossprod(reg, W[[g]][,j]*Z[[g]][,j]))
          eta.new[[g]][,j] <- sol[[g]][1,j] + cbind(A) %*% sol[[g]][2:(r+1),j] + FF[[g]] %*% sol[[g]][(r+2):(r+H[g]+1),j]
        }
      }

      if ("bernoulli" %in% family) {
        etainf <- log(muinf/(1 - muinf))
        indinf <- 1 * (eta.new[[g]][, family == "bernoulli"] < etainf)
        eta.new[[g]][, family == "bernoulli"] <- eta.new[[g]][, family == "bernoulli"] *
          (1 - indinf) + etainf * indinf
        indsup <- 1 * (eta.new[[g]][, family == "bernoulli"] > -etainf)
        eta.new[[g]][, family == "bernoulli"] <- eta.new[[g]][, family == "bernoulli"] * (1 - indsup) -
          etainf * indsup
        mu <- exp(eta.new[[g]][, family == "bernoulli"])/(1 + exp(eta.new[[g]][, family == "bernoulli"]))
        Z[[g]][, family == "bernoulli"] <- eta.new[[g]][, family == "bernoulli"] +
          (Y[, family == "bernoulli"] - mu)/(mu * (1 - mu))
        W[[g]][, family == "bernoulli"] <- mu * (1 - mu)
      }
      if ("binomial" %in% family) {
        etainf <- log(muinf/(1 - muinf))
        indinf <- 1 * (eta.new[[g]][, family == "binomial"] < etainf)
        eta.new[[g]][, family == "binomial"] <- eta.new[[g]][, family == "binomial"] * (1 - indinf) +
          etainf * indinf
        indsup <- 1 * (eta.new[[g]][, family == "binomial"] > -etainf)
        eta.new[[g]][, family == "binomial"] <- eta.new[[g]][, family == "binomial"] * (1 - indsup) -
          etainf * indsup
        mu <- exp(eta.new[[g]][, family == "binomial"])/(1 + exp(eta.new[[g]][, family == "binomial"]))
        Z[[g]][, family == "binomial"] <- eta.new[[g]][, family == "binomial"] +
          (Y[, family == "binomial"] - mu)/(mu * (1 - mu))
        W[[g]][, family == "binomial"] <- mu * (1 - mu) * size
      }
      if ("poisson" %in% family) {
        etainf <- log(muinf)
        indinf <- 1 * (eta.new[[g]][, family == "poisson"] < etainf)
        eta.new[[g]][, family == "poisson"] <- eta.new[[g]][, family == "poisson"] * (1 - indinf) +
          etainf * indinf
        if (is.null(offset)) {
          mu <- exp(eta.new[[g]][, family == "poisson"])
        } else {
          mu <- exp(eta.new[[g]][, family == "poisson"] + loffset)
        }
        Z[[g]][, family == "poisson"] <- eta.new[[g]][, family == "poisson"] +
          (Y[, family == "poisson"] - mu)/mu
        W[[g]][, family == "poisson"] <- mu
      }
      if ("gaussian" %in% family) {
        Z[[g]][, family == "gaussian"] <- Y[, family == "gaussian"]
      }
      W[[g]] <- apply(W[[g]], 2, function(x) x/sum(x))

      tol[g] <- mean((eta.old[[g]]-eta.new[[g]])^2)
      eta.old[[g]] <- eta.new[[g]]
    } # end loop on groups

    pg.dif <- abs(pg.old-pg.new)
    pg.old <- pg.new
    step <- step+1
  } # end loop EM

  return(list(pg=pg.old, postProbs=postProbs, W=W, eta=eta.new, Z=Z, sol=sol))
}
