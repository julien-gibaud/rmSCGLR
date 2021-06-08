#' @ Projected Iterated Normed Gradient
#'
#' @param Z matrix of the working variables
#' @param X matrix of the normalized covariates
#' @param AX matrix of the suplementary covariates
#' @param W matrix of weights
#' @param F matrix of components
#' @param u loading vector
#' @param method structural relevance criterion
#' @param poids weights on the working variables
#' @param E set of components
#'
#' @return the updated loading vectors
#' @export
#'

ping.RMSCGLR <- function(Z,X,AX,W,F,u,method,poids,E) {
  u_cur <- u
  h_cur <- hFunct.RMSCGLR(Z=Z,X=X,AX=AX,W=W,u=u_cur,method=method,poids=poids,E=E,F=F)$h
  u_new <- update_u.RMSCGLR(Z=Z,X=X,AX=AX,W=W,F=F,u=u_cur,method=method,poids=poids,E=E)
  h_new <- hFunct.RMSCGLR(Z=Z,X=X,AX=AX,W=W,u=u_new,method=method,poids=poids,E=E,F=F)$h
  ing_iter <- 0
  while((abs((h_new-h_cur)/h_cur)>method$epsilon)&(ing_iter<=method$maxiter)){
    u_cur <- u_new
    h_cur <- h_new
    u_new <- update_u.RMSCGLR(Z=Z,X=X,AX=AX,W=W,F=F,u=u_cur,method=method,poids=poids,E=E)
    h_new <- hFunct.RMSCGLR(Z=Z,X=X,AX=AX,W=W,u=u_new,method=method,poids=poids,E=E,F=F)$h
    ing_iter <- ing_iter+1
  }
  return(u_new)
}

update_u.RMSCGLR <- function(Z,X,AX,W,F,u,method,poids,E){
  out_h <- hFunct.RMSCGLR(Z=Z,X=X,AX=AX,W=W,u=u,method=method,poids=poids,E=E,F=F)
  if(is.null(F)) {
    m <- out_h$gradh/sqrt(sum(out_h$gradh^2))
  } else {
    C <- (crossprod(X,F))#/nrow(X))
    proj_C_ortho <- out_h$gradh - C%*%solve(crossprod(C),crossprod(C,out_h$gradh))
    m <- c(proj_C_ortho / sqrt(sum(proj_C_ortho^2)))
  }
  h_m <- hFunct.RMSCGLR(Z=Z,X=X,AX=AX,W=W,u=m,method=method,poids=poids,E=E,F=F)$h
  h_u <- out_h$h
  k <- 1
  # browser()
  while((h_m<h_u)&(k<method$bailout)){
    m <- c(u)+m
    m <- m/sqrt(sum(m^2))
    h_m <- hFunct.RMSCGLR(Z=Z,X=X,AX=AX,W=W,u=m,method=method,poids=poids,E=E,F=F)$h
    k <- k+1
  }
  if(k>method$bailout) print("ARGH")
  u_new <- m
  return(u_new)
}


hFunct.RMSCGLR <- function(Z,X,AX,W,u,method,poids,E,F)
{
  f <- c(X%*%u)
  psi <- 0
  gradpsi <- rep(0,length(u))
  if(!is.null(AX)){
    for(k in 1:ncol(Z)){
      # browser()
      AXtWkAX <- crossprod(AX,W[,k]*AX)
      projWkfAX <- c(AX%*%solve(AXtWkAX,crossprod(AX,W[,k]*f)))
      projWkforthoAX <- f - projWkfAX
      Zk <- wtScale(Z[,k],W[,k])#Z[,k] - sum(W[,k]*Z[,k])
      WZk <- W[,k]*Zk
      projWkZAX <- AX%*%solve(AXtWkAX,crossprod(AX,WZk))
      #calcul de psi
      scalsqpfZ <- sum(c(projWkforthoAX)*WZk)^2
      scalsqpfpf <- sum(c(projWkforthoAX)^2*W[,k])
      term1psi <- sum(scalsqpfZ/(scalsqpfpf))
      term2psi <- sum(WZk*projWkZAX)
      psi <- psi+poids[k]*(term1psi+term2psi)
      #calcul de grad de psi
      PiorthoPrimeWkZ <- WZk -  W[,k]*AX%*%solve(AXtWkAX,crossprod(AX,WZk))
      XprimeprojorthoWZ <- crossprod(X,PiorthoPrimeWkZ)
      term1 <- c(XprimeprojorthoWZ%*%crossprod(XprimeprojorthoWZ,u))/(scalsqpfpf)

      WprojWkOrthof <- W[,k]*projWkforthoAX
      term2 <-  scalsqpfZ*c(crossprod(X,WprojWkOrthof))/(scalsqpfpf^2)
      gradpsi <- gradpsi +poids[k]*(term1-term2)
    }
    gradpsi <- 2*gradpsi
  }else{
    for(k in 1:ncol(Z)){
      Zk <- wtScale(Z[,k],W[,k])#Z[,k] - sum(W[,k]*Z[,k])
      WZk <- W[,k]*Zk
      scalsqpfZ <- sum(c(f)*WZk)^2
      scalsqpfpf <- sum(c(f)^2*W[,k])
      #calcul de psi
      psi <- psi+poids[k]*sum(scalsqpfZ/(scalsqpfpf))
      #calcul de grad de psi
      XprimeWZ <- crossprod(X,WZk) #X'W_k Z_k
      term1 <- c(XprimeWZ%*%crossprod(XprimeWZ,u))/(scalsqpfpf)
      term2 <-  scalsqpfZ*c(crossprod(X,W[,k]*f))/(scalsqpfpf^2)
      gradpsi <- gradpsi +poids[k]*(term1-term2)
    }
    gradpsi <- 2*gradpsi
  }
  n <- nrow(X)
  # calcul phi Component Variance: cv
  if(method$phi=="cv") {
    phi <- c(crossprod(f))/n
    # calcul grad phi
    gradphi <- c(2*crossprod(X,f/n))
  } else {
    ### autre calcul de phi avec l>=1 : vpi: Variable Powered Inertia
    scalsqfX <- colSums(f*X/n)
    XtWX <- crossprod(X)/n
    phi <- (sum((scalsqfX^2)^method$l))^(1/method$l)
    # calcul de grad phi
    gradphi <- 2*phi^(1-method$l)*rowSums(XtWX%*%diag(scalsqfX)^(2*method$l-1))
  }
  if(!is.null(E) && method$t>0){
    varphi <- 0
    gradvarphi <- rep(0,length(f))
    H <- 1
    G <- length(E)
    if(!is.null(F)){
      H <- ncol(F)+1
      for(g in 1:G){
        for(i in 1:(H-1)){
          for(j in 1:ncol(E[[g]])){
            num <- c(crossprod(F[,i]/n, E[[g]][,j])^2)
            den <- sqrt(ncol(E[[g]]))*c(crossprod(F[,i])/n)*c(crossprod(E[[g]][,j])/n)
            varphi <- varphi + num/den
          }
        }
      }
    }
    for(g in 1:G){
      for(j in 1:ncol(E[[g]])){
        num <- c(crossprod(f/n, E[[g]][,j])^2)
        den <- sqrt(ncol(E[[g]]))*c(crossprod(f)/n)*c(crossprod(E[[g]][,j])/n)
        varphi <- varphi + num/den

        num1 <- c(crossprod(f/n, E[[g]][,j])*crossprod(f)/n)*c(E[[g]][,j])
        num2 <- c(crossprod(f/n, E[[g]][,j])^2)*f
        den1 <- sqrt(ncol(E[[g]]))*c(crossprod(f)/n)^2*c(crossprod(E[[g]][,j])/n)
        gradvarphi <- gradvarphi + (num1 - num2)/den1
      }
    }
    varphi <- abs(1-varphi/(sqrt(H)*length(E)))
    gradvarphi <- -2*crossprod(X/n, gradvarphi)/(sqrt(H)*length(E))

    # print(varphi)
    h <- (1-method$s-method$t)*log(psi)+method$s*log(phi)+method$t*log(varphi)
    gradh <- (1-method$s-method$t)*gradpsi/psi+method$s*gradphi/phi+method$t*gradvarphi/varphi
    return(list(h=h, gradh=gradh,psi=psi,gradpsi=gradpsi,phi=phi,gradphi=gradphi))
  }
  # calcul de h (s in R+)
  #h = log(psi)+method$s*log(phi)
  # print(log(phi))
  # print(log(psi))
  # print(log(varphi))
  #gradh=gradpsi/psi+method$s*gradphi/phi
  # calcul de h (s in [0..1])
  h <- (1-method$s)*log(psi)+method$s*log(phi)
  gradh <- (1-method$s)*gradpsi/psi+method$s*gradphi/phi
  return(list(h=h, gradh=gradh,psi=psi,gradpsi=gradpsi,phi=phi,gradphi=gradphi))
}
