#simulation_setting.R


#' @export
make_bandable <- function(p,rho=0.6,alpha=0.5,eps=10^(-4)){
  Sigma=diag(1,p)
  Sigma <- (rho*(abs(row(Sigma)-col(Sigma)))^(-alpha-1))
  diag(Sigma) <- 1

  return(pd_adjustment(Sigma,eps))
}

#' @export
make_expbanded <- function(p,k,rho=0.6,alpha=1,eps=10^(-4)){
  Sigma <- make_bandable(p,rho,alpha)

  Sigma0 <- FinCovRegularization::banding(Sigma,k)
  return(pd_adjustment(Sigma0,eps))
}

#' @export
make_linearbanded <- function(p,k,eps=10^(-4)){
  Sigma=diag(1,p)
  Sigma<-pmax(1-abs(row(Sigma)-col(Sigma))/(k+1),0)

  return(pd_adjustment(Sigma,eps))
}

#' @export
make_randombanded <- function(p,k,shape=5,scale=1,eps=10^(-4),seed=1){
  set.seed(seed)
  L=matrix(rnorm(p*p,0,1),p,p)
  Lind<-lower.tri(matrix(1,p,p)) & FinCovRegularization::banding(matrix(1,p,p),k)
  L[!Lind]<-0
  diag(L)<-1
  D<-diag(MCMCpack::rinvgamma(p,shape=shape,scale=scale))
  Sigma<-L%*%D%*%t(L)

  return(pd_adjustment(Sigma,eps))
}


#' @export
generate_data <- function(n,Sigma){
  return(mvtnorm::rmvnorm(n,sigma = Sigma))
}

