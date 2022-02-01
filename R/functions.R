#functions.R

#'
#'
#'
pd_adjustment <- function(Sigma,eps=10^(-4)){
  egmin <- min(eigen(Sigma)$value)
  if(egmin<eps){
    return(Sigma + diag(-egmin+eps,dim(Sigma)[1]))
  }
  return (Sigma)
}

#'
#'
#'
pd_adjustment_Matrix <- function(Sigma,eps=10^(-4)){
  egmin <- RSpectra::eigs(Sigma, 1, sigma =
                            -RSpectra::eigs(Sigma, 1, which = "LM")$values)$values
  if(egmin<eps){
    return(Sigma + Matrix::Diagonal(dim(Sigma)[1],-egmin+eps))
  }
  return (Sigma)
}
