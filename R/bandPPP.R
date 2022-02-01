#'
#' @export
#' @importFrom magrittr %>%
#' @examples
#'
#' library(magrittr)
#' k <- 4
#' eps <- 0.1
#' Sigma0 <- make_linearbanded(10,k)
#' X <- generate_data(50,Sigma0)
#' res2 <- bandPPP(X,4,mcmc.num=100)
#'
bandPPP <- function(X,k,eps=NULL,mcmc.num=100,hyper=NULL){
  library(magrittr)
  p <- dim(X)[2]
  n <- dim(X)[1]
  if(is.null(hyper)){
    hyper <- list(nu = p+k, A = diag(1,p))
  }
  if(is.null(eps)){
    eps <- (log(k))^2 * (k+ log(p))/n
  }
  param_initpost <- list(A=t(X)%*%X + hyper$A,nu=hyper$nu + n)

  purrr::rerun(mcmc.num,CholWishart::rInvWishart(1,param_initpost$nu,
                                                 param_initpost$A)[,,1] %>%
                 Matrix::Matrix(sparse = T) %>%
                 Matrix::band(k1 = -k,k2=k) %>% pd_adjustment_Matrix(eps=eps))

}

#'
#' @export
#' @importFrom magrittr %>%
#' @examples
#'
#' library(magrittr)
#' k <- 4
#' eps <- 0.1
#' Sigma0 <- make_linearbanded(10,k)
#' X <- generate_data(50,Sigma0)
#' res <- banddualPPP(X,4,mcmc.num=100)
#'
banddualPPP <- function(X,k,eps=NULL,mcmc.num=100,hyper=NULL){
  library(magrittr)
  p <- dim(X)[2]
  n <- dim(X)[1]
  if(is.null(hyper)){
    hyper <- list(nu = p+k, A = diag(1,p))
  }
  param_initpost <- list(A=t(X)%*%X + hyper$A,nu=hyper$nu + n)

  purrr::rerun(mcmc.num,rWishart(1,param_initpost$nu,
                                 Rfast::spdinv(param_initpost$A))[,,1] %>%
                 as.matrix(Matrix::band(k1 = -k,k2=k)) %>%
                 bandLA::dualMLE(k) %>%
                 as("symmetricMatrix"))

}
