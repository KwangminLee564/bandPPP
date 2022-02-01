#'
#' @export
#'
#' @examples
#' k <- 4
#' eps <- 0.5
#' Sigma0 <- make_linearbanded(10,k,eps=eps)
#' X <- generate_data(50,Sigma0)
#' Xnew <- generate_data(1,Sigma0)[1,]
#' postsample <- bandPPP(X,k=4,mcmc.num=10)
#' postsample2 <- banddualPPP(X,k=4,mcmc.num=10)
#'
#' pred_loglik(postsample2,Xnew)
#'
pred_loglik <- function(postsample,Xnew){
  p <- dim(postsample[[1]])[1]
  logpdf <- postsample %>%
    purrr::map(~mvnfast::dmvn(Xnew,mu=rep(0,p),.x,log = TRUE)) %>%
    unlist
  matrixStats::logSumExp(logpdf) - log(length(postsample))
}

#'
#' @export
#'
#' @examples
#' k <- 4
#' eps <- 0.5
#' Sigma0 <- make_linearbanded(10,k,eps=eps)
#' X <- generate_data(50,Sigma0)
#' Xnew <- generate_data(1,Sigma0)[1,]
#' postsample <- bandPPP(X,k=4,mcmc.num=10)
#' postsample2 <- banddualPPP(X,k=4,mcmc.num=10)
#'
#' pred_loglik(postsample2,Xnew)
#'
pred_loglik2 <- function(postsample,Xnew){
  p <- dim(postsample[[1]])[1]
  postsample %>%
    purrr::map(~sum(mvnfast::dmvn(Xnew,mu=rep(0,p),.x,log = TRUE)) ) %>%
    unlist %>% mean
}


#'
#' @export
#'
#' @examples
#' k <- 4
#' eps <- 0.5
#' Sigma0 <- make_linearbanded(10,k,eps=eps)
#' X <- generate_data(50,Sigma0)
#'
#'
loo_pred_loglik <- function(X,k,eps=NULL,mcmc.num=100,hyper=NULL){
  n <- dim(X)[1]

  1:n %>% purrr::map(~pred_loglik(bandPPP(X[-.x,],k=k,eps=eps,mcmc.num=mcmc.num,
                                          hyper=hyper),X[.x,])) %>% unlist

}

#'
#' @export
#'
#' @examples
#' k <- 4
#' eps <- 0.5
#' Sigma0 <- make_linearbanded(10,k,eps=eps)
#' X <- generate_data(12,Sigma0)
#' loo_pred_loglik_dual(X,k)
#'
loo_pred_loglik_dual <- function(X,k,eps=NULL,mcmc.num=100,hyper=NULL){
  n <- dim(X)[1]

  1:n %>% purrr::map(~pred_loglik(banddualPPP(X[-.x,],k=k,mcmc.num=mcmc.num,
                                              hyper=hyper),X[.x,])) %>% unlist

}


