

#'
#' @export
#'
#' @examples
#' p <- 10
#' k <- 5
#' eps <- 0.5
#' Sigma0 <- make_expbanded(p,k,eps=eps)
#' X <- generate_data(12,Sigma0)
#' kvec <- 4:6
#' epsvec <- c(0.1,0.01)
#' library(future)
#' plan(multisession, workers = 10)
#' z <- looCV_bandPPP(X,kvec,epsvec)
#'
looCV_bandPPP <- function(X,kvec,epsvec,mcmc.num=100,hyper=NULL){
  library(magrittr)
  library(Matrix)
  tuneparams <- expand.grid(k=kvec,eps=epsvec)
  furrr::future_map(1:dim(tuneparams)[1],
                    ~mean(loo_pred_loglik(X=X,k=tuneparams$k[.x],
                                          eps=tuneparams$eps[.x],
                                          mcmc.num=mcmc.num,hyper=hyper)),
                    .options = furrr::furrr_options(seed = TRUE)) %>% unlist %>%
    cbind(tuneparams,logpdf=.) %>% dplyr::arrange(desc(logpdf))
}

#'
#' @export
#'
#' @examples
#' p <- 10
#' k <- 5
#' eps <- 0.5
#' Sigma0 <- make_expbanded(p,k,eps=eps)
#' X <- generate_data(12,Sigma0)
#' kvec <- 4:6
#' epsvec <- c(0.1,0.01)
#' library(future)
#' plan(multisession, workers = 10)
#' z <- looCV_dualbandPPP(X,kvec,epsvec)
#'
looCV_dualbandPPP <- function(X,kvec,epsvec,mcmc.num=100,hyper=NULL){
  library(magrittr)
  library(Matrix)
  tuneparams <- expand.grid(k=kvec,eps=epsvec)
  furrr::future_map(1:dim(tuneparams)[1],
                    ~mean(loo_pred_loglik_dual(X=X,k=tuneparams$k[.x],
                                          eps=tuneparams$eps[.x],
                                          mcmc.num=mcmc.num,hyper=hyper)),
                    .options = furrr::furrr_options(seed = TRUE)) %>% unlist %>%
    cbind(tuneparams,logpdf=.) %>% dplyr::arrange(desc(logpdf))
}
