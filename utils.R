### generic functions for HS rat project 

ensure_directory <- function(directory){
  if(!dir.exists(directory)){
    dir.create(directory);
  }
}

subsample_chain <- function(x, burnin, thin){
  ind <- seq(burnin, ncol(x$sigma_sq), thin)
  return(list(alpha = x$alpha[,,ind,drop=FALSE], u = x$u[,,ind,drop=FALSE], sigma_sq = x$sigma_sq[,ind,drop=FALSE],
              sigma_sq_lambda = x$sigma_sq_lambda[,ind,drop=FALSE], sigma_sq_beta = x$sigma_sq_beta[,,ind,drop=FALSE]))
}