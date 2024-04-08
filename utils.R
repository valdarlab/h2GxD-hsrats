### generic functions for HS rat project 

ensure_directory <- function(directory){
  if(!dir.exists(directory)){
    dir.create(directory, recursive=TRUE);
  }
}

ifow <- function(test, yes, no){
  if(test){return(yes)} 
  no
}

igrep <- function(pattern, x, ..., value=FALSE, logical=TRUE){
# pass through method for grep that can return match as an indicator
# vector. From WVmisc package.
    if (!value & logical) {
        indices <- grep(pattern, x, value=value, ...)
        return (1:length(x) %in% indices)
    }
    grep(pattern, x, value=value, ...)
}


load.to.object <- function(fname){
  eh <- new.env()
  load(fname, envir=eh)
  x <- get(ls(eh)[1], envir=eh)
}

subsample_chain <- function(x, burnin, thin){
  ind <- seq(burnin, ncol(x$sigma_sq), thin)
  return(list(alpha = x$alpha[,,ind,drop=FALSE], u = x$u[,,ind,drop=FALSE], sigma_sq = x$sigma_sq[,ind,drop=FALSE],
              sigma_sq_lambda = x$sigma_sq_lambda[,ind,drop=FALSE], sigma_sq_beta = x$sigma_sq_beta[,,ind,drop=FALSE]))
}
