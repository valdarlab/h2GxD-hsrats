
# "sample_from_approx_prior" generates simulated values from the
# (approximate) prior distribution for the heritability model
# with a user-supplied number of variance terms (both genetic
# and non-genetic random effects).


# Arguments are
#      maxVal: a large maximum value to truncate the prior 
#              density in order to make the distribution proper 
#              (e.g. 100,000)
#      nBins: the number of bins to partition the prior density
#      numVar: the number of random effects (genetic + non-genetic)
#      sh: shape parameter for inverse Gamma prior
#      rt: rate parameter for inverse Gamma prior

sample_from_approx_prior<-function(maxVal,nBins,numVar,sh,rt){
  vals<-seq(0,maxVal,length=nBins+1)
  midpts<-(vals[-1]+vals[-length(vals)])/2
  prob_vals<-(midpts^(-(sh)-1)*exp(-rt/midpts))/sum(midpts^(-(sh)-1)*exp(-rt/midpts))
  var_vec<-matrix(NA,1e6,numVar+1)
  for(j in 1:ncol(var_vec)){
    var_vec[,j]<-sample(midpts,1e6,replace=TRUE,prob=prob_vals)
  }
  total_var<-rowSums(var_vec)
  return(var_vec/total_var)
}
