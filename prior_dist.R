### This script uses sampling to create prior distribution for 
### Bayes Factor analysis 

library(logspline)
source("sample_from_approx_prior.R") # Andrew's prior sampler code

prior <- c(sample_from_approx_prior(maxVal=1e5, nBins=1e6, sh=-1, rt=0, numVar=3)) # 4 million samples
# saveRDS(prior, file="priorsamples.RDS")
# prior <- readRDS("priorsamples.RDS")
lsp.prior <- logspline(prior, lbound=0, ubound=1)
saveRDS(lsp.prior, file="priorlsp.RDS")