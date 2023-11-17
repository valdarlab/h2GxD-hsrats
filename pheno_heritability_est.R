library(coda)
#library(pheatmap)
#library(regress)
library(MASS)
library(Matrix)

source('bayesian_LMM.R')
source('utils.R')
source('cmdline.R')

#-----------------------------------Setup--------------------------------------#
datapath <- cmdline.option("datapath")
matpath <- cmdline.option("matpath")
sex <- cmdline.option("sex")
pheno <- cmdline.option("pheno")

data <- read.csv(datapath)
data$cohort_f <- factor(data$Cohort)

G <- as.matrix(read.table(matpath, sep=" ", check.names = FALSE)) 

num_chains  <- cmdline.integer("num_chains", default=3)
num_samples <- cmdline.integer("num_samples", default=100000)
burnin      <- cmdline.integer("burnin", default=1000)
thin        <- cmdline.integer("thin", default=10) 


#----------------------------------Analysis------------------------------------#
ensure_directory("derived_data/posterior_samples")
outpath = paste("derived_data/posterior_samples/MCMC_", sex, "_", pheno, sep="")

fixed_form <- as.formula(ifelse(sex == "Male", paste(pheno, "~ Litter_Size + Study + z"), 
                                paste(pheno, "~ Litter_Size + z")))
MCMC_results <- MCMC_heritability_mixed_genetic_and_non_genetic(fixed_form = fixed_form,
                    random_form = ~ cohort_f-1, genetic_form = ~ z, dat = data, 
                    G = G, num_chains = num_chains, num_samples = num_samples, 
                    hyperpars_genetic = c(-1,0), hyperpars_random = c(-1,0), 
                    hyperpars_error = c(-1,0))
save(MCMC_results, file = outpath) # for use by SavageDicketGxT.Rmd

# thin 
MCMC_results_subsample <- subsample_chain(MCMC_results, burnin = burnin, thin = thin)

# Calculate median, standard deviation, HPD interval 
median_subsample_list <- mcmc.list()
for(cc in 1:3){
  intercept_coef<-MCMC_results_subsample$alpha[1,cc,]
  litter_size_coef<-MCMC_results_subsample$alpha[2,cc,]
  if(sex == "Male"){
    StudyThesis_coef<-MCMC_results_subsample$alpha[3,cc,]
    diet_coef<-MCMC_results_subsample$alpha[4,cc,]
  } else if (sex == "Female"){
    diet_coef<-MCMC_results_subsample$alpha[3,cc,]
  }
  heritability_genetic_component_a<-MCMC_results_subsample$sigma_sq_beta[1,cc,]/(MCMC_results_subsample$sigma_sq_beta[1,cc,]+MCMC_results_subsample$sigma_sq_beta[2,cc,]+MCMC_results_subsample$sigma_sq_lambda[cc,]+MCMC_results_subsample$sigma_sq[cc,])
  heritability_genetic_component_b<-MCMC_results_subsample$sigma_sq_beta[2,cc,]/(MCMC_results_subsample$sigma_sq_beta[1,cc,]+MCMC_results_subsample$sigma_sq_beta[2,cc,]+MCMC_results_subsample$sigma_sq_lambda[cc,]+MCMC_results_subsample$sigma_sq[cc,])
  heritability_cohort<-MCMC_results_subsample$sigma_sq_lambda[cc,]/(MCMC_results_subsample$sigma_sq_beta[1,cc,]+MCMC_results_subsample$sigma_sq_beta[2,cc,]+MCMC_results_subsample$sigma_sq_lambda[cc,]+MCMC_results_subsample$sigma_sq[cc,])
  
  temp_obj<-cbind(litter_size_coef,diet_coef,heritability_genetic_component_a,heritability_genetic_component_b,heritability_cohort)
  median_subsample_list[[cc]]<-mcmc(temp_obj)
}

hpdint <- HPDinterval(median_subsample_list[1])

#--------------------------------Output stats----------------------------------#
h2_stats <- data.frame(h2_compa_median = median(heritability_genetic_component_a),
                       h2_compa_sd = sd(heritability_genetic_component_a),
                       h2_compa_lower = hpdint[[1]][3,1],
                       h2_compa_upper = hpdint[[1]][3,2],
                       h2_compb_median = median(heritability_genetic_component_b),
                       h2_compb_sd = sd(heritability_genetic_component_b),
                       h2_compb_lower = hpdint[[1]][4,1],
                       h2_compb_upper = hpdint[[1]][4,2],
                       h2_cohort_lower = hpdint[[1]][5,1],
                       h2_cohort_upper = hpdint[[1]][5,2],
                       litsz_coef_lower = hpdint[[1]][1,1],
                       litsz_coef_upper = hpdint[[1]][1,2],
                       diet_coef_lower = hpdint[[1]][2,1],
                       diet_coef_upper = hpdint[[1]][2,2],
                       probability = attributes(hpdint[[1]])$Probability,
                       row.names = paste(sex, pheno, sep="_"))

ensure_directory("derived_data/heritability_stats")
outpath = paste("derived_data/heritability_stats/", sex, "_", pheno, ".csv", sep="")
write.csv(h2_stats, file = outpath)
