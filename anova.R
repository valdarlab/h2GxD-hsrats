library(lme4qtl)
library(pbkrtest) # Kenward-Roger method 
library(MASS)
source("utils.R")

dat <- read.csv("derived_data/transformed_data_all-G3package.csv")
K <- as.matrix(read.table("derived_data/relatedness_matrix_MetBehThesis.txt", 
                          sep=" ", check.names = FALSE))

pred_terms <- c('Study', 'Litter_Size', 'I(Litter_Size^2)', 'Sex', 'Diet', 'Sex:Diet')
pred_vars <- c('Study', 'Litter_Size', 'Sex', 'Diet') 
outcome_vars <- c("HarvWeight", "RetroFat_norm", "EpiFat_norm", "OmentalFat_norm",
                  "Total_AUC", "FastGluc", "FastIns", "ClosedJunc", "OpenJunc",
                  "REST_EPISODE_COUNT_5", "MOVEMENT_EPISODE_COUNT_5",
                  "VERTICAL_EPISODE_COUNT_5", "SWIM", "CLIMB", "FLOAT")

results.F <- data.frame(matrix(nrow = length(outcome_vars), ncol = length(pred_terms)), 
                      row.names = outcome_vars)
colnames(results.F) <- pred_terms
results.p <- results.F

for (i in 1:length(outcome_vars)){
  pheno.name <- outcome_vars[i]
  
  # Create null model
  # include first predictor variable in dataframe to avoid errors with first comparison
  k <- 1 # for iterating over predictor variables / adding columns to lm data 
  moddata <- dat[,c(pheno.name, 'Cohort', 'Access_ID', pred_vars[k])]  
  formula0 <- paste0(pheno.name, " ~ (1|Cohort) + (1|Access_ID)")
  mod0 <- relmatLmer(formula0, data = na.omit(moddata), relmat = list(Access_ID = K))
  
  formula1 <- formula0
  
  for (j in 1:length(pred_terms)){
  
    formula1 <- paste0(formula1, " + ", pred_terms[j])
    mod1 <- relmatLmer(formula1, data = na.omit(moddata), relmat = list(Access_ID = K))
    
    kr <- KRmodcomp(mod1, mod0)
    results.F[i,j] <- kr$test$stat[1]
    results.p[i,j] <- kr$test$p.value[1]
    
    # Create dataset and small model for next comparison 
    if((!is.na(pred_terms[j+1])) & (pred_terms[j+1]!='I(Litter_Size^2)') & (pred_terms[j+1]!='Sex:Diet')){
      k <- k+1
      moddata <- cbind(moddata, dat[,pred_vars[k]]) # add next col
      colnames(moddata)[k+3] <- pred_vars[k]
    }
    # current large model becomes small model for the next comparison (using new dataset
    # with NAs omitted for next variable as well)
    formula0 <- formula1 
    mod0 <- relmatLmer(formula0, data = na.omit(moddata), relmat = list(Access_ID = K))
  }
}

#-----------------------------P-value adjustments------------------------------#

# For each covariate (pred_terms)/column in results.p
# Rank p-values from smallest (most significant) to largest (least significant) 
### need to keep original ranks and order by original rank 
# Assign rank value 1-15
# Calculate adjusted p-value: original p-value * (17/rank)

for (i in 1:length(pred_terms)){
  # determine ranks 
  pdf <- data.frame(og.rank = seq(1,nrow(results.p)), pval = results.p[,i])
  pdf <- pdf[order(pdf$pval),]
  pdf$rank <- seq(1, nrow(pdf))
  
  # calculate adjusted p-value 
  pdf$adjpval <- pdf$pval*(15/pdf$rank) 
  
  # return to original order 
  pdf <- pdf[order(pdf$og.rank),]
  
  # replace pvals in results.p with adjusted pvals 
  results.p[,i] <- pdf$adjpval
}

#--------------------------------Save results----------------------------------#
ensure_directory("results")
write.csv(results.p, file = "results/anova-results.csv")



