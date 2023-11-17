library(lme4qtl)
library(pbkrtest) # Kenward-Roger method 
library(MASS)
source("utils.R")

#----------------------------------Load data-----------------------------------#
# phenotypes 
dat <- read.csv("derived_data/transformed_data_All-G3package.csv")
# kinship matrix 
K <- as.matrix(read.table("derived_data/relatedness_matrix_MaleFemale_PilotThesis.txt", 
                          sep=" ", check.names = FALSE))

#------------------------------Define variables--------------------------------#
### ADDED DIET_START_WEEK to pred_terms and pred_vars
pred_terms <- c('Study', 'Litter_Size', 'I(Litter_Size^2)', 'Diet_Start_Week', 'Sex', 'Diet', 'Sex:Diet')
pred_vars <- c('Study', 'Litter_Size', 'Diet_Start_Week', 'Sex', 'Diet') 
outcome_vars <- c("HarvWeight", "RetroFat_norm", "EpiFat_norm", "OmentalFat_norm",
                  "Total_AUC", "FastGluc", "FastIns", "Factor1", "ClosedJunc", 
                  "OpenJunc", "Factor3", "SWIM", "CLIMB", "FLOAT", "Factor4",
                  "REST_EPISODE_COUNT_5", "MOVEMENT_EPISODE_COUNT_5",
                  "VERTICAL_EPISODE_COUNT_5", "Factor2")


#-------------------------------2-way ANOVA------------------------------------#
# All data were pooled and analyzed for batch effects and age at diet start 
# within and between studies via 2-way ANOVA  
# twoway <- data.frame(phenotype = outcome_vars,
#                      diet_start_p = rep(NA, length(outcome_vars)),
#                      cohort_p = rep(NA, length(outcome_vars)))
# 
# for (i in 1:length(outcome_vars)){
#   pheno.name <- outcome_vars[i]
#   form <- paste0(pheno.name, " ~ Diet_Start_Week + Cohort")
#   aov <- anova(lm(form, data = dat))
#   twoway[i,'diet_start_p'] <- aov$`Pr(>F)`[1]
#   twoway[i,'cohort_p'] <- aov$`Pr(>F)`[2]
# }
# 
# ensure_directory("results")
# write_csv(twoway, "results/twowayanova.csv")


#--------------------------Sex-combined analysis-------------------------------#
# ANOVA results dataframe (F values and P values)
results.F <- data.frame(matrix(nrow = length(outcome_vars), 
                               ncol = length(pred_terms)), 
                        row.names = outcome_vars)
colnames(results.F) <- pred_terms

# need columns for unadjusted and adjusted p-values
results.p <- cbind(results.F, results.F)  
# rename adjusted p-val columns 
colnames(results.p)[8:14] <- paste0(colnames(results.p)[8:14], "-adj")

# heritability dataframe 
h2 <- data.frame(phenotype = outcome_vars, 
                 heritability = vector(length = length(outcome_vars)))

# run LMM
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
    # Don't need to add data if this is last term or next term is Litter.Size^2, Sex:Diet
    if ((!is.na(pred_terms[j+1])) & (pred_terms[j+1]!='I(Litter_Size^2)') & (pred_terms[j+1]!='Sex:Diet')){
      k <- k+1
      moddata <- cbind(moddata, dat[,pred_vars[k]]) # add next col
      colnames(moddata)[k+3] <- pred_vars[k]
    }
    
    # if there are terms left to add, current large model becomes small model for 
    # the next comparison (using new dataset with NAs omitted for next variable as well)
    if (!is.na(pred_terms[j+1])){
      formula0 <- formula1 
      mod0 <- relmatLmer(formula0, data = na.omit(moddata), relmat = list(Access_ID = K))
    } else { # this is the last term, so calculate heritability from the full model 
      h2$heritability[i] <- VarProp(mod1)$prop[1]
    }
  }
}

#-----------------------------P-value correction-------------------------------#
# Rank p-values from smallest (most significant) to largest (least significant) 
# Assign rank value 1-19
# Calculate adjusted p-value: original p-value * (17/rank)

for (i in 1:length(pred_terms)){
  # determine ranks (keep originals so you can reorder by them)
  pdf <- data.frame(og.rank = seq(1,nrow(results.p)), pval = results.p[,i])
  pdf <- pdf[order(pdf$pval),] # rank by significance 
  pdf$rank <- seq(1, nrow(pdf))
  
  # calculate adjusted p-value 
  pdf$adjpval <- pdf$pval*(nrow(results.p)/pdf$rank) 
  
  # return to original order 
  pdf <- pdf[order(pdf$og.rank),]
  
  # replace pvals in results.p with adjusted pvals 
  results.p[,i+7] <- pdf$adjpval  
}

#--------------------------------Save results----------------------------------#
write.csv(results.p, file = "results/covariate-analysis.csv")
#write.csv(results.F, file = "results/anova-results-F.csv")
write.csv(h2, file = 'results/lmm-h2-estimates.csv')



#-------------------------Sex-stratified analysis------------------------------#
res <- data.frame(phenotype = outcome_vars,
                  diet_P_male = rep(NA, length(outcome_vars)),
                  diet_adjP_male = rep(NA, length(outcome_vars)),
                  diet_P_female = rep(NA, length(outcome_vars)),
                  diet_adjP_female = rep(NA, length(outcome_vars)))

#----------------------------------Male----------------------------------------#
datM <- dat[which(dat$Sex == 'M'),]

for (i in 1:length(outcome_vars)){
  pheno.name <- outcome_vars[i]
  
  ### ADDED DIET_START_WEEK 
  formula0 <- paste0(pheno.name, " ~ Study + Litter_Size + I(Litter_Size^2) + Diet_Start_Week + (1|Cohort) + (1|Access_ID)")
  mod0 <- relmatLmer(formula0, data = datM, relmat = list(Access_ID = K))
  
  ### ADDED DIET_START_WEEK
  formula1 <- paste0(pheno.name, " ~ Diet + Study + Litter_Size + I(Litter_Size^2) + Diet_Start_Week + (1|Cohort) + (1|Access_ID)")
  mod1 <- relmatLmer(formula1, data = datM, relmat = list(Access_ID = K))
  
  kr <- KRmodcomp(mod1, mod0)
  #Fval <- kr$test$stat[1]
  res[i, 'diet_P_male'] <- kr$test$p.value[1]
}

# P-value correction 
pdf <- data.frame(og.rank = seq(1,nrow(res)), pval = res$diet_P_male) 
pdf <- pdf[order(pdf$pval),] # rank by significance 
pdf$rank <- seq(1, nrow(pdf)) 
pdf$adjpval <- pdf$pval*(nrow(res)/pdf$rank) # calculate adjusted p-value
pdf <- pdf[order(pdf$og.rank),] # return to original order 
res[,'diet_adjP_male'] <- pdf$adjpval # add to res 


#---------------------------------Female---------------------------------------#
datF <- dat[which(dat$Sex == 'F'),]

for (i in 1:length(outcome_vars)){
  pheno.name <- outcome_vars[i]
  
  # no study predictor bc all females are from one study 
  ### ADDED DIET_START_WEEK
  formula0 <- paste0(pheno.name, " ~ Litter_Size + I(Litter_Size^2) + Diet_Start_Week + (1|Cohort) + (1|Access_ID)")
  mod0 <- relmatLmer(formula0, data = datF, relmat = list(Access_ID = K))
  
  ### ADDED DIET_START_WEEK
  formula1 <- paste0(pheno.name, " ~ Diet + Litter_Size + I(Litter_Size^2) + Diet_Start_Week + (1|Cohort) + (1|Access_ID)")
  mod1 <- relmatLmer(formula1, data = datF, relmat = list(Access_ID = K))
  
  kr <- KRmodcomp(mod1, mod0)
  #Fval <- kr$test$stat[1]
  res[i, 'diet_P_female'] <- kr$test$p.value[1]
}

# P-value correction 
pdf <- data.frame(og.rank = seq(1,nrow(res)), pval = res$diet_P_female) 
pdf <- pdf[order(pdf$pval),] # rank by significance 
pdf$rank <- seq(1, nrow(pdf)) 
pdf$adjpval <- pdf$pval*(nrow(res)/pdf$rank) # calculate adjusted p-value
pdf <- pdf[order(pdf$og.rank),] # return to original order 
res[,'diet_adjP_female'] <- pdf$adjpval # add to res 


#--------------------------------Save results----------------------------------#
write.csv(res, file = "results/diet-effects.csv")




