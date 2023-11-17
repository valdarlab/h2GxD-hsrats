### This script transforms raw data based on pheno_transformations.csv
### and creates a new spreadsheet with transformed data 
library(tidyverse)
source("utils.R")

metadata <- c("Study", "Parent_Pair", "Cohort", "Access_ID", "Sex", "Diet", 
              "Diet_Start_Week", "Litter_Size")

tx <- read_csv("source_data/pheno_transformations.csv")
phenotypes <- tx$Phenotype

#------------------------------Load/filter data--------------------------------#
dat <- read_csv("source_data/raw_data_all-G3package.csv")
dat <- dat %>% dplyr::select(c(all_of(metadata), all_of(phenotypes)))

#----------------------------Summarize raw data--------------------------------#
data.sum <- data.frame(Sex = c(rep("Male", length(phenotypes)), rep("Female", length(phenotypes))),
                       Phenotype = rep(phenotypes, 2),
                       HFD.mean = rep(NA, length = length(phenotypes)*2),
                       HFD.sd = rep(NA, length = length(phenotypes)*2),
                       LFD.mean = rep(NA, length = length(phenotypes)*2),
                       LFD.sd = rep(NA, length = length(phenotypes)*2))

for (i in 1:length(phenotypes)){
  pheno.m.hfd <- dat %>% filter(Sex == 'M', Diet == "HFD") %>% select(all_of(phenotypes[i])) 
  pheno.m.hfd.mean <- mean(pheno.m.hfd[[1]], na.rm = TRUE)
  pheno.m.hfd.sd <- sd(pheno.m.hfd[[1]], na.rm = TRUE)
  
  pheno.m.lfd <- dat %>% filter(Sex == 'M', Diet == "LFD") %>% select(all_of(phenotypes[i])) 
  pheno.m.lfd.mean <- mean(pheno.m.lfd[[1]], na.rm = TRUE)
  pheno.m.lfd.sd <- sd(pheno.m.lfd[[1]], na.rm = TRUE)
  
  pheno.f.hfd <- dat %>% filter(Sex == 'F', Diet == "HFD") %>% select(all_of(phenotypes[i])) 
  pheno.f.hfd.mean <- mean(pheno.f.hfd[[1]], na.rm = TRUE)
  pheno.f.hfd.sd <- sd(pheno.f.hfd[[1]], na.rm = TRUE)
  
  pheno.f.lfd <- dat %>% filter(Sex == 'F', Diet == "LFD") %>% select(all_of(phenotypes[i])) 
  pheno.f.lfd.mean <- mean(pheno.f.lfd[[1]], na.rm = TRUE)
  pheno.f.lfd.sd <- sd(pheno.f.lfd[[1]], na.rm = TRUE)

  m.data <- c(pheno.m.hfd.mean, pheno.m.hfd.sd, pheno.m.lfd.mean, pheno.m.lfd.sd)
  f.data <- c(pheno.f.hfd.mean, pheno.f.hfd.sd, pheno.f.lfd.mean, pheno.f.lfd.sd)
  
  data.sum[which(data.sum$Sex == 'Male' & data.sum$Phenotype == phenotypes[i]),c(3:6)] <- m.data
  data.sum[which(data.sum$Sex == 'Female' & data.sum$Phenotype == phenotypes[i]),c(3:6)] <- f.data
}

ensure_directory("results")
write_csv(data.sum, file = "results/phenotype_summary.csv")

#------------------------------Transform data----------------------------------#
for (i in 1:length(phenotypes)){
  phenotype <- tx$Phenotype[i]
  transformation <- tx$Transformation[i]
  
  if (transformation == "natlog"){
    dat[phenotype] <- log(dat[phenotype]+1)
  } else if (transformation == "inv"){
    dat[phenotype] <- dat[phenotype]^(-1)
  } else if (transformation == "squared"){
    dat[phenotype] <- dat[phenotype]^2
  } else if (transformation == "sqrt"){
    dat[phenotype] <- sqrt(dat[phenotype])
  } else if (transformation == "cbrt"){
    dat[phenotype] <- dat[phenotype]^(1/3)
  } else if (transformation == "raw"){
    # do nothing 
  }
}

#---------------------------Define other variables-----------------------------#
dat$z <- ifelse(dat$Diet == "HFD", 1, -1)

#-----------------------------Save to csv by sex-------------------------------#
mdat <- dat %>% filter(Sex == "M")
fdat <- dat %>% filter(Sex == 'F')

ensure_directory("derived_data")
write_csv(dat, file = "derived_data/transformed_data_All-G3package.csv")
write_csv(mdat, file = "derived_data/transformed_data_Male-G3package.csv")
write_csv(fdat, file = "derived_data/transformed_data_Female-G3package.csv")

