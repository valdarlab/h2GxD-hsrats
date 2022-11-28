### This script transforms raw data based on pheno_transformations.csv
### and creates a new spreadsheet with transformed data 
library(tidyverse)
source("utils.R")

metadata <- c("Study", "Parent_Pair", "Cohort", "Access_ID", "Sex", "Diet", 
              "Diet_Start_Week", "Litter_Size")

tx <- read_csv("source_data/pheno_transformations.csv")
phenotypes <- tx$Phenotype

#---------------------------------Load data------------------------------------#
dat <- read_csv("source_data/raw_data_all-G3package.csv")

#--------------------------------Filter data-----------------------------------#
dat <- dat %>% dplyr::select(c(all_of(metadata), all_of(phenotypes)))
mdat <- dat %>% filter(Sex == "M")
fdat <- dat %>% filter(Sex == 'F')

#------------------------------Transform data----------------------------------#
for (i in 1:length(phenotypes)){
  phenotype <- tx$Phenotype[i]
  transformation <- tx$Transformation[i]
  
  if (transformation == "natlog"){
    dat[phenotype] <- log(dat[phenotype]+1)
    fdat[phenotype] <- log(fdat[phenotype]+1)
    mdat[phenotype] <- log(mdat[phenotype]+1)
  } else if (transformation == "inv"){
    dat[phenotype] <- dat[phenotype]^(-1)
    fdat[phenotype] <- fdat[phenotype]^(-1)
    mdat[phenotype] <- mdat[phenotype]^(-1)
  } else if (transformation == "squared"){
    dat[phenotype] <- dat[phenotype]^2
    fdat[phenotype] <- fdat[phenotype]^2
    mdat[phenotype] <- mdat[phenotype]^2
  } else if (transformation == "sqrt"){
    dat[phenotype] <- sqrt(dat[phenotype])
    fdat[phenotype] <- sqrt(fdat[phenotype])
    mdat[phenotype] <- sqrt(mdat[phenotype])
  } else if (transformation == "cbrt"){
    dat[phenotype] <- dat[phenotype]^(1/3)
    fdat[phenotype] <- fdat[phenotype]^(1/3)
    mdat[phenotype] <- mdat[phenotype]^(1/3)
  } else if (transformation == "raw"){
    # do nothing 
  }
}

#---------------------------Define other variables-----------------------------#
dat$z <- ifelse(dat$Diet == "HFD", 1, -1)
mdat$z <- ifelse(mdat$Diet == "HFD", 1, -1)
fdat$z <- ifelse(fdat$Diet == "HFD", 1, -1)

#--------------------------------Save to csv-----------------------------------#
ensure_directory("derived_data")
write_csv(dat, file = "derived_data/transformed_data_all-G3package.csv")
write_csv(mdat, file = "derived_data/transformed_data_male-G3package.csv")
write_csv(fdat, file = "derived_data/transformed_data_female-G3package.csv")

