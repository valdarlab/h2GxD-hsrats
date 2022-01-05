### This script transforms raw data based on phenotype_transformations.xlsx 
### and creates a new spreadsheet with transformed data 
library(readxl)
library(tidyverse)
source("utils.R")

metadata <- c("Study", "Parent_Pair", "Cohort", "Access_ID", "Sex", "Diet", 
              "Diet_Start_Week", "Litter_Size")
phenotypes <- c("Total_AUC", "HarvWeight", "RetroFat_norm", "FastGluc",
                "FastIns", "SWIM", "CLIMB", "FLOAT", "REST_EPISODE_COUNT_5",
                "MOVEMENT_EPISODE_COUNT_5", "VERTICAL_EPISODE_COUNT_5")

#---------------------------------Load data------------------------------------#
female_data <- read_csv("source_data/raw_data_female-G3package.csv",
                        col_types = list(Sex = col_character()))
colnames(female_data) <- gsub(" ", "_", colnames(female_data))

male_data <- read_csv("source_data/raw_data_male-G3package.csv")
colnames(male_data) <- gsub(" ", "_", colnames(male_data))

transform <- read_xlsx("source_data/phenotype_transformations.xlsx")
transform$Phenotype <- gsub(" ", "_", transform$Phenotype)

#--------------------------------Filter data-----------------------------------#
female_data <- female_data %>% dplyr::select(c(all_of(metadata), all_of(phenotypes)))
male_data <- male_data %>% dplyr::select(c(all_of(metadata), all_of(phenotypes)))

#------------------------------Transform data----------------------------------#
# loop over columns in male_data and female_data 
for (sex in c("M", "F")){
  data_tbl <- ifelse(sex == "M", "male_data", "female_data")
  data_tbl <- get(data_tbl)
  for (i in 1:length(phenotypes)){
    # get phenotype and transformation type 
    phenotype <- phenotypes[i]
    transform_type <- transform %>% 
      filter(Phenotype == phenotype, Sex == sex) %>% 
      dplyr::select(Transformation)
    
    # transform and replace data 
    if (transform_type == "nat log"){
      data_tbl[phenotype] <- log(data_tbl[[phenotype]]+1)
    } else if (transform_type == "inverse"){
      data_tbl[phenotype] <- (data_tbl[[phenotype]])^(-1)
    } else if (transform_type == "squared"){
      data_tbl[phenotype] <- (data_tbl[[phenotype]])^(2)
    } else if (transform_type == "sq root"){
      data_tbl[phenotype] <- sqrt(data_tbl[[phenotype]])
    } else if (transform_type == "raw"){
      # do nothing 
    }
  }
  if (sex == "M"){
    male_data <- data_tbl 
  } else { # sex == "F"
    female_data <- data_tbl
  }
}

#---------------------------Define other variables-----------------------------#
female_data$z <- ifelse(female_data$Diet == "HFD", 1, -1)
male_data$z <- ifelse(male_data$Diet == "HFD", 1, -1)

#--------------------------------Save to csv-----------------------------------#
ensure_directory("derived_data")
write_csv(male_data, file = "derived_data/transformed_data_male-G3package.csv")
write_csv(female_data, file = "derived_data/transformed_data_female-G3package.csv")

