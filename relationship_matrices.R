library(AGHmatrix) # to get kinship matrix
# as used in https://onunicornsandgenes.blog/2019/10/13/using-r-animal-model-with-simulated-data/
source("utils.R")

#--------------------------------Data setup------------------------------------#
males.ped <- read.delim(sep=" ", "source_data/MetBehThesisMalesPedigree.txt")
females.ped <- read.delim(sep=" ", "source_data/MetBehThesisFemalesPedigree.txt")                                                                     
male.dat <- read.csv("derived_data/transformed_data_male-G3package.csv")
female.dat <- read.csv("derived_data/transformed_data_female-G3package.csv")

#-------------------------Calculate kinship matrices---------------------------#
all.ped <- rbind(males.ped, females.ped)
# prepare ped file for AGHmatrix library
strict.ped <- all.ped[, c(1, 3, 4)] 
strict.ped$sire[is.na(strict.ped$sire)] <- 0
strict.ped$dam[is.na(strict.ped$dam)] <- 0       
strict.ped <- strict.ped[order(strict.ped$id), ]
strict.ped <- strict.ped[!duplicated(strict.ped$id), ]

# convert missing parents to NA
missing.sires <- setdiff(strict.ped$sire, c(strict.ped$id, 0))
missing.dams <- setdiff(strict.ped$dam, c(strict.ped$id, 0))
strict.ped$sire[strict.ped$sire %in% missing.sires] <- 0
strict.ped$dam[strict.ped$dam %in% missing.dams] <- 0
cat("Missing id for sires {", paste(missing.sires, collapse=", "), "}, converting to founders\n")
cat("Missing id for dams {", paste(missing.dams, collapse=", "), "}, converting to founders\n")

K <- Amatrix(strict.ped) # calc additive relationship matrix for entire pedigree
hist(diag(K)) # show diagonal elements

#-------------------------Make sex-specific matrices---------------------------#
ii <- rownames(K) %in% male.dat$Access_ID
K.male <- K[ii, ii]
# heatmap(K.male, main=paste0("Males, n=", sum(ii)))
# range(diag(K.male))

ii <- rownames(K) %in% female.dat$Access_ID
K.female <- K[ii, ii]
# heatmap(K.female, main=paste0("Females, n=", sum(ii)))
# range(diag(K.female))

#----------------------------Save kinship matrices-----------------------------#
ensure_directory("derived_data")
write.table(K.male, file="derived_data/relatedness_matrix_MetBehThesis_male.txt")
write.table(K.female, file="derived_data/relatedness_matrix_MetBehThesis_female.txt")

