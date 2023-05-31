########################  Creating a relatedness matrix using QTlRel kinship()  ##################### 
# Below is example code for how to create a relatedness matrix using QTlRel kinship() function. For #
# additional information, please see the package documentation at this address:                     #
#      https://cran.r-project.org/web/packages/QTLRel/QTLRel.pdf                                    #
#####################################################################################################
library(dplyr)
library(QTLRel)

# Read in the csv file for each generation.
for (i in c(0:34)){
  obj_name <- paste0('gen', i)
  assign(obj_name, read.delim(file = paste0('source_data/generations/generation', i, '.csv'), 
                              sep = ",", header = TRUE, stringsAsFactors = FALSE))
}

# load gen 35
gen35 <- read.delim(file = "source_data/gen35_MaleFemale_MetPilot_BehPilot_Thesis.csv", 
                           sep= ",", header = TRUE, stringsAsFactors = FALSE)

# create list with all generations (except last one) 
gens <- list(gen0, gen1, gen2, gen3, gen4, gen5, gen6, gen7, gen8, gen9, gen10,
             gen11, gen12, gen13, gen14, gen15, gen16, gen17, gen18, gen19, gen20,
             gen21, gen22, gen23, gen24, gen25, gen26, gen27, gen28, gen29, gen30,
             gen31, gen32, gen33, gen34, gen35)

# pick out the columns you want: kinship() requires an ID, Sex, Sire, Dam, and Generation.
colselect <- function(df){df[,c(1,3,7,8,11)]}
gens <- lapply(X = gens, FUN = colselect) 

# change column names 
changeNames <- function(df){names(df) <- c("id", "sex", "sire", "dam", "generation"); return(df)}
gens <- lapply(X = gens, FUN = changeNames) 

# change all generations from 0, 100, 200, 300, ... to 0, 1, 2, 3, ... 
for (i in seq_along(gens)){ 
  gen_num <- i-1
  gens[[i]]$generation <- gen_num
}

# create male/female-specific gen35 df  
gen35 <- gens[[36]]
gen35M <- gen35[gen35$sex == 'M',]
gen35F <- gen35[gen35$sex == 'F',]

# create male/female-specific df lists 
gensM <- append(gens[1:35], list(gen35M))
gensF <- append(gens[1:35], list(gen35F))

# bind all of the generations together to create a pedigree.
pedigree <- bind_rows(gens)
pedigreeM <- bind_rows(gensM)
pedigreeF <- bind_rows(gensF)

# save pedigree into a file so that you can start from the full pedigree instead 
# of from the individual generations. 
write.table(pedigree, file = 'derived_data/MaleFemale_PilotThesis_Pedigree.txt')
write.table(pedigreeM, file = 'derived_data/Male_PilotThesis_Pedigree.txt')
write.table(pedigreeF, file = 'derived_data/Female_PilotThesis_Pedigree.txt')

# create the relatedness matrix using kinship().
kinship_matrix <- as.data.frame(kinship(pedigree, pedigree$id))
kinship_matrixM <- as.data.frame(kinship(pedigreeM, pedigreeM$id))
kinship_matrixF <- as.data.frame(kinship(pedigreeF, pedigreeF$id))

# save full kinship matrices
# write.table(kinship_matrix, file = 'derived_data/full_relatedness_matrix_MaleFemale_PilotThesis.txt')
# write.table(kinship_matrixM, file = 'derived_data/full_relatedness_matrix_Male_PilotThesis.txt')
# write.table(kinship_matrixF, file = 'derived_data/full_relatedness_matrix_Female_PilotThesis.txt')

# since you created the pedigree that was fed to kinship(), you know the exact order 
# of all of the animals in the matrix. You can subset the relatedness matrix and 
# extract only the animals/generation that you are interested in. 
n = 512 # 513 total rats
kinship_subset <- kinship_matrix[(nrow(kinship_matrix)-n):nrow(kinship_matrix), 
                                 (ncol(kinship_matrix)-n):ncol(kinship_matrix)]
n = 324 # 325 male rats 
kinship_subsetM <- kinship_matrixM[(nrow(kinship_matrixM)-n):nrow(kinship_matrixM), 
                                   (ncol(kinship_matrixM)-n):ncol(kinship_matrixM)]
n = 187 # 188 female rats
kinship_subsetF <- kinship_matrixF[(nrow(kinship_matrixF)-n):nrow(kinship_matrixF), 
                                   (ncol(kinship_matrixF)-n):ncol(kinship_matrixF)]

# save subsetted kinship matrices 
write.table(kinship_subset*2, 'derived_data/relatedness_matrix_MaleFemale_PilotThesis.txt')
write.table(kinship_subsetM*2, 'derived_data/relatedness_matrix_Male_PilotThesis.txt')
write.table(kinship_subsetF*2, 'derived_data/relatedness_matrix_Female_PilotThesis.txt')

