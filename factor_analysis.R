# Factor analysis
library(tidyverse)
library(missMDA)
library(factoextra)

### DELETE THIS 
setwd('/Users/ellenrisemberg/Documents/ValdarFerris/GxD0531')

#---------------------------------Load Data------------------------------------#
dat <- read.csv('derived_data/transformed_data_All-G3package.csv')

phenos <- c('HarvWeight', 'RetroFat_norm', 'EpiFat_norm', 'OmentalFat_norm', 
            'Total_AUC', 'FastGluc', 'ClosedJunc', 'OpenJunc', 'SWIM',
            'CLIMB', 'FLOAT', 'REST_EPISODE_COUNT_5', 'MOVEMENT_EPISODE_COUNT_5',
            'VERTICAL_EPISODE_COUNT_5')

fadat <- dat %>% select(phenos)

#--------------------------Imputation/Iterative PCA----------------------------#
ncp <- estim_ncpPCA(fadat, scale = TRUE) 
comp_data <- imputePCA(fadat, method = "regularized", ncp = ncp$ncp, scale = TRUE)

# run iterative PCA on imputed dataset - determine number of factors to use in FA 
ipca <- prcomp(comp_data$completeObs, scale = TRUE)
#print(fviz_eig(ipca)) # print scree plot 
ggsave(fviz_eig(ipca), filename = 'results/plots/scree-plot.pdf')

#-------------------------------Factor Analysis--------------------------------#
fares <- factanal(x = comp_data$completeObs, factors = 4, rotation = "varimax")
write.csv(fares$loadings, file = 'results/factor-loadings.csv')

#------------------------------Factor Phenotypes-------------------------------# 
for (i in 1:4){
  load <- loadings(fares)[,i]
  cname <- paste0('Factor', i)
  dat[,cname] <- rowSums(sweep(fadat, MARGIN = 2, load, FUN = '*'), na.rm = TRUE)
}


#---------------------------------Save by Sex----------------------------------# 
mdat <- dat %>% filter(Sex == "M")
fdat <- dat %>% filter(Sex == 'F')

write_csv(dat, file = "derived_data/transformed_data_All-G3package.csv")
write_csv(mdat, file = "derived_data/transformed_data_Male-G3package.csv")
write_csv(fdat, file = "derived_data/transformed_data_Female-G3package.csv")



