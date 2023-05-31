# Factor analysis
library(tidyverse)

### DELETE THIS 
setwd('/Users/ellenrisemberg/Documents/ValdarFerris/GxD0522')

#dat <- read.csv("Copy of All Data MaleFemale RawData merged for sex LMM.csv")
dat <- read.csv('source_data/raw_data_all-G3package.csv')
#dat <- read.csv('derived_data/transformed_data_All-G3package.csv')

fadat <- dat %>% select('HarvWeight', 'RetroFat_norm', 'EpiFat_norm', 'OmentalFat_norm', 
                        'Total_AUC', 'FastGluc', 'ClosedJunc', 'OpenJunc', 'SWIM',
                        'CLIMB', 'FLOAT', 'REST_EPISODE_COUNT_5', 'MOVEMENT_EPISODE_COUNT_5',
                        'VERTICAL_EPISODE_COUNT_5')

# fadat <- dat %>% select(23, 24, 25, 26, 21, 20, 27, 28, 34, 32, 33, 29, 30, 31) 
# 20:21, 23:34
# 20-21: Total AUC, FastGluc, no FastIns (bc missing data?)
# 23-34: HarvWeight, RetroFat, EpiFat, OmentalFat, ClosedJunc, OpenJunc, 
# REST, MOVEMENT, VERTICAL, SWIM, CLIMB, FLOAT

fadato <- na.omit(fadat)

# principal components analysis 
res <- prcomp(fadato, scale = TRUE)
screeplot(res)

# calculate total variance explained by each principal component
var_explained = res$sdev^2 / sum(res$sdev^2)

# create scree plot
qplot(c(1:14), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

# factor analysis 
fares <- factanal(x = fadato, factors = 4, rotation = "varimax")

loadings(fares)

print(fares, cutoff=0.00001)

write_csv(loadings(fares), "results/fa-results-unnorm.csv")
