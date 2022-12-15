# Hierarchical Clustering

library(dplyr)
library(factoextra)
library(NbClust)
library(clusterpval)
library(fastcluster)
source('~/Documents/Psychologie/MA/Code/R/GeneralParameters.R')

# Loading Data
setwd(params$PathBrainData)

AI_ASD_SA <- read.delim('./ASD/AI_ASD_SA.txt', header = T, sep = '\t', dec = '.')
AI_ASD_CT <- read.delim('./ASD/AI_ASD_CT.txt', header = T, sep = '\t', dec = '.')

colnames(AI_ASD_SA) <- colnames(AI_ASD_SA) %>%
  paste('SA', sep = '_')
colnames(AI_ASD_CT) <- colnames(AI_ASD_CT) %>%
  paste('CT', sep = '_')

dataAI <- cbind(AI_ASD_SA, AI_ASD_CT)

# Standardize data
scale(dataAI)

# Compute dissimilarity matrix
distance <- dist(dataAI, method = params$Distance)

# Linkage
clust <- hclust(d = distance, method = params$ClustMethod)
dend <- fviz_dend(clust, show_labels = F, lwd = 0.5)
print(dend)

# Determine optimal number of clusters
Nclust <- NbClust(data = dataAI, method = params$ClustMethod, index = 'alllong')
table(Nclust$Best.nc[1,]) # best number of clusters = 3

dendK <- fviz_dend(clust, show_labels = F, lwd = 0.5, k = 3, k_colors = params$paletteASD, main = '', ggtheme = theme_classic(base_family = 'sans'))
print(dendK)

# Cutting tree into groups
groups <- cutree(clust, k = 3)
table(groups)

data_clust <- dataAI
data_clust$clust <- groups

# Writing clustered data file
write.table(data_clust, file = './data_clustered.txt', sep = '\t', dec = '.', row.names = T, col.names = T)


# Computing Gao correction (https://www.lucylgao.com/clusterpval/)
X <- as.matrix(data_clust[, -149])

Gao1v2 <- test_hier_clusters_exact(X, link = 'ward.D', hcl = clust, K = 3, k1 = 1, k2 = 2)
Gao1v2$pval
Gao1v3 <- test_hier_clusters_exact(X, link = 'ward.D', hcl = clust, K = 3, k1 = 1, k2 = 3)
Gao1v3$pval
Gao2v3 <- test_hier_clusters_exact(X, link = 'ward.D', hcl = clust, K = 3, k1 = 2, k2 = 3)
Gao2v3$pval

# Computing Gao correction without ASD-III (n = 3)
dataNO3 <- data_clust %>%
  filter(clust != 3)
dataNO3 <- dataNO3 %>%
  select(-c('clust'))

distNO3 <- dist(dataNO3, method = params$Distance)
clustNO3 <- hclust(d = distNO3, method = params$ClustMethod)
NclustNO3 <- NbClust(data = dataNO3, method = params$ClustMethod, index = 'alllong')
fviz_dend(clustNO3, show_labels = F, lwd = 0.5, k = 2, k_colors = params$paletteASD)
groupsNO3 <- cutree(clustNO3, k = 2)
table(groupsNO3)

data_clustNO3 <- dataNO3
data_clustNO3$clust <- groupsNO3

X2 <- as.matrix(data_clustNO3[, -(75:149)])
test_hier_clusters_exact(X2, link = 'ward.D', hcl = clustNO3, K = 2, k1 = 1, k2 = 2)
X3 <- as.matrix(data_clustNO3[, -c(1:74, 149)])
test_hier_clusters_exact(X3, link = 'ward.D', hcl = clustNO3, K = 2, k1 = 1, k2 = 2)

set.seed(123)
test_complete_hier_clusters_approx(X2, K=2, k1=1, k2=2, ndraws=10000, hcl=clustNO3)
set.seed(123)
test_complete_hier_clusters_approx(X3, K=2, k1=1, k2=2, ndraws=10000, hcl=clustNO3)

