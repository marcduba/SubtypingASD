# Combining Braindata & Behavioral data

library(dplyr)
library(plyr)
library(stringr)
source('~/Documents/Psychologie/MA/Code/R/GeneralParameters.R')

setwd(params$PathBrainData)
dataASD <- read.delim('data_clustered.txt', header = T, sep = '\t', dec = '.')
AI_HC_SA <- read.delim('./HC/AI_HC_SA.txt', header = T, sep = '\t', dec = '.')
colnames(AI_HC_SA) <- colnames(AI_HC_SA) %>%
  paste('SA', sep = '_')
AI_HC_CT <- read.delim('./HC/AI_HC_CT.txt', header = T, sep = '\t', dec = '.')
colnames(AI_HC_CT) <- colnames(AI_HC_CT) %>%
  paste('CT', sep = '_')
dataHC <- cbind(AI_HC_SA, AI_HC_CT)

# Saving final EID lists
IDasd <- rownames(dataASD)
IDasd <- IDasd %>%
  str_replace_all('sub-', '') %>%
  str_replace_all('_acq-HCP_T1w', '') %>%
  str_replace_all('_acq-HCP_run-02_T1w', '') %>%
  str_replace_all('_acq-HCP_run-01_T1w', '') %>%
  str_replace_all('_acq-VNav_T1w', '')

IDasd <- data.frame(IDasd)
write.table(IDasd, file = '~/Documents/Psychologie/MA/Data/IDlists/ID_ASD_final.txt', sep = '\t', quote = F, row.names = F, col.names = F)

IDhc <- rownames(dataHC)
IDhc <- IDhc %>%
  str_replace_all('sub-', '') %>%
  str_replace_all('_acq-HCP_T1w', '') %>%
  str_replace_all('_acq-HCP_run-02_T1w', '') %>%
  str_replace_all('_acq-HCP_run-01_T1w', '') %>%
  str_replace_all('_acq-VNav_T1w', '')

IDhc <- data.frame(IDhc)
write.table(IDhc, file = '~/Documents/Psychologie/MA/Data/IDlists/ID_HC_final.txt', sep = '\t', quote = F, row.names = F, col.names = F)

# Combining behavioral and asymmetry data
BD <- read.csv(params$BasicDataPath)

dataASD['diag'] <- 'ASD'
dataASD['EID'] <- IDasd
dataHC['diag'] <- 'HC'
dataHC['clust'] <- 0
dataHC['EID'] <- IDhc

braindata <- rbind(dataASD, dataHC)

allData <- inner_join(braindata, BD, by = 'EID')

# Writing data file
setwd(params$WorkingDirectory)
write.csv(allData, file = './Data/allData.csv', sep = '\t', dec = '.', row.names = T, col.names = T)
