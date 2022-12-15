# Selecting subjects based on behavioral data

library(dplyr)

source('~/Documents/Psychologie/MA/Code/R/GeneralParameters.R')

BD <- read.csv(params$BasicDataPath)

nTotalBehav <- nrow(BD)

# Filtering for ASD Diagnosis
bdASD <- BD %>% 
  filter(DX_01 == params$ASDdiag | 
           DX_02 == params$ASDdiag |
           DX_03 == params$ASDdiag |
           DX_04 == params$ASDdiag |
           DX_05 == params$ASDdiag |
           DX_06 == params$ASDdiag | 
           DX_07 == params$ASDdiag |
           DX_08 == params$ASDdiag |
           DX_09 == params$ASDdiag |
           DX_10 == params$ASDdiag)

nDiagASD <- nrow(bdASD)

# Filtering for behavioral measures
bdBehavASD <- bdASD %>% 
  filter(ASSQ.ASSQ_Total >= 0 & SRS.SRS_Total >= 0 & RBS.RBS_Total >= 0 & SCQ.SCQ_Total >= 0)

nBehavASD <- nrow(bdBehavASD)

# Filtering for right-handedness
bdHandednessASD <- bdBehavASD %>%
  filter(EHQ.EHQ_Total >= params$EHQCutoff)

nHandednessASD <- nrow(bdHandednessASD)

bdHandednessASD['Diag'] <- 'ASD'


# Filtering for control subjects
bdHC <- BD %>% 
  filter(DX_01 == params$HCdiag)

nDiagHC <- nrow(bdHC)

# Filtering for behavioral measures
bdBehavHC <- bdHC %>%
  filter(ASSQ.ASSQ_Total >= 0 & SRS.SRS_Total >= 0 & RBS.RBS_Total >= 0 & SCQ.SCQ_Total >= 0)

nBehavHC <- nrow(bdBehavHC)

# Filtering for right-handedness
bdHandednessHC <- bdBehavHC %>% 
  filter(EHQ.EHQ_Total >= params$EHQCutoff)

nHandednessHC <- nrow(bdHandednessHC)

bdHandednessHC['Diag'] <- 'HC'

# Combining & saving Dataset
allData <- rbind(bdHandednessASD, bdHandednessHC)

# Lists of IDs for ASD & HC subjects
id_ASD <- data.frame(allData$EID[which(allData$Diag == 'ASD')])
colnames(id_ASD) <- 'EID'
id_HC <- data.frame(allData$EID[which(allData$Diag == 'HC')])
colnames(id_HC) <- 'EID'

write.table(id_ASD, file = '/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/id_ASD.txt', row.names = F, col.names = F, quote = F)
write.table(id_HC, file = '/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/id_HC.txt', row.names = F, col.names = F, quote = F)

