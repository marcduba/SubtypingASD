##################################
# Matching MRI & Behavioral Data #
##################################

library(dplyr)
library(stringr)
source('~/Documents/Psychologie/MA/Code/R/GeneralParameters.R')

setwd(params$WorkingDirectory)

# Loading ID lists
idASD <- read.delim('./Data/IDlists/id_ASD.txt', header = F)
colnames(idASD) <- 'EID'
idHC <- read.delim('./Data/IDlists/id_HC.txt', header = F)
colnames(idHC) <- 'EID'

IDcb <- read.delim('./Data/IDlists/ID_CB.txt', header = F) 
IDcb$Site <- 'CB'
IDcu <- read.delim('./Data/IDlists/ID_CU.txt', header = F)
IDcu$Site <- 'CU'
IDru <- read.delim('./Data/IDlists/ID_RU.txt', header = F)
IDru$Site <- 'RU'
IDsi <- read.delim('./Data/IDlists/ID_SI.txt', header = F)
IDsi$Site <- 'SI'

IDmri <- rbind(IDcb, IDcu, IDru, IDsi) # all IDs from Methlab-server

listIDmri <- lapply(IDmri$V1, IDsplitFunction) # creates list with all info from IDs

dfIDmri <- as.data.frame(t(sapply(listIDmri, '[', 1:3)))
dfIDmri <- cbind(IDmri$V1, dfIDmri, IDmri$Site)

colnames(dfIDmri) <- c('Full ID', 'EID', 'Protocol', 'Run', 'Site')

mriTotaln2 <- nrow(dfIDmri)

# Preference for Run 02
table(dfIDmri$Run)
dfRun01 <- dfIDmri %>%
  filter(Run == 1)
dfRun02 <- dfIDmri %>%
  filter(Run == 2)

dfRuns <- rbind(dfRun02, dfRun01)
dfRuns$Run

# Preference for Protocol
table(dfRuns$Protocol)
nHCP <- nrow(filter(dfRuns, Protocol == params$Protocols$HCP))
nVNav <- nrow(filter(dfRuns, Protocol == params$Protocols$VNav))
nVNavN <- nrow(filter(dfRuns, Protocol == params$Protocols$VNavN))
nOther <- nrow(filter(dfRuns, Protocol == params$Protocols$Other))

mriHCP <- dfRuns[dfRuns$Protocol == params$Protocols$HCP, ]
mriVNav <- dfRuns[dfRuns$Protocol == params$Protocols$VNav, ]
mriVNavN <- dfRuns[dfRuns$Protocol == params$Protocols$VNavN, ]
mriOther <- dfRuns[dfRuns$Protocol == params$Protocols$Other, ]

dfmri <- rbind(mriHCP, rbind(mriVNav, rbind(mriVNavN, mriOther)))
dfmri$Protocol

mri <- dfmri %>%
  filter(!duplicated(EID))
nTotalmri <- nrow(mri)

# Matching IDs from behavioral & MRI data
ASDavailable <- inner_join(mri, idASD, by = 'EID')
HCavailable <- inner_join(mri, idHC, by = 'EID')

table(ASDavailable$Site)
ASD_CB <- ASDavailable %>%
  filter(Site == 'CB')
ASD_CU <- ASDavailable %>%
  filter(Site == 'CU')
ASD_RU <- ASDavailable %>%
  filter(Site == 'RU')

table(HCavailable$Site)
HC_CB <- HCavailable %>%
  filter(Site == 'CB')
HC_CU <- HCavailable %>%
  filter(Site == 'CU')
HC_RU <- HCavailable %>%
  filter(Site == 'RU')

# saving full IDs for stats extraction
setwd('~/Documents/Psychologie/MA/Data/IDlists/')
write.table(ASD_CB$`Full ID`, file = 'ASD_CB.txt', quote = F, sep = '', row.names = F, col.names = F)
write.table(ASD_CU$`Full ID`, file = 'ASD_CU.txt', quote = F, sep = '', row.names = F, col.names = F)
write.table(ASD_RU$`Full ID`, file = 'ASD_RU.txt', quote = F, sep = '', row.names = F, col.names = F)

write.table(HC_CB$`Full ID`, file = 'HC_CB.txt', quote = F, sep = '', row.names = F, col.names = F)
write.table(HC_CU$`Full ID`, file = 'HC_CU.txt', quote = F, sep = '', row.names = F, col.names = F)
write.table(HC_RU$`Full ID`, file = 'HC_RU.txt', quote = F, sep = '', row.names = F, col.names = F)

# saving protocol for further analysis
protocol <- rbind(ASD_CB[, c('EID', 'Protocol')],
                       ASD_CU[, c('EID', 'Protocol')], 
                       ASD_RU[, c('EID', 'Protocol')], 
                       HC_CB[, c('EID', 'Protocol')], 
                       HC_CU[, c('EID', 'Protocol')], 
                       HC_RU[, c('EID', 'Protocol')])
write.table(protocol, file = 'Protocol.txt', quote = F, sep = '\t', row.names = F, col.names = F)
