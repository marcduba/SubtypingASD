# Clustercomparison behavioral variables

library(dplyr)
library(plyr)
library(tidyr)
library(stringr)
library(jtools)
library(Hmisc)
library(effsize)
library(ggplot2)
library(ggpubr)
library(rempsyc)
source('~/Documents/Psychologie/MA/Code/R/GeneralParameters.R')
load('~/Documents/Psychologie/MA/Code/R/SignificantRegions.Rda')

setwd(params$WorkingDirectory)
allData <- read.csv('./Data/allData.csv')

# Load site & protocol file
setwd('~/Documents/Psychologie/MA/Data/IDlists/')
protocol <- read.table('Protocol.txt')
colnames(protocol) <- c('EID', 'Protocol')

# combine both dataframes
allData <- inner_join(allData, protocol, by = 'EID')


## Descriptive Statistics
params$BehavMeasures
params$SRSsubscales

desc_stats_function <- function(measure) {
  Ns <- list(
    N_Ctrls = round(length(na.omit(allData[allData[, 'diag'] == 'HC', measure])), digits = 0),
    N_ASD = round(length(na.omit(allData[allData[, 'diag'] == 'ASD', measure])), digits = 0), 
    N_ASD1 = round(length(na.omit(allData[allData[, 'clust'] == 1, measure])), digits = 0),
    N_ASD2 = round(length(na.omit(allData[allData[, 'clust'] == 2, measure])), digits = 0),
    N_ASD3 = round(length(na.omit(allData[allData[, 'clust'] == 3, measure])), digits = 0)
  )
  means <- list(
    Total_Ctrls = format(round(mean(allData[allData[, 'diag'] == 'HC', measure], na.rm = T), digits = 2), nsmall = 2),
    Total_ASD = format(round(mean(allData[allData[, 'diag'] == 'ASD', measure], na.rm = T), digits = 2), nsmall = 2),
    ASD1 = format(round(mean(allData[allData[,'clust'] == 1, measure], na.rm = T), digits = 2), nsmall = 2),
    ASD2 = format(round(mean(allData[allData[,'clust'] == 2, measure], na.rm = T), digits = 2), nsmall = 2), 
    ASD3 = format(round(mean(allData[allData[,'clust'] == 3, measure], na.rm = T), digits = 2), nsmall = 2)
  )
  SDs <- list(
    Total_Ctrls = format(round(sd(allData[allData[, 'diag'] == 'HC', measure], na.rm = T), digits = 2), nsmall = 2),
    Total_ASD = format(round(sd(allData[allData[, 'diag'] == 'ASD', measure], na.rm = T), digits = 2), nsmall = 2), 
    ASD1 = format(round(sd(allData[allData[, 'clust'] == 1, measure], na.rm = T), digits = 2), nsmall = 2), 
    ASD2 = format(round(sd(allData[allData[, 'clust'] == 2, measure], na.rm = T), digits = 2), nsmall = 2), 
    ASD3 = format(round(sd(allData[allData[, 'clust'] == 3, measure], na.rm = T), digits = 2), nsmall = 2)
  )
  results <- list(N = Ns, Mean = means, SD = SDs)
  return(results)
}

statsASSQ <- desc_stats_function(params$BehavMeasures[1])
statsSRS <- desc_stats_function(params$BehavMeasures[2])
statsAWR <- desc_stats_function(params$SRSsubscales[1])
statsCOG <- desc_stats_function(params$SRSsubscales[2])
statsCOM <- desc_stats_function(params$SRSsubscales[3])
statsMOT <- desc_stats_function(params$SRSsubscales[5])
statsRRB <- desc_stats_function(params$SRSsubscales[6])
statsRBS <- desc_stats_function(params$BehavMeasures[3])
statsSCQ <- desc_stats_function(params$BehavMeasures[4])
statsAge <- desc_stats_function(params$BehavMeasures[5])

# Gender 
N_male_total <- table(allData$Basic_Demos.Sex)[1]
N_fem_total <- table(allData$Basic_Demos.Sex)[2]
N_maleHC <- table(allData$Basic_Demos.Sex[which(allData$diag == 'HC')])[1]
N_maleASD <- table(allData$Basic_Demos.Sex[which(allData$diag == 'ASD')])[1]
N_maleASD1 <- table(allData$Basic_Demos.Sex[which(allData$clust == 1)])[1]
N_maleASD2 <- table(allData$Basic_Demos.Sex[which(allData$clust == 2)])[1]
N_maleASD3 <- table(allData$Basic_Demos.Sex[which(allData$clust == 3)])[1]

n_sexHC <- paste(table(allData$Basic_Demos.Sex[which(allData$diag == 'HC')])[1], table(allData$Basic_Demos.Sex[which(allData$diag == 'HC')])[2], sep = '/')
n_sexASD <- paste(table(allData$Basic_Demos.Sex[which(allData$diag == 'ASD')])[1], table(allData$Basic_Demos.Sex[which(allData$diag == 'ASD')])[2], sep = '/')
n_sexASD1 <- paste(table(allData$Basic_Demos.Sex[which(allData$clust == 1)])[1], table(allData$Basic_Demos.Sex[which(allData$clust == 1)])[2], sep = '/')
n_sexASD2 <- paste(table(allData$Basic_Demos.Sex[which(allData$clust == 2)])[1], table(allData$Basic_Demos.Sex[which(allData$clust == 2)])[2], sep = '/')
n_sexASD3 <- paste(table(allData$Basic_Demos.Sex[which(allData$clust == 3)])[1], table(allData$Basic_Demos.Sex[which(allData$clust == 3)])[2], sep = '/')

statsTable <- data.frame(Measure = c('', 'Sex', 'Age', 'ASSQ', 'SRS', ' AWR', ' COG', ' COM', ' MOT', ' RRB', 'RBS', 'SCQ'), 
                         N.Ctrls = c('N (m/f)', n_sexHC, statsAge$N$N_Ctrls, statsASSQ$N$N_Ctrls, statsSRS$N$N_Ctrls, statsAWR$N$N_Ctrls, statsCOG$N$N_Ctrls, statsCOM$N$N_Ctrls, statsMOT$N$N_Ctrls, statsRRB$N$N_Ctrls, statsRBS$N$N_Ctrls, statsSCQ$N$N_Ctrls), 
                         Mean.Ctrls = c('M', '-', statsAge$Mean$Total_Ctrls, statsASSQ$Mean$Total_Ctrls, statsSRS$Mean$Total_Ctrls, statsAWR$Mean$Total_Ctrls, statsCOG$Mean$Total_Ctrls, statsCOM$Mean$Total_Ctrls, statsMOT$Mean$Total_Ctrls, statsRRB$Mean$Total_Ctrls, statsRBS$Mean$Total_Ctrls, statsSCQ$Mean$Total_Ctrls),
                         SD.HC = c('SD', '-', statsAge$SD$Total_Ctrls, statsASSQ$SD$Total_Ctrls, statsSRS$SD$Total_Ctrls, statsAWR$SD$Total_Ctrls, statsCOG$SD$Total_Ctrls, statsCOM$SD$Total_Ctrls, statsMOT$SD$Total_Ctrls, statsRRB$SD$Total_Ctrls, statsRBS$SD$Total_Ctrls, statsSCQ$SD$Total_Ctrls),
                         N.ASD = c('N (m/f)', n_sexASD, statsAge$N$N_ASD, statsASSQ$N$N_ASD, statsSRS$N$N_ASD, statsAWR$N$N_ASD, statsCOG$N$N_ASD, statsCOM$N$N_ASD, statsMOT$N$N_ASD, statsRRB$N$N_ASD, statsRBS$N$N_ASD, statsSCQ$N$N_ASD),
                         Mean.ASD = c('M', '-', statsAge$Mean$Total_ASD, statsASSQ$Mean$Total_ASD, statsSRS$Mean$Total_ASD, statsAWR$Mean$Total_ASD, statsCOG$Mean$Total_ASD, statsCOM$Mean$Total_ASD, statsMOT$Mean$Total_ASD, statsRRB$Mean$Total_ASD, statsRBS$Mean$Total_ASD, statsSCQ$Mean$Total_ASD),
                         SD.ASD = c('SD', '-', statsAge$SD$Total_ASD, statsASSQ$SD$Total_ASD, statsSRS$SD$Total_ASD, statsAWR$SD$Total_ASD, statsCOG$SD$Total_ASD, statsCOM$SD$Total_ASD, statsMOT$SD$Total_ASD, statsRRB$SD$Total_ASD, statsRBS$SD$Total_ASD, statsSCQ$SD$Total_ASD),
                         N.ASD1 = c('N (m/f)', n_sexASD1, statsAge$N$N_ASD1, statsASSQ$N$N_ASD1, statsSRS$N$N_ASD1, statsAWR$N$N_ASD1, statsCOG$N$N_ASD1, statsCOM$N$N_ASD1, statsMOT$N$N_ASD1, statsRRB$N$N_ASD1, statsRBS$N$N_ASD1, statsSCQ$N$N_ASD1),
                         Mean.ASD1 = c('M', '-', statsAge$Mean$ASD1, statsASSQ$Mean$ASD1, statsSRS$Mean$ASD1, statsAWR$Mean$ASD1, statsCOG$Mean$ASD1, statsCOM$Mean$ASD1, statsMOT$Mean$ASD1, statsRRB$Mean$ASD1, statsRBS$Mean$ASD1, statsSCQ$Mean$ASD1),
                         SD.ASD1 = c('SD', '-', statsAge$SD$ASD1, statsASSQ$SD$ASD1, statsSRS$SD$ASD1, statsAWR$SD$ASD1, statsCOG$SD$ASD1, statsCOM$SD$ASD1, statsMOT$SD$ASD1, statsRRB$SD$ASD1, statsRBS$SD$ASD1, statsSCQ$SD$ASD1),
                         N.ASD2 = c('N (m/f)', n_sexASD2, statsAge$N$N_ASD2, statsASSQ$N$N_ASD2, statsSRS$N$N_ASD2, statsAWR$N$N_ASD2, statsCOG$N$N_ASD2, statsCOM$N$N_ASD2, statsMOT$N$N_ASD2, statsRRB$N$N_ASD2,  statsRBS$N$N_ASD2, statsSCQ$N$N_ASD2), 
                         Mean.ASD2 = c('M', '-', statsAge$Mean$ASD2, statsASSQ$Mean$ASD2, statsSRS$Mean$ASD2, statsAWR$Mean$ASD2, statsCOG$Mean$ASD2, statsCOM$Mean$ASD2, statsMOT$Mean$ASD2, statsRRB$Mean$ASD2, statsRBS$Mean$ASD2, statsSCQ$Mean$ASD2), 
                         SD.ASD2 = c('SD', '-', statsAge$SD$ASD2, statsASSQ$SD$ASD2, statsSRS$SD$ASD2, statsAWR$SD$ASD2, statsCOG$SD$ASD2, statsCOM$SD$ASD2, statsMOT$SD$ASD2, statsRRB$SD$ASD2, statsRBS$SD$ASD2, statsSCQ$SD$ASD2), 
                         N.ASD3 = c('N (m/f)', n_sexASD3, statsAge$N$N_ASD3, statsASSQ$N$N_ASD3, statsSRS$N$N_ASD3, statsAWR$N$N_ASD3, statsCOG$N$N_ASD3, statsCOM$N$N_ASD3, statsMOT$N$N_ASD3, statsRRB$N$N_ASD3, statsRBS$N$N_ASD3, statsSCQ$N$N_ASD3), 
                         Mean.ASD3 = c('M', '-', statsAge$Mean$ASD3, statsASSQ$Mean$ASD3, statsSRS$Mean$ASD3, statsAWR$Mean$ASD3, statsCOG$Mean$ASD3, statsCOM$Mean$ASD3, statsMOT$Mean$ASD3, statsRRB$Mean$ASD3, statsRBS$Mean$ASD3, statsSCQ$Mean$ASD3), 
                         SD.ASD3 = c('SD', '-', statsAge$SD$ASD3, statsASSQ$SD$ASD3, statsSRS$SD$ASD3, statsAWR$SD$ASD3, statsCOG$SD$ASD3, statsCOM$SD$ASD3, statsMOT$SD$ASD3, statsRRB$SD$ASD3, statsRBS$SD$ASD3, statsSCQ$SD$ASD3))

colnames(statsTable) <- c('Measure', '', 'Controls', '', '', 'Full ASD', '', '', 'ASD-I', '', '', 'ASD-II', '', '', 'ASD-III', '')
write.csv(statsTable, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/descriptive_statistics.csv')
write.csv(statsTable[, 1:7], file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/descriptive_statistics_ASDvC.csv')

write.table(statsTable[c(1, 3:12),], file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/descriptive_statisticsCA12.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')

# Create df with ASD1 & ASD2
dfASD <- allData %>%
  filter(clust == 1 | clust == 2) %>%
  select('clust', 'Site', 'Protocol', params$BehavMeasures, params$SRSsubscales[c(1:3, 5, 6)])

dfASD$group <- as.factor(ifelse(dfASD$clust == 1, 'ASD-I', ifelse(dfASD$clust == 2, 'ASD-II', 'ASD-III')))
dfASD$Basic_Demos.Sex <- as.factor(ifelse(dfASD$Basic_Demos.Sex == 0, 'M', 'F'))
dfASD$Site <- as.factor(dfASD$Site)
dfASD$Protocol <- as.factor(dfASD$Protocol)
colnames(dfASD) <- c('clust', 'Site', 'Protocol', 'ASSQ', 'SRS', 'RBS', 'SCQ', 'Age', 'Sex', 'AWR', 'COG', 'COM', 'MOT', 'RRB', 'group')

# Create df with Controls & ASD1 & ASD2
dfASDC <- allData %>%
  filter(clust == 1 | clust == 2 | clust == 0) %>%
  select('diag', 'clust', 'Site', 'Protocol', params$BehavMeasures, params$SRSsubscales[c(1:3, 5, 6)])

dfASDC$group <- as.factor(ifelse(dfASDC$clust == 1, 'ASD-I', ifelse(dfASDC$clust == 2, 'ASD-II', 'Controls')))
dfASDC$diag <- as.factor(dfASDC$diag)
dfASDC$Basic_Demos.Sex <- as.factor(ifelse(dfASDC$Basic_Demos.Sex == 0, 'M', 'F'))
dfASDC$Site <- as.factor(dfASDC$Site)
dfASDC$Protocol <- as.factor(dfASDC$Protocol)
colnames(dfASDC) <- c('diag', 'clust', 'Site', 'Protocol', 'ASSQ', 'SRS', 'RBS', 'SCQ', 'Age', 'Sex', 'AWR', 'COG', 'COM', 'MOT', 'RRB', 'group')

# Create df wit ASD1 & ASD2 & ASD3
dfASD3 <- allData %>%
  filter(clust == 1 | clust == 2 | clust == 3 | clust == 0) %>%
  select('clust', 'Site', 'Protocol', params$BehavMeasures, params$SRSsubscales[c(1:3, 5, 6)])

dfASD3$group <- as.factor(ifelse(dfASD3$clust == 1, 'ASD-I', ifelse(dfASD3$clust == 2, 'ASD-II', ifelse(dfASD3$clust == 3, 'ASD-III', 'Controls'))))
dfASD3$Basic_Demos.Sex <- as.factor(ifelse(dfASD3$Basic_Demos.Sex == 0, 'M', 'F'))
dfASD3$Site <- as.factor(dfASD3$Site)
dfASD3$Protocol <- as.factor(dfASD3$Protocol)
colnames(dfASD3) <- c('clust', 'Site', 'Protocol', 'ASSQ', 'SRS', 'RBS', 'SCQ', 'Age', 'Sex', 'AWR', 'COG', 'COM', 'MOT', 'RRB', 'group')

# Plots
plot_function <- function(df, measure) {
  density.plot <- ggplot(df, aes(x = df[, measure], col = df[, 'group'], group = df[, 'group'])) + 
    geom_density() + 
    theme_apa() + 
    labs(x = measure, y = 'Density') +
    theme(legend.position = 'none', 
          text = element_text(family = 'sans'), 
          axis.text.x = element_text(margin = margin(t=5)), 
          axis.title.x = element_text(margin = margin(t=10))) +
    scale_color_manual(values = c(params$paletteASD[1], params$paletteASD[3], params$paletteASD[2], 'black'))
  
  violin.plot <- ggplot(df, aes(x = df[, 'group'], y = df[, measure], col = df[, 'group'])) + 
    geom_jitter() +
    geom_violin(trim = F, width = 0.8, alpha = 0.5) + 
    geom_boxplot(width = 0.3, color = 'black') + 
    theme_apa() + 
    labs(x = 'Group', y = measure) + 
    theme(legend.position = 'none', 
          text = element_text(family = 'sans'), 
          axis.text.x = element_text(margin = margin(t=5)), 
          axis.title.x = element_text(margin = margin(t=10))) + 
    scale_color_manual(values = c(params$paletteASD[1], params$paletteASD[3], params$paletteASD[2], 'grey'))
  
  plot.list <- list(violin.plot, density.plot)
  plot <- ggarrange(plotlist = plot.list, align = 'h')
  results <- list(violin.plot, density.plot, plot)
  
  return(results)
}

plotASSQ.ASD <- plot_function(df = dfASD, measure = 'ASSQ')
plotASSQ.Ctr <- plot_function(df = dfASDC, measure = 'ASSQ')
plotASSQ.ASD3 <- plot_function(df = dfASD3, measure = 'ASSQ')

plotSRS.ASD <- plot_function(df = dfASD, measure = 'SRS')
plotSRS.Ctr <- plot_function(df = dfASDC, measure = 'SRS')
plotSRS.ASD3 <- plot_function(df = dfASD3, measure = 'SRS')

plotAWR.ASD <- plot_function(df = dfASD, measure = 'AWR')
plotAWR.Ctr <- plot_function(df = dfASD, measure = 'AWR')

plotCOG.ASD <- plot_function(df = dfASD, measure = 'COG')
plotCOG.Ctr <- plot_function(df = dfASD, measure = 'COG')

plotCOM.ASD <- plot_function(df = dfASD, measure = 'COM')
plotCOM.Ctr <- plot_function(df = dfASD, measure = 'COM')

plotMOT.ASD <- plot_function(df = dfASD, measure = 'MOT')
plotMOT.Ctr <- plot_function(df = dfASD, measure = 'MOT')

plotRRB.ASD <- plot_function(df = dfASD, measure = 'RRB')
plotRRB.Ctr <- plot_function(df = dfASD, measure = 'RRB')

plotRBS.ASD <- plot_function(df = dfASD, measure = 'RBS')
plotRBS.Ctr <- plot_function(df = dfASDC, measure = 'RBS')
plotRBS.ASD3 <- plot_function(df = dfASD3, measure = 'RBS')

plotSCQ.ASD <- plot_function(df = dfASD, measure = 'SCQ')
plotSCQ.Ctr <- plot_function(df = dfASDC, measure = 'SCQ')
plotSCQ.ASD3 <- plot_function(df = dfASD3, measure = 'SCQ')

plotAGE.ASD <- plot_function(df = dfASD, measure = 'Age')
plotAGE.Ctr <- plot_function(df = dfASDC, measure = 'Age')
plotAGE.ASD3 <- plot_function(df = dfASD3, measure = 'Age')

violinASSQ <- plotASSQ.Ctr[[1]]
violinSRS <- plotSRS.Ctr[[1]]
violinAWR  <- plotAWR.Ctr[[1]]
violinCOG <- plotCOG.Ctr[[1]]
violinCOM <- plotCOM.Ctr[[1]]
violinMOT <- plotMOT.Ctr[[1]]
violinRRB <- plotRRB.Ctr[[1]]
violinRBS <- plotRBS.Ctr[[1]]
violinSCQ <- plotSCQ.Ctr[[1]]

# all total scores / ASD1, ASD2, Controls
plotlistviolins <- list(violinASSQ, violinSRS, violinRBS, violinSCQ)
violinsCtr <- ggarrange(ncol = 2, nrow = 2, plotlist = plotlistviolins)
ggsave(violinsCtr, filename = '~/Documents/Psychologie/MA/Results/Behavioral/Figures/violins.png', width = 1200, height = 1000, units = 'px', scale = 3)

# all total scores and SRS subscales / ASD1, ASD2
plotlistASD <- list(plotASSQ.ASD[[3]], plotSRS.ASD[[3]], plotRBS.ASD[[3]], plotSCQ.ASD[[3]])
plotlistSRS <- list(plotAWR.ASD[[3]], plotCOG.ASD[[3]], plotCOM.ASD[[3]], plotMOT.ASD[[3]], plotRRB.ASD[[3]])

plotASD <- ggarrange(ncol = 2, nrow = 2, plotlist = plotlistASD)
plotSRS <- ggarrange(ncol = 2, nrow = 3, plotlist = plotlistSRS)
ggsave(plotSRS, filename = '~/Documents/Psychologie/MA/Results/Behavioral/Figures/SRSsubscales.png', width = 1200, height = 1200, units = 'px', scale = 3)

# all total scores / ASD1, ASD2, Controls
plotlistCtr <- list(plotASSQ.Ctr[[3]], plotSRS.Ctr[[3]], plotRBS.Ctr[[3]], plotSCQ.Ctr[[3]], plotAGE.Ctr[[3]])
plotCtrl <- ggarrange(ncol = 2, nrow = 3, plotlist = plotlistCtr)

# all total scores / ASD1, ASD2, ASD3, Controls
plotlistASD3 <- list(plotASSQ.ASD3[[3]], plotSRS.ASD3[[3]], plotRBS.ASD3[[3]], plotSCQ.ASD3[[3]])
plotASD3 <- ggarrange(ncol = 2, nrow = 2, plotlist = plotlistASD3)

# Distributions of Sex, Protocol and Site
n_protHC <- paste(table(dfASD3$Protocol[which(dfASD3$group == 'Controls')])[1], table(dfASD3$Protocol[which(dfASD3$group == 'Controls')])[2], sep = '/')
n_protASD <- paste(table(dfASD3$Protocol[which(dfASD3$group != 'Controls')])[1], table(dfASD3$Protocol[which(dfASD3$group != 'Controls')])[2], sep = '/')
n_protASD1 <- paste(table(dfASD3$Protocol[which(dfASD3$clust == 1)])[1], table(dfASD3$Protocol[which(dfASD3$clust == 1)])[2], sep = '/')
n_protASD2 <- paste(table(dfASD3$Protocol[which(dfASD3$clust == 2)])[1], table(dfASD3$Protocol[which(dfASD3$clust == 2)])[2], sep = '/')
n_protASD3 <- paste(table(dfASD3$Protocol[which(dfASD3$clust == 3)])[1], table(dfASD3$Protocol[which(dfASD3$clust == 3)])[2], sep = '/')

n_siteHC <- list(CBIC = table(dfASD3$Site[which(dfASD3$group == 'Controls')])[1], RUBIC = table(dfASD3$Site[which(dfASD3$group == 'Controls')])[2], CUNY = table(dfASD3$Site[which(dfASD3$group == 'Controls')])[3])
n_siteASD <- list(CBIC = table(dfASD3$Site[which(dfASD3$group != 'Controls')])[1], RUBIC = table(dfASD3$Site[which(dfASD3$group != 'Controls')])[2], CUNY = table(dfASD3$Site[which(dfASD3$group != 'Controls')])[3])
n_siteASD1 <- list(CBIC = table(dfASD3$Site[which(dfASD3$clust == 1)])[1], RUBIC = table(dfASD3$Site[which(dfASD3$clust == 1)])[2], CUNY = table(dfASD3$Site[which(dfASD3$clust == 1)])[3])
n_siteASD2 <- list(CBIC = table(dfASD3$Site[which(dfASD3$clust == 2)])[1], RUBIC = table(dfASD3$Site[which(dfASD3$clust == 2)])[2], CUNY = table(dfASD3$Site[which(dfASD3$clust == 2)])[3])
n_siteASD3 <- list(CBIC = table(dfASD3$Site[which(dfASD3$clust == 3)])[1], RUBIC = table(dfASD3$Site[which(dfASD3$clust == 3)])[2], CUNY = table(dfASD3$Site[which(dfASD3$clust == 3)])[3])

demos_table <- data.frame(Group = c('Controls', 'Full ASD', 'ASD-I', 'ASD-II', 'ASD-III'), 
                             'Sex (M/F)' = c(n_sexHC, n_sexASD, n_sexASD1, n_sexASD2, n_sexASD3), 
                             'Protocol (HCP/VNav)' = c(n_protHC, n_protASD, n_protASD1, n_protASD2, n_protASD3), 
                             'CBIC' = c(n_siteHC$CBIC, n_siteASD$CBIC, n_siteASD1$CBIC, n_siteASD2$CBIC, n_siteASD3$CBIC),
                             'RUBIC' = c(n_siteHC$RUBIC, n_siteASD$RUBIC, n_siteASD1$RUBIC, n_siteASD2$RUBIC, n_siteASD3$RUBIC),
                             'CUNY' = c(n_siteHC$CUNY, n_siteASD$CUNY, n_siteASD1$CUNY, n_siteASD2$CUNY, n_siteASD3$CUNY))

write.table(demos_table, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/demos.txt', row.names = F, col.names = T, quote = F, sep = '\t')

# Chi-square test for gender and site
gender1v2 <- chisq.test(table(dfASD$group, dfASD$Sex))
genderAvC <- chisq.test(table(dfASDC$diag, dfASDC$Sex))
print(c(gender1v2$statistic, gender1v2$parameter, gender1v2$p.value))
print(c(genderAvC$statistic, genderAvC$parameter, genderAvC$p.value))

site1v2 <- chisq.test(table(dfASD$group, dfASD$Site))
siteAvC <- chisq.test(table(dfASDC$diag, dfASDC$Site))
print(c(site1v2$statistic, site1v2$parameter, site1v2$p.value))
print(c(siteAvC$statistic, siteAvC$parameter, siteAvC$p.value))

# Fishers exact test for protocol
protocol1v2 <- fisher.test(table(dfASD$group, dfASD$Protocol))
protocolAvC <- fisher.test(table(dfASDC$diag, dfASDC$Protocol))
print(c(protocol1v2$estimate, protocol1v2$p.value, protocol1v2$alternative))
print(c(protocolAvC$estimate, protocolAvC$p.value, protocolAvC$alternative))


# Correlations between behavioral measures
corrmatrix_function_behav <- function(df1, grp) {
  corrmat <- rcorr(as.matrix(df1[df1['diag'] == grp, c(params$BehavMeasures[c(6,5,1,3,4,2)], params$SRSsubscales[c(1, 2, 3, 5, 6)])]) , type = 'pearson')
  corrmatR <- round(corrmat$r, digits = 2)
  corrmatP <- round(corrmat$P, digits = 5)
  
  # create significance stars
  sig.stars <- ifelse(corrmatP < .0001, "***", ifelse(corrmatP < .001, "*** ", ifelse(corrmatP < .01, "**  ", ifelse(corrmatP < .05, "*   ", "    "))))
  
  # display correlations with stars
  corrmat.sig <- matrix(paste(corrmatR, sig.stars, sep = ''), ncol = ncol(corrmatR))
  diag(corrmat.sig) <- paste(diag(corrmatR), '', sep = '')
  rownames(corrmat.sig) <- paste(c(1:11), colnames(corrmatR), sep = '. ')
  colnames(corrmat.sig) <- paste(colnames(corrmatR), '', sep = '')
  
  # remove upper part of matrix
  Rmat <- as.matrix(corrmat.sig)
  Rmat[upper.tri(Rmat, diag = T)] <- ''
  Rmat <- as.data.frame(Rmat)
  
  # remove empty column
  Rmat <- cbind(Rmat[1:length(Rmat)-1])
  
  matrices <- list(corrTable = Rmat, Hmap = corrmatR)
  return(matrices)
}

# only ASD subjects
cmatASD <- corrmatrix_function_behav(allData, 'ASD')
colnames(cmatASD$corrTable) <- paste(c(1:10), '', sep = '.')

# ASD and Control subjects
cmatC <- corrmatrix_function_behav(allData, c('ASD', 'HC'))
colnames(cmatC$corrTable) <- paste(c(1:10), '', sep = '.')

write.table(cmatASD$corrTable, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/corrTableASD.txt', row.names = T, col.names = T, quote = F, dec = '.', sep = '\t')
write.table(cmatC$corrTable, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/corrTableCtrl.txt', row.names = T, col.names = T, quote = F, dec = '.', sep = '\t')

# Checking distribution of behavioral variables
sm.density.compare(allData$ASSQ.ASSQ_Total, allData$clust, xlab = 'ASSQ Total', col = c('black', params$paletteASD[1], params$paletteASD[3], params$paletteASD[2]), lty = c('solid', 'solid', 'solid', 'solid'))
qqnorm(allData$ASSQ.ASSQ_Total[which(allData$clust == 1) ])
qqline(allData$ASSQ.ASSQ_Total[which(allData$clust == 1)], col = 2)
qqnorm(allData$ASSQ.ASSQ_Total[which(allData$clust == 2) ])
qqline(allData$ASSQ.ASSQ_Total[which(allData$clust == 2)], col = 2)

sm.density.compare(allData$SRS.SRS_Total, allData$clust, xlab = 'SRS Total', col = c('black', params$paletteASD[1], params$paletteASD[3], params$paletteASD[2]), lty = c('solid', 'solid', 'solid', 'solid'))
qqnorm(allData$SRS.SRS_Total[which(allData$clust == 1)])
qqline(allData$SRS.SRS_Total[which(allData$clust == 1)], col = 2)
qqnorm(allData$SRS.SRS_Total[which(allData$clust == 2)])
qqline(allData$SRS.SRS_Total[which(allData$clust == 2)], col = 2)

sm.density.compare(allData$RBS.RBS_Total, allData$clust, xlab = 'RBS Total', col = c('black', params$paletteASD[1], params$paletteASD[3], params$paletteASD[2]), lty = c('solid', 'solid', 'solid', 'solid'))
qqnorm(allData$RBS.RBS_Total[which(allData$clust == 1)])
qqline(allData$RBS.RBS_Total[which(allData$clust == 1)], col = 2)
qqnorm(allData$RBS.RBS_Total[which(allData$clust == 2)])
qqline(allData$RBS.RBS_Total[which(allData$clust == 2)], col = 2)

sm.density.compare(allData$SCQ.SCQ_Total, allData$clust, xlab = 'SCQ Total', col = c('black', params$paletteASD[1], params$paletteASD[3], params$paletteASD[2]), lty = c('solid', 'solid', 'solid', 'solid'))
qqnorm(allData$SCQ.SCQ_Total[which(allData$clust == 1)])
qqline(allData$SCQ.SCQ_Total[which(allData$clust == 1)], col = 2)
qqnorm(allData$SCQ.SCQ_Total[which(allData$clust == 2)])
qqline(allData$SCQ.SCQ_Total[which(allData$clust == 2)], col = 2)

# Nyholt correction
nyholt_behav <- getMeff(dfASD[, c('ASSQ', 'SRS', 'RBS', 'SCQ', 'Age', 'AWR', 'COG', 'COM', 'MOT', 'RRB')])

### Hypothesis 2: Predict behavioral measures based on brain regions & cluster

# Creating data frames with significant FSbrainregionnames
significantRegions$CT$ASD1v2
significantRegions$SA$ASD1v2

regionnames_CT <- significantRegions$CT$ASD1v2
regionnames_SA <- significantRegions$SA$ASD1v2
regionnames_CT <- paste0(regionnames_CT, '_CT', sep = '')
regionnames_SA <- paste0(regionnames_SA, '_SA', sep = '')
regionnames_CT <- str_replace_all(regionnames_CT, '[-]', '.')
regionnames_SA <- str_replace_all(regionnames_SA, '[-]', '.')

glmData <- allData %>%
  filter(clust == 1 | clust == 2) %>%
  select(clust, params$BehavMeasures[6:1], all_of(regionnames_CT), all_of(regionnames_SA))
colnames(glmData) <- colnames(glmData) %>%
  str_replace_all('Basic_Demos.', '') %>%
  str_replace_all('.SCQ_Total', '') %>%
  str_replace_all('.RBS_Total', '') %>%
  str_replace_all('.SRS_Total', '')%>%
  str_replace_all('.ASSQ_Total', '')

glmData$clust <- as.factor(glmData$clust)
glmData$Sex <- as.factor(glmData$Sex)

getMeff(glmData[, 3:10])
nyholt.bonf.p(glmData[, 3:10], 0.05)
nyholt_glm <- getMeff(glmData[, 3:10])

## Regression Models
dfASSQ <- glmData %>%
  select(-c(SRS, RBS, SCQ))
dfSRS <- glmData %>%
  select(-c(ASSQ, RBS, SCQ))
dfRBS <- glmData %>%
  select(-c(ASSQ, SRS, SCQ))
dfSCQ <- glmData %>%
  select(-c(ASSQ, SRS, RBS))

glmASSQ <- lm(data = dfASSQ, formula = ASSQ ~  . + clust:.)
glmSRS <- lm(data = dfSRS, formula = SRS ~ . + clust:.)
glmRBS <- lm(data = dfRBS, formula = RBS ~ . + clust:.)
glmSCQ <- lm(data = dfSCQ, formula = SCQ ~ . + clust:.)

summary(glmASSQ)
glance(glmASSQ)
car::vif(glmASSQ)

summary(glmSRS)
glance(glmSRS)
car::vif(glmSRS)

summary(glmRBS)
glance(glmRBS)
car::vif(glmRBS)

summary(glmSCQ)
glance(glmSCQ)
car::vif(glmSCQ)

regression_table <- function(glm) {
  coefficients <- data.frame(summary(glm)$coefficients)
  coefficients$p.cor <- coefficients$Pr...t.. * nyholt_glm
  significant <- data.frame(coefficients[which(coefficients$p.cor <= .05),])
  
  # creating significance stars
  sig.stars <- ifelse(coefficients$p.cor < 0.0001, '***', 
                      ifelse(coefficients$p.cor < 0.001, '*** ', 
                             ifelse(coefficients$p.cor < 0.01, '**  ', 
                                    ifelse(coefficients$p.cor < 0.05, '*   ', '    '))))
  # display coefficients with stars
  coefficients.stars <- coefficients
  coefficients.stars[1:3] <- format(round(coefficients.stars[1:3], digits = 2), nsmall = 2)
  coefficients.stars[4:5] <- format(round(coefficients.stars[4:5], digits = 4), nsmall = 4)
  coefficients.stars$p.cor <- paste(coefficients.stars$p.cor, sig.stars, sep = '')
  colnames(coefficients.stars) <- c('Estimate', 'SE', 't.value', 'p.value', 'p.corr')
  rownames(coefficients.stars) <- rownames(coefficients.stars) %>%
    str_replace_all(':', ' * ') %>%
    str_replace_all('clust2', 'Group') %>%
    str_replace_all('1', '')
  
  results <- list(coef = coefficients, sig = significant, coef.stars = coefficients.stars, aic = glm$aic)
  return(results)
}

regTabASSQ <- regression_table(glmASSQ)
print(regTabASSQ$coef.stars)
print(regTabASSQ$sig)
regTabSRS <- regression_table(glmSRS)
print(regTabSRS$coef.stars)
print(regTabSRS$sig)
regTabRBS <- regression_table(glmRBS)
print(regTabRBS$coef.stars)
print(regTabRBS$sig)
regTabSCQ <- regression_table(glmSCQ)
print(regTabSCQ$coef.stars)
print(regTabSCQ$sig)

# saving model outputs
write.table(regTabASSQ$coef.stars, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/ASSQ.txt', row.names = T, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(regTabSRS$coef.stars, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/SRS.txt', row.names = T, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(regTabRBS$coef.stars, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/RBS.txt', row.names = T, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(regTabSCQ$coef.stars, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/SCQ.txt', row.names = T, col.names = T, quote = F, sep = '\t', dec = '.')

# Table for comorbidities
dx <- c("DX_01", "DX_02", "DX_03", "DX_04", "DX_05", "DX_06", "DX_07", "DX_08", "DX_09", "DX_10")

comorb_ASD1 <- data.frame(table(unlist(apply(allData[allData$clust == 1, dx], MARGIN = 1, FUN = unique))))
comorb_ASD1$Perc <- (comorb_ASD1$Freq / nrow(allData[allData$clust == 1,])) * 100
colnames(comorb_ASD1) <- c('Var1', 'ASD-I', 'Perc1')
comorb_ASD2 <- data.frame(table(unlist(apply(allData[allData$clust == 2, dx], MARGIN = 1, FUN = unique))))
comorb_ASD2$Perc <- (comorb_ASD2$Freq / nrow(allData[allData$clust == 2,])) * 100
colnames(comorb_ASD2) <- c('Var1', 'ASD-II', 'Perc2')

comorbData <- join(comorb_ASD1, comorb_ASD2, by = 'Var1', type = 'full')
comorbData <- comorbData[order(comorbData$`ASD-I`, decreasing = T),]
comorbData <- comorbData[-1, ]

write.table(comorbData, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/Comorbidities.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(head(comorbData, 10), file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/Comorbidities_head.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')

## t-Tests for comparison of group means

ttest_function <- function(measure) {
  ttest1v2 = t.test(allData[allData['clust'] == 1, measure], allData[allData['clust'] == 2, measure])
  ttestASDvC = t.test(allData[allData['diag'] == 'HC', measure], allData[allData['diag'] == 'ASD', measure])
  
  ttestResults <- list(
    t = list(t1v2 = as.numeric(format(round(ttest1v2$statistic, 2), nsmall = 2)), tASDvC = as.numeric(format(round(ttestASDvC$statistic, 2), nsmall = 2))),
    p = list(p1v2 = ttest1v2$p.value, p1v2cor = ttest1v2$p.value * nyholt_behav, pASDvC = ttestASDvC$p.value, pASDvCcor = ttestASDvC$p.value * nyholt_behav),
    df = list(df1v2 = format(round(ttest1v2$parameter, 2), nsmall = 2), dfASDvC = format(round(ttestASDvC$parameter, 2), nsmall = 2))
  )
  
  results <- list(t = ttestResults$t, p = ttestResults$p, df = ttestResults$df)
  return(results)
}

ttestASSQ <- ttest_function(params$BehavMeasures[1])
ttestSRS <- ttest_function(params$BehavMeasures[2])
ttestAWR <- ttest_function(params$SRSsubscales[1])
ttestCOG <- ttest_function(params$SRSsubscales[2])
ttestCOM <- ttest_function(params$SRSsubscales[3])
ttestMOT <- ttest_function(params$SRSsubscales[5])
ttestRRB <- ttest_function(params$SRSsubscales[6])
ttestRBS <- ttest_function(params$BehavMeasures[3])
ttestSCQ <- ttest_function(params$BehavMeasures[4])
ttestAge <- ttest_function(params$BehavMeasures[5])

cohensD_function <- function(measure) {
  cohensD1v2 = cohen.d(allData[allData['clust'] == 1, measure], allData[allData['clust'] == 2, measure])
  cohensDASDvC = cohen.d(allData[allData['diag'] == 'HC', measure], allData[allData['diag'] == 'ASD', measure])
  
  return(list(D1v2 = format(round(cohensD1v2$estimate, digits = 2), nsmall = 2), DASDvC = format(round(cohensDASDvC$estimate, digits = 2), nsmall = 2)))
}

cohensDASSQ <- cohensD_function(params$BehavMeasures[1])
cohensDSRS <- cohensD_function(params$BehavMeasures[2])
cohensDAWR <- cohensD_function(params$SRSsubscales[1])
cohensDCOG <- cohensD_function(params$SRSsubscales[2])
cohensDCOM <- cohensD_function(params$SRSsubscales[3])
cohensDMOT <- cohensD_function(params$SRSsubscales[5])
cohensDRRB <- cohensD_function(params$SRSsubscales[6])
cohensDRBS <- cohensD_function(params$BehavMeasures[3])
cohensDSCQ <- cohensD_function(params$BehavMeasures[4])
cohensDAge <- cohensD_function(params$BehavMeasures[5])

# Tables
# full ASD v C
ttestTableASDvC <- data.frame(Measure = c('Age', 'ASSQ', 'SRS', 'AWR', 'COG', 'COM', 'MOT', 'RRB', 'RBS', 'SCQ'), 
                            Mean.Ctrl = c(statsAge$Mean$Total_Ctrls, statsASSQ$Mean$Total_Ctrls, statsSRS$Mean$Total_Ctrls, statsAWR$Mean$Total_Ctrls, statsCOG$Mean$Total_Ctrls, statsCOM$Mean$Total_Ctrls, statsMOT$Mean$Total_Ctrls, statsRRB$Mean$Total_Ctrls, statsRBS$Mean$Total_Ctrls, statsSCQ$Mean$Total_Ctrls),
                            SD.Ctrl = c(statsAge$SD$Total_Ctrls, statsASSQ$SD$Total_Ctrls, statsSRS$SD$Total_Ctrls, statsAWR$SD$Total_Ctrls, statsCOG$SD$Total_Ctrls, statsCOM$SD$Total_Ctrls, statsMOT$SD$Total_Ctrls, statsRRB$SD$Total_Ctrls, statsRBS$SD$Total_Ctrls, statsSCQ$SD$Total_Ctrls),
                            Mean.ASD = c(statsAge$Mean$Total_ASD, statsASSQ$Mean$Total_Ctrls, statsSRS$Mean$Total_ASD, statsAWR$Mean$Total_ASD, statsCOG$Mean$Total_ASD, statsCOM$Mean$Total_ASD, statsMOT$Mean$Total_ASD, statsRRB$Mean$Total_ASD, statsRBS$Mean$Total_ASD, statsSCQ$Mean$Total_ASD),
                            SD.ASD = c(statsAge$SD$Total_ASD, statsASSQ$SD$Total_ASD, statsSRS$SD$Total_ASD, statsAWR$SD$Total_ASD, statsCOG$SD$Total_ASD, statsCOM$SD$Total_ASD, statsMOT$SD$Total_ASD, statsRRB$SD$Total_ASD, statsRBS$SD$Total_ASD, statsSCQ$SD$Total_ASD),
                            CohensD = c(cohensDAge$DASDvC, cohensDASSQ$DASDvC, cohensDSRS$DASDvC, cohensDAWR$DASDvC, cohensDCOG$DASDvC, cohensDCOM$DASDvC, cohensDMOT$DASDvC, cohensDRRB$DASDvC, cohensDRBS$DASDvC, cohensDSCQ$DASDvC),
                            t = c(ttestAge$t$tASDvC, ttestASSQ$t$tASDvC, ttestSRS$t$tASDvC, ttestAWR$t$tASDvC, ttestCOG$t$tASDvC, ttestCOM$t$tASDvC, ttestMOT$t$tASDvC, ttestRRB$t$tASDvC, ttestRBS$t$tASDvC, ttestSCQ$t$tASDvC),
                            df = c(ttestAge$df$dfASDvC, ttestASSQ$df$dfASDvC, ttestSRS$df$dfASDvC, ttestAWR$df$dfASDvC, ttestCOG$df$dfASDvC, ttestCOM$df$dfASDvC, ttestMOT$df$dfASDvC, ttestRRB$df$dfASDvC, ttestRBS$df$dfASDvC, ttestSCQ$df$dfASDvC),
                            p = c(ttestAge$p$pASDvC, ttestASSQ$p$pASDvC, ttestSRS$p$pASDvC, ttestAWR$p$pASDvC, ttestCOG$p$pASDvC, ttestCOM$p$pASDvC, ttestMOT$p$pASDvC, ttestRRB$p$pASDvC, ttestRBS$p$pASDvC, ttestSCQ$p$pASDvC), 
                            p.corr = c(ttestAge$p$pASDvCcor, ttestASSQ$p$pASDvCcor, ttestSRS$p$pASDvCcor, ttestAWR$p$pASDvCcor, ttestCOG$p$pASDvCcor, ttestCOM$p$pASDvCcor, ttestMOT$p$pASDvCcor, ttestRRB$p$pASDvCcor, ttestRBS$p$pASDvCcor, ttestSCQ$p$pASDvCcor))

colnames(ttestTableASDvC) <- c('Measure', 'Controls M', 'SD', 'ASD M','SD', 'd', 't', 'df', 'p', 'p*corr')

ttestTableASDvC$p <- ifelse(ttestTableASDvC$p < .0001, format(round(ttestTableASDvC$p, 3), nsmall = 3), 
                                    ifelse(ttestTableASDvC$p < .001, format(round(ttestTableASDvC$p, 3), nsmall = 3),
                                           ifelse(ttestTableASDvC$p < .01, paste(round(ttestTableASDvC$p, 3), '** ', sep = ''), 
                                                  ifelse(ttestTableASDvC$p < .05, round(ttestTableASDvC$p, 3), round(ttestTableASDvC$p, 2)))))

ttestTableASDvC$`p*corr` <- ifelse(ttestTableASDvC$`p*corr` < .0001, '< .001***', 
                               ifelse(ttestTableASDvC$`p*corr` < .001, '< .001***',
                                      ifelse(ttestTableASDvC$`p*corr` < .01, paste(round(ttestTableASDvC$`p*corr`, 4), '** ', sep = ''), 
                                             ifelse(ttestTableASDvC$`p*corr` < .05, paste(round(ttestTableASDvC$`p*corr`, 4), '*  ', sep = ''), paste(round(ttestTableASDvC$`p*corr`, 3), '   ', sep = '')))))


write.table(ttestTableASDvC, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/tTestASDvC.txt', row.names = F, quote = F, col.names = T, sep = '\t', dec = '.')

# ASD1 v ASD2
ttestTable1v2 <- data.frame(Measure = c('Age', 'ASSQ', 'SRS', 'AWR', 'COG', 'COM', 'MOT', 'RRB', 'RBS', 'SCQ'), 
                         Mean.ASD1 = c(statsAge$Mean$ASD1, statsASSQ$Mean$ASD1, statsSRS$Mean$ASD1, statsAWR$Mean$ASD1, statsCOG$Mean$ASD1, statsCOM$Mean$ASD1, statsMOT$Mean$ASD1, statsRRB$Mean$ASD1, statsRBS$Mean$ASD1, statsSCQ$Mean$ASD1),
                         SD.ASD1 = c(statsAge$SD$ASD1, statsASSQ$SD$ASD1, statsSRS$SD$ASD1, statsAWR$SD$ASD1, statsCOG$SD$ASD1, statsCOM$SD$ASD1, statsMOT$SD$ASD1, statsRRB$SD$ASD1, statsRBS$SD$ASD1, statsSCQ$SD$ASD1),
                         Mean.ASD2 = c(statsAge$Mean$ASD2, statsASSQ$Mean$ASD2, statsSRS$Mean$ASD2, statsAWR$Mean$ASD2, statsCOG$Mean$ASD2, statsCOM$Mean$ASD2, statsMOT$Mean$ASD2, statsRRB$Mean$ASD2, statsRBS$Mean$ASD2, statsSCQ$Mean$ASD2),
                         SD.ASD2 = c(statsAge$SD$ASD2, statsASSQ$SD$ASD2, statsSRS$SD$ASD2, statsAWR$SD$ASD2, statsCOG$SD$ASD2, statsCOM$SD$ASD2, statsMOT$SD$ASD2, statsRRB$SD$ASD2, statsRBS$SD$ASD2, statsSCQ$SD$ASD2),
                         CohensD = c(cohensDAge$D1v2, cohensDASSQ$D1v2, cohensDSRS$D1v2, cohensDAWR$D1v2, cohensDCOG$D1v2, cohensDCOM$D1v2, cohensDMOT$D1v2, cohensDRRB$D1v2, cohensDRBS$D1v2, cohensDSCQ$D1v2),
                         t = c(ttestAge$t$t1v2, ttestASSQ$t$t1v2, ttestSRS$t$t1v2, ttestAWR$t$t1v2, ttestCOG$t$t1v2, ttestCOM$t$t1v2, ttestMOT$t$t1v2, ttestRRB$t$t1v2, ttestRBS$t$t1v2, ttestSCQ$t$t1v2),
                         df = c(ttestAge$df$df1v2, ttestASSQ$df$df1v2, ttestSRS$df$df1v2, ttestAWR$df$df1v2, ttestCOG$df$df1v2, ttestCOM$df$df1v2, ttestMOT$df$df1v2, ttestRRB$df$df1v2, ttestRBS$df$df1v2, ttestSCQ$df$df1v2),
                         p = c(ttestAge$p$p1v2, ttestASSQ$p$p1v2, ttestSRS$p$p1v2, ttestAWR$p$p1v2, ttestCOG$p$p1v2, ttestCOM$p$p1v2, ttestMOT$p$p1v2, ttestRRB$p$p1v2, ttestRBS$p$p1v2, ttestSCQ$p$p1v2),
                         p.corr = c(ttestAge$p$p1v2cor, ttestASSQ$p$p1v2cor, ttestSRS$p$p1v2cor, ttestAWR$p$p1v2cor, ttestCOG$p$p1v2cor, ttestCOM$p$p1v2cor, ttestMOT$p$p1v2cor, ttestRRB$p$p1v2cor, ttestRBS$p$p1v2cor, ttestSCQ$p$p1v2cor))

colnames(ttestTable1v2) <- c('Measure', 'ASD-I M', 'SD', 'ASD-II M', 'SD', 'd', 't', 'df', 'p', 'p*corr')

ttestTable1v2$p <- ifelse(ttestTable1v2$p < .0001, format(round(ttestTable1v2$p, 3), nsmall = 3), 
                            ifelse(ttestTable1v2$p < .001, format(round(ttestTable1v2$p, 3), nsmall = 3),
                                   ifelse(ttestTable1v2$p < .01, paste(round(ttestTable1v2$p, 3), '** ', sep = ''), 
                                          ifelse(ttestTable1v2$p < .05, round(ttestTable1v2$p, 3), round(ttestTable1v2$p, 2)))))

ttestTable1v2$`p*corr` <- ifelse(ttestTable1v2$`p*corr` < .0001, '< .001***', 
                                   ifelse(ttestTable1v2$`p*corr` < .001, '< .001***',
                                          ifelse(ttestTable1v2$`p*corr` < .01, paste(round(ttestTable1v2$`p*corr`, 4), '** ', sep = ''), 
                                                 ifelse(ttestTable1v2$`p*corr` < .05, paste(round(ttestTable1v2$`p*corr`, 4), '*  ', sep = ''), round(ttestTable1v2$`p*corr`, 3)))))

write.table(ttestTable1v2, file = '~/Documents/Psychologie/MA/Results/Behavioral/Tables/tTest1v2.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')

