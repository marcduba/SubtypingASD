### Hypothesis 1: COMPARISON OF ASYMMETRY INDEXES BTW. CLUSTERS/CONTROLS

library(dplyr)
library(tidyr)
library(stringr)
library(jtools)
library(ggplot2)
library(ggpubr)
library(lsr)
library(effsize)
library(fsbrain)
library(rgl)
library(freesurferformats)
library(xtable)
library(data.table)
library(Hmisc)
library(clusterpval)
library(sm)
source('~/Documents/Psychologie/MA/Code/R/GeneralParameters.R')

setwd(params$FreeSurferDir)
template_subject = params$template_sub
template_subject_dir = params$template_sub_dir

# Load data
setwd(params$WorkingDirectory)
allData <- read.csv('./Data/allData.csv')

# Separating datase into variables that are needed for analysis
VarNames <- colnames(allData[2:152])
VarNames <- c(VarNames, params$BehavMeasures[5], params$BehavMeasures[6])
braindata <- allData[, VarNames]
braindata$diag <-  as.factor(braindata$diag)
braindata$clust <- as.factor(braindata$clust)
braindata$Basic_Demos.Sex <- as.factor(braindata$Basic_Demos.Sex)

# Getting FreeSurfer regionnames of Destrieux atlas
FSregionnames <- get.atlas.region.names(atlas = 'aparc.a2009s', 
                                      template_subjects_dir = template_subject_dir, 
                                      template_subject = template_subject, 
                                      hemi = 'lh')

FSregionnames <- FSregionnames[FSregionnames != 'Unknown']
FSregionnames <- FSregionnames[FSregionnames != 'Medial_wall']

# Separating df for CT & SA
allSA <- braindata %>% 
  select(c('EID', 'diag', 'clust', ends_with('_SA')))
colnames(allSA) <- colnames(allSA) %>%
  str_replace_all('_SA', '') %>%
  str_replace_all('[.]','-')

allCT <- braindata %>%
  select(c('EID', 'diag', 'clust', ends_with('_CT')))
colnames(allCT) <- colnames(allCT) %>%
  str_replace_all('_CT', '') %>%
  str_replace_all('[.]', '-')

# Creating subsets for groups
CT_HC <- allCT %>%
  filter(diag == 'HC')
CT_ASD <- allCT %>%
  filter(diag == 'ASD')
CT_ASD1 <- allCT %>%
  filter(clust == 1)
CT_ASD2 <- allCT %>%
  filter(clust == 2)
CT_ASD3 <- allCT %>% 
  filter(clust == 3)
CT_12C <- allCT %>%
  filter(clust == 0 | clust == 1 | clust == 2)

SA_HC <- allSA %>%
  filter(diag == 'HC')
SA_ASD <- allSA %>% 
  filter(diag == 'ASD')
SA_ASD1 <- allSA %>%
  filter(clust == 1)
SA_ASD2 <- allSA %>% 
  filter(clust == 2)
SA_ASD3 <- allSA %>% 
  filter(clust == 3)
SA_12C <- allSA %>%
  filter(clust == 0 | clust == 1 | clust == 2)


# Correlations between features (CT & SA regions)
corrmatrix_between_function <- function(dfCT, dfSA) {
  df1 = dfCT %>%
    select(-c('EID', 'diag', 'clust'))
  colnames(df1) = colnames(df1) %>%
    paste('_CT', sep = '')
  df2 = dfSA %>%
    select(-c('EID', 'diag', 'clust'))
  colnames(df2) = colnames(df2) %>%
    paste('_SA', sep = '')
  matrix1 = as.matrix(df1)
  matrix2 = as.matrix(df2)
  
  corrmat = rcorr(matrix1, matrix2, type = 'pearson')
  corrmatR = round(corrmat$r, digits = 2)
  corrmatP = round(corrmat$P, digits = 4)
  
  sig.stars = ifelse(corrmatP < .0001, "***", ifelse(corrmatP < .001, "*** ", ifelse(corrmatP < .01, "**  ", ifelse(corrmatP < .05, "*   ", "    "))))
  
  corrmat.sig = data.frame(matrix(paste(corrmatR, sig.stars, sep = ''), ncol = ncol(corrmatR)))
  between = data.frame(corrmat.sig[75:148, 1:74])
  rownames(between) = rownames(corrmatR)[75:148]
  colnames(between) = colnames(corrmatR)[1:74]
  diag = diag(as.matrix(between))
  diagmatrix = data.frame(colnames(between), diag)
  colnames(diagmatrix) = c('Region', 'r')
  diagmatrix$Region = diagmatrix$Region %>%
    str_replace_all('_CT', '')

  matrices = list(b = between, d = diag, mat = diagmatrix)
  return(matrices)
}

corrCT.SA <- corrmatrix_between_function(allCT, allSA)
corrCT.SA$mat

setwd(params$WorkingDirectory)
write.table(corrCT.SA$mat, file = './Results/Anatomical/Tables/corrCT.SA.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')


# Checking distribution of brain regions
# CT
regionsCT <- c(colnames(CT_12C[4:77]))
qqCT <- list()

for(i in regionsCT) {
  x <- CT_12C[,i]
  qqCT[[i]] <- qqnorm(x)
  qqCT[[i]][[i]] <- qqline(x, col = 2)
  Sys.sleep(1)

}

# SA
regionsSA <- c(colnames(SA_12C[4:77]))
qqSA <- list()

for(i in regionsSA) {
  x <- SA_12C[,i]
  qqSA[[i]] <- qqnorm(x)
  qqSA[[i]][[i]] <- qqline(x, col = 2)
  Sys.sleep(1)
}

# Nyholt correction
getMeff(allCT[,4:77])
nyholt.bonf.p(allCT[, 4:77], 0.05)
nyholt.bonf_CT <- getMeff(allCT[,4:77])

getMeff(allSA[, 4:77])
nyholt.bonf.p(allSA[, 4:77], 0.05)
nyholt.bonf_SA <- getMeff(allSA[,4:77])

### tTest for differences in group means of brain regions

# Function for tTest
ttest_brain_function <- function(data, clust1, clust2, nyholt) {
  dfXY = data %>% 
    select(c(FSregionnames, 'clust'))
  meanX = NULL
  meanY = NULL
  sdX = NULL
  sdY = NULL
  t.value = NULL
  degf = NULL
  p.value = NULL
  cohen.d = NULL
  
  for(i in 1:(ncol(dfXY)-1)) {
    meanX = as.data.frame(rbind(meanX, t.test(dfXY[dfXY['clust'] == clust1, i], dfXY[dfXY['clust'] == clust2, i])$estimate[1]))
    meanY = as.data.frame(rbind(meanY, t.test(dfXY[dfXY['clust'] == clust1, i], dfXY[dfXY['clust'] == clust2, i])$estimate[2]))
    sdX = as.data.frame(rbind(sdX, sd(dfXY[dfXY['clust' == clust1, i]])))
    sdY = as.data.frame(rbind(sdY, sd(dfXY[dfXY['clust' == clust2, i]])))
    t.value = as.data.frame(rbind(t.value, t.test(dfXY[dfXY['clust'] == clust1, i], dfXY[dfXY['clust'] == clust2, i])$statistic))
    degf = as.data.frame(rbind(degf, t.test(dfXY[dfXY['clust'] == clust1, i], dfXY[dfXY['clust'] == clust2, i])$parameter))
    p.value = as.data.frame(rbind(p.value, t.test(dfXY[dfXY['clust'] == clust1, i], dfXY[dfXY['clust'] == clust2, i])$p.value))
    cohen.d = as.data.frame(rbind(cohen.d, cohen.d(dfXY[dfXY['clust'] == clust1, i], dfXY[dfXY['clust'] == clust2, i])$estimate))
    
  }
  table = data.frame(FSregionnames, meanX, sdX, meanY, sdY, cohen.d, t.value, degf, p.value)
  
  colnames(table) = c('FSregion', 'MeanX', 'SDX', 'MeanY', 'SDY', 'CohenD', 't', 'df', 'p')
  table$p.cor = table$p * nyholt
  table = table[1:74, ]
  table_sig = table[which(table$p.cor <= 0.05), ]
  
  ttestData = list(fullTable = table, sigTable = table_sig, df = dfXY)
  
  return(ttestData)
}

# Creating function for scatter plot
limsCT <- c(-.3, .3)
limsSA <- c(-1, 1)

scatterplot_func <- function(df, table, grp1, grp2, col1, col2, lims) {
  # Prepare dataframe
  dfAB <- df %>%
    dplyr::select(c('clust', matches(table$sigTable$FSregion))) %>%
    filter(clust == grp1 | clust == grp2)
  dfAB$clust <- ifelse(dfAB$clust == 1, 'ASD-I', ifelse(dfAB$clust == 2, 'ASD-II', 'Controls'))
  dfAB$clust <- as.factor(dfAB$clust)
  
  dfCD <- melt(dfAB, id.vars = 'clust')
  
  # Create plot
  plot <- ggplot(dfCD, aes(x = variable, y = value, fill = clust)) +
    geom_jitter(aes(col = clust, fill = clust), position = position_dodge2(width = 0.75), alpha = 1) +
    geom_boxplot(fatten = NULL, position = position_dodge2(width = 0.75), alpha = 0.1, outlier.alpha = 0) + 
    stat_summary(fun.y = mean, geom = 'errorbar', aes(ymax = ..y.., ymin = ..y.., width = 0.75), position = position_dodge2(width = 0.75)) +
    theme_apa(remove.x.gridlines = F) + 
    xlab('Regions') + 
    ylab('Mean AI') + 
    lims(y = lims) +
    coord_flip() +
    scale_color_manual(values = c(col1, col2)) +
    scale_fill_manual(values = c(col1, col2)) + 
    theme(text = element_text(family = 'sans'))
  
  return(plot)
}

## Cortical Thickness

# all ASD vs Controls
dfXY = allCT %>% 
  select(c(FSregionnames, 'diag'))
meanX = NULL
meanY = NULL
sdX = NULL
sdY = NULL
t.value = NULL
degf = NULL
p.value = NULL
cohen.d = NULL

for(i in 1:(ncol(dfXY)-1)) {
  meanX = as.data.frame(rbind(meanX, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$estimate[1]))
  meanY = as.data.frame(rbind(meanY, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$estimate[2]))
  sdX = as.data.frame(rbind(sd(dfXY['diag' == 'HC', i], na.rm = T)))
  sdY = as.data.frame(rbind(sd(dfXY['diag' == 'ASD', i], na.rm = T)))
  t.value = as.data.frame(rbind(t.value, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$statistic))
  degf = as.data.frame(rbind(degf, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$parameter))
  p.value = as.data.frame(rbind(p.value, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$p.value))
  cohen.d = as.data.frame(rbind(cohen.d, cohen.d(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$estimate))
  
}
ttestCT_AvC = data.frame(FSregionnames, meanX, sdX, meanY, sdY, cohen.d, t.value, degf, p.value)
colnames(ttestCT_AvC) = c('FSregionnames', 'MeanX', 'SDX', 'MeanY', 'SDY', 'CohenD', 't', 'df', 'p')
ttestCT_AvC$pcor = ttestCT_AvC$p * nyholt.bonf_CT
ttestCT_AvCsig = ttestCT_AvC[which(ttestCT_AvC$pcor <= 0.05), ]

# ASD1 v HC CT
ttestCT_1vC <- ttest_brain_function(allCT, 0, 1, nyholt.bonf_CT)
ttestCT_1vC$sigTable
ttestCT_1vC$fullTable
scatterCT_1vC <- scatterplot_func(allCT, ttestCT_1vC, 1, 0, params$paletteASD[1], 'grey', limsCT)

# ASD2 v HC CT
ttestCT_2vC <- ttest_brain_function(allCT, 0, 2, nyholt.bonf_CT)
ttestCT_2vC$sigTable
ttestCT_2vC$fullTable
scatterCT_2vC <- scatterplot_func(allCT, ttestCT_2vC, 2, 0, params$paletteASD[3], 'grey', limsCT)

# ASD1 v ASD2 CT
ttestCT_1v2 <- ttest_brain_function(allCT, 1, 2, nyholt.bonf_CT)
ttestCT_1v2$sigTable
ttestCT_1v2$fullTable
scatterCT_1v2 <- scatterplot_func(allCT, ttestCT_1v2, 1, 2, params$paletteASD[1], params$paletteASD[3], limsCT)


## Surface Area

# all ASD vs Controls
dfXY = allSA %>% 
  select(c(FSregionnames, 'diag'))
meanX = NULL
meanY = NULL
t.value = NULL
degf = NULL
p.value = NULL
cohen.d = NULL

for(i in 1:(ncol(dfXY)-1)) {
  meanX = as.data.frame(rbind(meanX, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$estimate[1]))
  meanY = as.data.frame(rbind(meanY, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$estimate[2]))
  t.value = as.data.frame(rbind(t.value, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$statistic))
  degf = as.data.frame(rbind(degf, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$parameter))
  p.value = as.data.frame(rbind(p.value, t.test(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$p.value))
  cohen.d = as.data.frame(rbind(cohen.d, cohen.d(dfXY[dfXY['diag'] == 'HC', i], dfXY[dfXY['diag'] == 'ASD', i])$estimate))
  
}
ttestSA_AvC = data.frame(FSregionnames, meanX, meanY, cohen.d, t.value, degf, p.value)
colnames(ttestSA_AvC) = c('FSregionnames', 'MeanX', 'MeanY', 'CohenD', 't', 'df', 'p')
ttestSA_AvC$pcor = ttestSA_AvC$p * nyholt.bonf_SA
ttestSA_AvCsig = ttestSA_AvC[which(ttestSA_AvC$pcor <= 0.05), ]

# ASD1 v HC SA
ttestSA_1vC <- ttest_brain_function(allSA, 0, 1, nyholt.bonf_SA)
ttestSA_1vC$sigTable
ttestSA_1vC$fullTable
scatterSA_1vC <- scatterplot_func(allSA, ttestSA_1vC, 1, 0, params$paletteASD[1], 'grey', limsSA)

# ASD2 v HC SA
ttestSA_2vC <- ttest_brain_function(allSA, 0, 2, nyholt.bonf_SA)
ttestSA_2vC$sigTable
ttestSA_2vC$fullTable
scatterSA_2vC <- scatterplot_func(allSA, ttestSA_2vC, 2, 0, params$paletteASD[3], 'grey', limsSA)

# ASD1 v ASD2 SA
ttestSA_1v2 <- ttest_brain_function(allSA, 1, 2, nyholt.bonf_SA)
ttestSA_1v2$sigTable
ttestSA_1v2$fullTable
scatterSA_1v2 <- scatterplot_func(allSA, ttestSA_1v2, 1, 2, params$paletteASD[1], params$paletteASD[3], limsSA)

# Saving list of significant brain regions to .rda file
significantRegions = list(
  CT = list('ASD1vC' = ttestCT_1vC$sigTable$FSregion, 'ASD2vC' = ttestCT_2vC$sigTable$FSregion, 'ASD1v2' = ttestCT_1v2$sigTable$FSregion, 'ASDvC' = ttestCT_AvCsig$FSregion),
  SA = list('ASD1vC' = ttestSA_1vC$sigTable$FSregion, 'ASD2vC' = ttestSA_2vC$sigTable$FSregion, 'ASD1v2' = ttestSA_1v2$sigTable$FSregion, 'ASDvC' = ttestSA_AvCsig$FSregion)
)
save(significantRegions, file = '~/Documents/Psychologie/MA/Code/R/SignificantRegions.Rda')

# Saving tTest Tables and Scatter Plots
formating_function <- function(table, grp1, grp2) {
  table[c(2:3, 7:8)] <- format(round(table[c(2:3, 7:8)], digits = 4), nsmall = 4, scientific = F)
  table[4:6] <- format(round(table[4:6], digits = 2), nsmall = 2, scientific = F)
  colnames(table) <- c('Region', grp1, grp2, 'Cohen`s D', 't', 'df', 'p', 'pcor')
  
  return(table)
}

niceCT_AvC <- formating_function(ttestCT_AvC, 'Mean Controls', 'Mean ASD')
niceCT_1vC <- formating_function(ttestCT_1vC$fullTable, 'Mean Controls', 'Mean ASD-I')
niceCT_1vCsig <- formating_function(ttestCT_1vC$sigTable, 'Mean Controls', 'Mean ASD-I')
niceCT_2vC <- formating_function(ttestCT_2vC$fullTable, 'Mean Controls', 'Mean ASD-II')
niceCT_2vCsig <- formating_function(ttestCT_2vC$sigTable, 'Mean Controls', 'Mean ASD-II')
niceCT_1v2 <- formating_function(ttestCT_1v2$fullTable, 'Mean ASD-I', 'Mean ASD-II')
niceCT_1v2sig <- formating_function(ttestCT_1v2$sigTable, 'Mean ASD-I', 'Mean ASD-II')

niceSA_AvC <- formating_function(ttestSA_AvC, 'Mean Controls', 'Mean ASD')
niceSA_1vC <- formating_function(ttestSA_1vC$fullTable, 'Mean Controls', 'Mean ASD-I')
niceSA_1vCsig <- formating_function(ttestSA_1vC$sigTable, 'Mean Controls', 'Mean ASD-I')
niceSA_2vC <- formating_function(ttestSA_2vC$fullTable, 'Mean Controls', 'Mean ASD-II')
niceSA_2vCsig <- formating_function(ttestSA_2vC$sigTable, 'Mean Controls', 'Mean ASD-II')
niceSA_1v2 <- formating_function(ttestSA_1v2$fullTable, 'Mean ASD-I', 'Mean ASD-II')
niceSA_1v2sig <- formating_function(ttestSA_1v2$sigTable, 'Mean ASD-I', 'Mean ASD-II')

# full results from tTests
write.table(niceCT_AvC, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestCT_AvC.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(niceCT_1vC, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestCT_1vC.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(niceCT_2vC, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestCT_2vC.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(niceCT_1v2, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestCT_1v2.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')

write.table(niceSA_AvC, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestSA_AvC.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(niceSA_1vC, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestSA_1vC.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(niceSA_2vC, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestSA_2vC.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(niceSA_1v2, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestSA_1v2.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')

# significant results from tTests
write.table(niceCT_1v2sig, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestCT_1v2sig.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')

write.table(niceSA_1vCsig, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestSA_1vCsig.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(niceSA_2vCsig, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestSA_2vCsig.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')
write.table(niceSA_1v2sig, file = '~/Documents/Psychologie/MA/Results/Anatomical/Tables/ttestSA_1v2sig.txt', row.names = F, col.names = T, quote = F, sep = '\t', dec = '.')

# Scatter plots
ggsave(scatterCT_1v2, filename = 'scatterCT_1v2.png', path = '~/Documents/Psychologie/MA/Results/Anatomical/Figures/', device = 'png', width = 1200, height = 600, units = 'px', dpi = 150)

ggsave(scatterSA_1vC, filename = 'scatterSA_1vC.png', path = '~/Documents/Psychologie/MA/Results/Anatomical/Figures/', device = 'png', width = 1200, height = 600, units = 'px', dpi = 150)
ggsave(scatterSA_2vC, filename = 'scatterSA_2vC.png', path = '~/Documents/Psychologie/MA/Results/Anatomical/Figures/', device = 'png', width = 1200, height = 600, units = 'px', dpi = 150)
ggsave(scatterSA_1v2, filename = 'scatterSA_1v2.png', path = '~/Documents/Psychologie/MA/Results/Anatomical/Figures/', device = 'png', width = 1200, height = 600, units = 'px', dpi = 150)


## Visualization of Asymmetry

# Setting visualisation options
rgl_options = list('windowRect' = c(20, 20, 1800, 1200))

# Setting color functions
colFn_raw = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(10, name="RdBu")))

# MEAN ASYMMETRY FOR EACH GROUP FOR ALL REGIONS
mean_AI_function <- function(df) {
  raw = data.frame(colMeans(df[, 4:77], na.rm = T))
  colnames(raw) = 'Mean_AI'
  
  return(raw)
}

rawCT_HC <- mean_AI_function(CT_HC)
rawCT_ASD <- mean_AI_function(CT_ASD)
rawCT_ASD1 <- mean_AI_function(CT_ASD1)
rawCT_ASD2 <- mean_AI_function(CT_ASD2)
rawCT_ASD3 <- mean_AI_function(CT_ASD3)

rawSA_HC <- mean_AI_function(SA_HC)
rawSA_ASD <- mean_AI_function(SA_ASD)
rawSA_ASD1 <- mean_AI_function(SA_ASD1)
rawSA_ASD2 <- mean_AI_function(SA_ASD2)
rawSA_ASD3 <- mean_AI_function(SA_ASD3)

range(c(rawCT_HC$Mean_AI, rawCT_ASD$Mean_AI, rawCT_ASD1$Mean_AI, rawCT_ASD2$Mean_AI, rawCT_ASD3$Mean_AI))
rangeCT <- c(-0.35, 0.35)

range(c(rawSA_HC$Mean_AI, rawSA_ASD$Mean_AI, rawSA_ASD1$Mean_AI, rawSA_ASD2$Mean_AI, rawSA_ASD3$Mean_AI))
rangeSA <- c(-0.55, 0.55) 

# Cortical Thickness
mkc_rawCT = list('colFn' = colFn_raw, 'range' = rangeCT, 'symm' = T)

# HC
regions_CTHC <- c(rawCT_HC$Mean_AI)
names(regions_CTHC) <- FSregionnames

bplot_CTHC <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                           subject_id = template_subject, 
                                           atlas = params$atlas, 
                                           lh_region_value_list = regions_CTHC,
                                           rh_region_value_list = -(regions_CTHC),
                                           surface = 'white', 
                                           draw_colorbar = T, 
                                           makecmap_options = mkc_rawCT, 
                                           views = c('sd_lateral_lh', 'sd_medial_lh'), 
                                           rgloptions = rgl_options, 
                                           bg = 'curv')
# ASD
regions_CTASD <- c(rawCT_ASD$Mean_AI)
names(regions_CTASD) <- FSregionnames

bplot_CTASD <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                           subject_id = template_subject, 
                                           atlas = params$atlas, 
                                           lh_region_value_list = regions_CTASD,
                                           rh_region_value_list = -(regions_CTASD),
                                           surface = 'white', 
                                           draw_colorbar = T, 
                                           makecmap_options = mkc_rawCT, 
                                           views = c('sd_lateral_lh', 'sd_medial_lh'),
                                           rgloptions = rgl_options,
                                           bg = 'curv')
# ASD1
regions_CTASD1 <- c(rawCT_ASD1$Mean_AI)
names(regions_CTASD1) <- FSregionnames

bplot_CTASD1 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                            subject_id = template_subject, 
                                            atlas = params$atlas, 
                                            lh_region_value_list = regions_CTASD1,
                                            rh_region_value_list = -(regions_CTASD1),
                                            surface = 'white', 
                                            draw_colorbar = T, 
                                            makecmap_options = mkc_rawCT, 
                                            views = c('sd_lateral_lh', 'sd_medial_lh'),
                                            rgloptions = rgl_options,
                                            bg = 'curv')
# ASD2
regions_CTASD2 <- c(rawCT_ASD2$Mean_AI)
names(regions_CTASD2) <- FSregionnames

bplot_CTASD2 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                            subject_id = template_subject, 
                                            atlas = params$atlas, 
                                            lh_region_value_list = regions_CTASD2,
                                            rh_region_value_list = -(regions_CTASD2),
                                            surface = 'white', 
                                            draw_colorbar = T, 
                                            makecmap_options = mkc_rawCT, 
                                            views = c('sd_lateral_lh', 'sd_medial_lh'),
                                            rgloptions = rgl_options,
                                            bg = 'curv')
# ASD3
regions_CTASD3 <- c(rawCT_ASD3$Mean_AI)
names(regions_CTASD3) <- FSregionnames

bplot_CTASD3 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                            subject_id = template_subject, 
                                            atlas = params$atlas, 
                                            lh_region_value_list = regions_CTASD3,
                                            rh_region_value_list = -(regions_CTASD3),
                                            surface = 'white', 
                                            draw_colorbar = T, 
                                            makecmap_options = mkc_rawCT, 
                                            views = c('sd_lateral_lh', 'sd_medial_lh'),
                                            rgloptions = rgl_options,
                                            bg = 'curv')

## Surface Area
mkc_rawSA = list('colFn' = colFn_raw, 'range' = rangeSA, 'symm' = T)

# HC
regions_SAHC <- c(rawSA_HC$Mean_AI)
names(regions_SAHC) <- FSregionnames

bplot_SAHC <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                            subject_id = template_subject, 
                                            atlas = params$atlas, 
                                            lh_region_value_list = regions_SAHC,
                                            rh_region_value_list = -(regions_SAHC),
                                            surface = 'white', 
                                            draw_colorbar = T, 
                                            makecmap_options = mkc_rawSA, 
                                            views = c('sd_lateral_lh', 'sd_medial_lh'),
                                            rgloptions = rgl_options,
                                            bg = 'curv')
# ASD
regions_SAASD <- c(rawSA_ASD$Mean_AI)
names(regions_SAASD) <- FSregionnames

bplot_SAASD <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                           subject_id = template_subject, 
                                           atlas = params$atlas, 
                                           lh_region_value_list = regions_SAASD,
                                           rh_region_value_list = -(regions_SAASD),
                                           surface = 'white', 
                                           draw_colorbar = T, 
                                           makecmap_options = mkc_rawSA, 
                                           views = c('sd_lateral_lh', 'sd_medial_lh'), 
                                           rgloptions = rgl_options,
                                           bg = 'curv')
# ASD1
regions_SAASD1 <- c(rawSA_ASD1$Mean_AI)
names(regions_SAASD1) <- FSregionnames

bplot_SAASD1 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                           subject_id = template_subject, 
                                           atlas = params$atlas, 
                                           lh_region_value_list = regions_SAASD1,
                                           rh_region_value_list = -(regions_SAASD1),
                                           surface = 'white', 
                                           draw_colorbar = T, 
                                           makecmap_options = mkc_rawSA, 
                                           views = c('sd_lateral_lh', 'sd_medial_lh'),
                                           rgloptions = rgl_options,
                                           bg = 'curv')
# ASD2
regions_SAASD2 <- c(rawSA_ASD2$Mean_AI)
names(regions_SAASD2) <- FSregionnames

bplot_SAASD2 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                           subject_id = template_subject, 
                                           atlas = params$atlas, 
                                           lh_region_value_list = regions_SAASD2,
                                           rh_region_value_list = -(regions_SAASD2),
                                           surface = 'white', 
                                           draw_colorbar = T, 
                                           makecmap_options = mkc_rawSA, 
                                           views = c('sd_lateral_lh', 'sd_medial_lh'),
                                           rgloptions = rgl_options,
                                           bg = 'curv')
# ASD3
regions_SAASD3 <- c(rawSA_ASD3$Mean_AI)
names(regions_SAASD3) <- FSregionnames

bplot_SAASD3 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                           subject_id = template_subject, 
                                           atlas = params$atlas, 
                                           lh_region_value_list = regions_SAASD3,
                                           rh_region_value_list = -(regions_SAASD3),
                                           surface = 'white', 
                                           draw_colorbar = T, 
                                           makecmap_options = mkc_rawSA, 
                                           views = c('sd_lateral_lh', 'sd_medial_lh'),
                                           rgloptions = rgl_options,
                                           bg = 'curv')


# separate legend
legendRawCT <- coloredmesh.plot.colorbar.separate(bplot_CTASD, show = T, image.plot_extra_options = list('horizontal' = T, 'legend.width' = 2, 'legend.mar' = 14.1, 'legend.only' = T, 'col' = colFn_raw('n' = 11L), 'legend.lab' = 'Mean Asymmetry Index (L > R)', 'legend.cex' = 1.2, 'legend.line' = 4))
legendRawSA <- coloredmesh.plot.colorbar.separate(bplot_SAASD, show = T, image.plot_extra_options = list('horizontal' = T, 'legend.width' = 2, 'legend.mar' = 14.1, 'legend.only' = T, 'col' = colFn_raw('n' = 11L), 'legend.lab' = 'Mean Asymmetry Index (L > R)', 'legend.cex' = 1.2, 'legend.line' = 3))

# Exporting brainplots as PNG
export(bplot_CTHC, draw_colorbar = F, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_HC.png')
export(bplot_CTASD, draw_colorbar = F, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_ASD.png')
export(bplot_CTASD1, draw_colorbar = F, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_ASD1.png')
export(bplot_CTASD2, draw_colorbar = F, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_ASD2.png')
export(bplot_CTASD3, draw_colorbar = T, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_ASD3.png')

export(bplot_SAHC, draw_colorbar = F, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_HC.png')
export(bplot_SAASD, draw_colorbar = F, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_ASD.png')
export(bplot_SAASD1, draw_colorbar = F, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_ASD1.png')
export(bplot_SAASD2, draw_colorbar = F, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_ASD2.png')
export(bplot_SAASD3, draw_colorbar = T, colorbar_legend = 'Mean Asymmetry Indexes (L>R)', view_angles = c('sd_lateral_lh', 'sd_medial_lh'), output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_ASD3.png')

# HC
paths1 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_HC.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_HC.png')
arrange.brainview.images.grid(num_per_row = 2, brainview_images = paths1, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/HC.png')

# ASD
paths2 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_ASD.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_ASD.png')
arrange.brainview.images.grid(num_per_row = 2, brainview_images = paths2, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/ASD.png')

# ASD1
paths3 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_ASD1.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_ASD1.png')
arrange.brainview.images.grid(num_per_row = 2, brainview_images = paths3, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/ASD1.png')

# ASD2
paths4 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_ASD2.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_ASD2.png')
arrange.brainview.images.grid(num_per_row = 2, brainview_images = paths4, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/ASD2.png')

# ASD3
paths5 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/CT_ASD3.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/SA_ASD3.png')
arrange.brainview.images.grid(num_per_row = 2, brainview_images = paths5, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/raw_AI/ASD3.png')



# VISUALISATION OF SIGNIFICANT BRAIN REGIONS

# Setting color function options
mkcCT <- list('colFn' = colFn_raw, 'range' = c(-.04, .04), 'symm' = T)
mkcSA <- list('colFn' = colFn_raw, 'range' = c(-.3, .3), 'symm' = T)

#creating arbitrary regionlist for rh (will not be diplayed)
emptyRegionList <- list('G_postcentral' = 0.00001)

## Cortical Thickness

# ASD v HC 
ttestCT_AvCsig
# --> no significant regions

# ASD1 v HC
ttestCT_1vC$sigTable
# --> no significant regions

# ASD2 v HC 
ttestCT_2vC$sigTable
# --> no significant regions

# ASD1 v ASD2
# ASD1
ttestCT_1v2$sigTable
regions_CT_ASD1v2_1 <- c(ttestCT_1v2$sigTable$MeanX)
names(regions_CT_ASD1v2_1) <- c(ttestCT_1v2$sigTable$FSregion)

bplot_CTASD1v2_1 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                subject_id = template_subject, 
                                                atlas = params$atlas, 
                                                lh_region_value_list = regions_CT_ASD1v2_1, 
                                                rh_region_value_list = emptyRegionList, 
                                                surface = 'white', 
                                                draw_colorbar = T, 
                                                makecmap_options = mkcCT, 
                                                views = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'),
                                                rgloptions = rgl_options,
                                                bg= 'curv')
# ASD2
ttestCT_1v2$sigTable
regions_CT_ASD1v2_2 <- c(ttestCT_1v2$sigTable$MeanY)
names(regions_CT_ASD1v2_2) <- c(ttestCT_1v2$sigTable$FSregion)

bplot_CTASD1v2_2 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                 subject_id = template_subject, 
                                                 atlas = params$atlas, 
                                                 lh_region_value_list = regions_CT_ASD1v2_2, 
                                                 rh_region_value_list = emptyRegionList, 
                                                 surface = 'white', 
                                                 draw_colorbar = T, 
                                                 makecmap_options = mkcCT, 
                                                 views = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'),
                                                 rgloptions = rgl_options,
                                                 bg= 'curv')

## Surface Area

# ASD v HC
ttestSA_AvCsig
# --> no significant regions

# ASD1 v HC
ttestSA_1vC$sigTable
regions_SA_ASD1vHC <- c(ttestSA_1vC$sigTable$MeanY)
names(regions_SA_ASD1vHC) <- c(ttestSA_1vC$sigTable$FSregion)

bplot_SAASD1vHC <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                subject_id = template_subject, 
                                                atlas = params$atlas, 
                                                lh_region_value_list = regions_SA_ASD1vHC, 
                                                rh_region_value_list = emptyRegionList, 
                                                surface = 'white', 
                                                draw_colorbar = T, 
                                                makecmap_options = mkcSA, 
                                                views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                rgloptions = rgl_options,
                                                bg= 'curv')

# ASD2 v HC
ttestSA_2vC$sigTable
regions_SA_ASD2vHC <- c(ttestSA_2vC$sigTable$MeanY)
names(regions_SA_ASD2vHC) <- c(ttestSA_2vC$sigTable$FSregion)

bplot_SAASD2vHC <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                subject_id = template_subject, 
                                                atlas = params$atlas, 
                                                lh_region_value_list = regions_SA_ASD2vHC, 
                                                rh_region_value_list = emptyRegionList, 
                                                surface = 'white', 
                                                draw_colorbar = T, 
                                                makecmap_options = mkcSA, 
                                                views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                rgloptions = rgl_options,
                                                bg= 'curv')

# ASD1/2 v HC (plot for HC)
ttestSA_2vC$sigTable
regions_SA_HC <- c(ttestSA_2vC$sigTable$MeanX)
names(regions_SA_HC) <- c(ttestSA_2vC$sigTable$FSregion)

bplot_SAHCsig <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                subject_id = template_subject, 
                                                atlas = params$atlas, 
                                                lh_region_value_list = regions_SA_HC, 
                                                rh_region_value_list = emptyRegionList, 
                                                surface = 'white', 
                                                draw_colorbar = T, 
                                                makecmap_options = mkcSA, 
                                                views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                rgloptions = rgl_options,
                                                bg= 'curv')

# ASD1 v ASD2
# ASD1
ttestSA_1v2$sigTable
regions_SA_ASD1v2_1 <- c(ttestSA_1v2$sigTable$MeanX)
names(regions_SA_ASD1v2_1) <- c(ttestSA_1v2$sigTable$FSregion)

bplot_SAASD1v2_1 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                 subject_id = template_subject, 
                                                 atlas = params$atlas, 
                                                 lh_region_value_list = regions_SA_ASD1v2_1, 
                                                 rh_region_value_list = emptyRegionList, 
                                                 surface = 'white', 
                                                 draw_colorbar = T, 
                                                 makecmap_options = mkcSA, 
                                                 views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                 rgloptions = rgl_options,
                                                 bg= 'curv')
# ASD2
ttestSA_1v2$sigTable
regions_SA_ASD1v2_2 <- c(ttestSA_1v2$sigTable$MeanY)
names(regions_SA_ASD1v2_2) <- c(ttestSA_1v2$sigTable$FSregion)

bplot_SAASD1v2_2 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                 subject_id = template_subject, 
                                                 atlas = params$atlas, 
                                                 lh_region_value_list = regions_SA_ASD1v2_2, 
                                                 rh_region_value_list = emptyRegionList, 
                                                 surface = 'white', 
                                                 draw_colorbar = T, 
                                                 makecmap_options = mkcSA, 
                                                 views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                 rgloptions = rgl_options,
                                                 bg= 'curv')

# separate legend
legendMeanCT <- coloredmesh.plot.colorbar.separate(bplot_CTASD1v2_1, show = T, image.plot_extra_options = list('horizontal' = T, 'legend.width' = 2, 'legend.mar' = 14.1, 'legend.only' = T, 'col' = colFn_raw('n' = 10L), 'legend.lab' = 'Mean Asymmetry Index (L > R)', 'legend.cex' = 1.1, 'legend.line' = 2.5), png_options = list(filename = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/cbar_CT.png', res = 100))
legendMeanSA <- coloredmesh.plot.colorbar.separate(bplot_SAASD1v2_1, show = T, image.plot_extra_options = list('horizontal' = T, 'legend.width' = 2, 'legend.mar' = 14.1, 'legend.only' = T, 'col' = colFn_raw('n' = 10L), 'legend.lab' = 'Mean Asymmetry Index (L > R)', 'legend.cex' = 1.1, 'legend.line' = 2.5), png_options = list(filename = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/cbar_SA.png', res = 100))
legendMeanSAvert <- coloredmesh.plot.colorbar.separate(bplot_SAASD1v2_1, show = T, image.plot_extra_options = list('horizontal' = F, 'legend.width' = 2, 'legend.mar' = 14.1, 'legend.only' = T, 'col' = colFn_raw('n' = 10L), 'legend.lab' = 'Mean Asymmetry Index (L > R)', 'legend.cex' = 1.1, 'legend.line' = 3), png_options = list(filename = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/cbar_SA_v.png', res = 100))

# Exporting brainplots as PNG
export(bplot_CTASD1v2_1, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/CT_ASD1v2_1.png')
export(bplot_CTASD1v2_2, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/CT_ASD1v2_2.png')

export(bplot_SAASD1vHC, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_ASD1vHC.png')
export(bplot_SAASD2vHC, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_ASD2vHC.png')
export(bplot_SAASD1v2_1, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_ASD1v2_1.png')
export(bplot_SAASD1v2_2, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_ASD1v2_2.png')

export(bplot_SAHCsig, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_HCsig.png')

# Comparison: ASD1/2 v C
paths1 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_HCsig.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_ASD1vHC.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_ASD2vHC.png')
arrange.brainview.images.grid(num_per_row = 1, brainview_images = paths1, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/ASDvCsig.png')

# Comparison: ASD1 v ASD2
paths2 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/CT_ASD1v2_1.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_ASD1v2_1.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/CT_ASD1v2_2.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/SA_ASD1v2_2.png')
arrange.brainview.images.grid(num_per_row = 2, brainview_images = paths2,  output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Significant/ASD1v2sig.png')



### VISUALISATION OF DIFFERENCE MEASURE

## Cortical Thickness
mkcCT_dif_abs <- list('colFn' = colFn_raw, 'range' = c(-.05, .05), 'symm' = T)

# CT ASD1 v HC (HC - ASD1)
ttestCT_1vC$sigTable
# --> no significant regions

# CT ASD1 v HC (HC - ASD2)
ttestCT_2vC$sigTable
# --> no significant regions

# CT ASD1 v ASD2 (ASD1 - ASD2)
ttestCT_1v2$sigTable
regionsCT_dif_1v2_1 <- c(abs(ttestCT_1v2$sigTable$MeanX) - abs(ttestCT_1v2$sigTable$MeanY))
names(regionsCT_dif_1v2_1) <- c(ttestCT_1v2$sigTable$FSregion)

bplotCT_dif_1v2_1 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                subject_id = template_subject, 
                                                atlas = params$atlas, 
                                                lh_region_value_list = regionsCT_dif_1v2_1, 
                                                rh_region_value_list = emptyRegionList, 
                                                surface = 'white', 
                                                draw_colorbar = T, 
                                                makecmap_options = mkcCT_dif_abs, 
                                                views = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'),
                                                rgloptions = rgl_options,
                                                bg = 'curv')

# CT ASD1 v ASD2 (ASD2 - ASD1)
ttestCT_1v2$sigTable
regionsCT_dif_1v2_2 <- c(abs(ttestCT_1v2$sigTable$MeanY) - abs(ttestCT_1v2$sigTable$MeanX))
names(regionsCT_dif_1v2_2) <- c(ttestCT_1v2$sigTable$FSregion)

bplotCT_dif_1v2_2 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                subject_id = template_subject, 
                                                atlas = params$atlas, 
                                                lh_region_value_list = regionsCT_dif_1v2_2, 
                                                rh_region_value_list = emptyRegionList, 
                                                surface = 'white', 
                                                draw_colorbar = T, 
                                                makecmap_options = mkcCT_dif_abs, 
                                                views = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'),
                                                rgloptions = rgl_options,
                                                bg = 'curv')

## Surface Area
mkcSA_dif_abs <- list('colFn' = colFn_raw, 'range' = c(-.2, .2), 'symm' = T)

# SA ASD-I v HC (ASD1 - HC)
ttestSA_1vC$sigTable
regionsSA_dif_1vHC <- c(abs(ttestSA_1vC$sigTable$MeanY) - abs(ttestSA_1vC$sigTable$MeanX))
names(regionsSA_dif_1vHC) <- c(ttestSA_1vC$sigTable$FSregion)

bplotSA_dif_1vHC <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                 subject_id = template_subject, 
                                                 atlas = params$atlas, 
                                                 lh_region_value_list = regionsSA_dif_1vHC, 
                                                 rh_region_value_list = emptyRegionList, 
                                                 surface = 'white', 
                                                 draw_colorbar = T, 
                                                 makecmap_options = mkcSA_dif_abs, 
                                                 views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                 rgloptions = rgl_options,
                                                 bg = 'curv')

# SA ASD-II v HC (ASD2 - HC)
ttestSA_2vC$sigTable
regionsSA_dif_2vHC <- c(abs(ttestSA_2vC$sigTable$MeanY) - abs(ttestSA_2vC$sigTable$MeanX))
names(regionsSA_dif_2vHC) <- c(ttestSA_2vC$sigTable$FSregion)

bplotSA_dif_2vHC <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                 subject_id = template_subject, 
                                                 atlas = params$atlas, 
                                                 lh_region_value_list = regionsSA_dif_2vHC, 
                                                 rh_region_value_list = emptyRegionList, 
                                                 surface = 'white', 
                                                 draw_colorbar = T, 
                                                 makecmap_options = mkcSA_dif_abs, 
                                                 views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                 rgloptions = rgl_options,
                                                 bg = 'curv')

# SA ASD-I v ASD-II (ASD1 - ASD2)
ttestSA_1v2$sigTable
regionsSA_dif_1v2_1 <- c(abs(ttestSA_1v2$sigTable$MeanX) - abs(ttestSA_1v2$sigTable$MeanY))
names(regionsSA_dif_1v2_1) <- c(ttestSA_1v2$sigTable$FSregion)

bplotSA_dif_1v2_1 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                subject_id = template_subject, 
                                                atlas = params$atlas, 
                                                lh_region_value_list = regionsSA_dif_1v2_1, 
                                                rh_region_value_list = emptyRegionList, 
                                                surface = 'white', 
                                                draw_colorbar = T, 
                                                makecmap_options = mkcSA_dif_abs, 
                                                views = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'),
                                                rgloptions = rgl_options,
                                                bg = 'curv')

# SA ASD-I v ASD-II (ASD2 - ASD1)
ttestSA_1v2$sigTable
regionsSA_dif_1v2_2 <- c(abs(ttestSA_1v2$sigTable$MeanY) - abs(ttestSA_1v2$sigTable$MeanX))
names(regionsSA_dif_1v2_2) <- c(ttestSA_1v2$sigTable$FSregion)

bplotSA_dif_1v2_2 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                  subject_id = template_subject, 
                                                  atlas = params$atlas, 
                                                  lh_region_value_list = regionsSA_dif_1v2_2, 
                                                  rh_region_value_list = emptyRegionList, 
                                                  surface = 'white', 
                                                  draw_colorbar = T, 
                                                  makecmap_options = mkcSA_dif_abs, 
                                                  views = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'),
                                                  rgloptions = rgl_options,
                                                  bg = 'curv')

# separate legend
legendDiffCT <- coloredmesh.plot.colorbar.separate(bplotCT_dif_1v2_1, png_options = list(filename = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/cbar_CT.png', res = 100), show = T, image.plot_extra_options = list('legend.width' = 2, 'legend.mar' = 14.1, 'legend.only' = T, 'col' = colFn_raw('n' = 10L), 'legend.lab' = 'Absolute Difference', 'legend.cex' = 1.1, 'legend.line' = 4))
legendDiffSA <- coloredmesh.plot.colorbar.separate(bplotSA_dif_1v2_1, png_options = list(filename = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/cbar_SA.png', res = 100), show = T, image.plot_extra_options = list('legend.width' = 2, 'legend.mar' = 14.1, 'legend.only' = T, 'col' = colFn_raw('n' = 10L), 'legend.lab' = 'Absolute Difference', 'legend.cex' = 1.1, 'legend.line' = 4))

# Exporting brainplots as PNG
export(bplotCT_dif_1v2_1, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/CT_ASD1v2_1.png')
export(bplotCT_dif_1v2_2, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/CT_ASD1v2_2.png')

export(bplotSA_dif_1vHC, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/SA_ASD1vHC.png')
export(bplotSA_dif_2vHC, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/SA_ASD2vHC.png')
export(bplotSA_dif_1v2_1, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/SA_ASD1v2_1.png')
export(bplotSA_dif_1v2_2, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/SA_ASD1v2_2.png')

# Comparison: ASD1/2 v C
paths1 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/SA_ASD1vHC.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/SA_ASD2vHC.png')
arrange.brainview.images.grid(num_per_row = 1, brainview_images = paths1, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/ASDvC.png')

# Comparison: ASD1 v ASD2
paths2 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/CT_ASD1v2_1.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/SA_ASD1v2_1.png')
arrange.brainview.images.grid(num_per_row = 1, brainview_images = paths2, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/Difference/ASD1v2.png')


### VISUALISATION OF COHENS D

mkc_cohenCT <- list('colFn' = colFn_raw, 'range' = c(-1, 1), 'symm' = T)
mkc_cohenSA <- list('colFn' = colFn_raw, 'range' = c(-2, 2), 'symm' = T)

## Cortical Thickness

# ASD1 v HC
ttestCT_1vC$sigTable
# --> no significant regions

# ASD2 v HC
ttestCT_2vC$sigTable
# --> no significant regions

# ASD1 v ASD2
ttestCT_1v2$sigTable
regions_cohen_CT1v2 <- c(ttestCT_1v2$sigTable$CohenD)
names(regions_cohen_CT1v2) <- c(ttestCT_1v2$sigTable$FSregion)

bplotCT_cohen1v2 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                 subject_id = template_subject, 
                                                 atlas = params$atlas,
                                                 lh_region_value_list = regions_cohen_CT1v2, 
                                                 rh_region_value_list = emptyRegionList, 
                                                 surface = 'white', 
                                                 draw_colorbar = T, 
                                                 makecmap_options = mkc_cohenCT, 
                                                 views = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'),
                                                 rgloptions = rgl_options,
                                                 bg = 'curv')

## Surface Area
# ASD1 v HC
ttestSA_1vC$sigTable
regions_cohen_SA1vC <- c(ttestSA_1vC$sigTable$CohenD) 
names(regions_cohen_SA1vC) <- c(ttestSA_1vC$sigTable$FSregion)

bplotSA_cohen1vC <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                 subject_id = template_subject, 
                                                 atlas = params$atlas,
                                                 lh_region_value_list = regions_cohen_SA1vC, 
                                                 rh_region_value_list = emptyRegionList, 
                                                 surface = 'white', 
                                                 draw_colorbar = T, 
                                                 makecmap_options = mkc_cohenSA, 
                                                 views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                 rgloptions = rgl_options,
                                                 bg = 'curv')

# ASD2 v HC
ttestSA_2vC$sigTable
regions_cohen_SA2vC <- c(ttestSA_2vC$sigTable$CohenD)
names(regions_cohen_SA2vC) <- c(ttestSA_2vC$sigTable$FSregion)

bplotSA_cohen2vC <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                 subject_id = template_subject, 
                                                 atlas = params$atlas,
                                                 lh_region_value_list = regions_cohen_SA2vC, 
                                                 rh_region_value_list = emptyRegionList, 
                                                 surface = 'white', 
                                                 draw_colorbar = T, 
                                                 makecmap_options = mkc_cohenSA, 
                                                 views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                 rgloptions = rgl_options,
                                                 bg = 'curv')

# ASD1 v ASD2
ttestSA_1v2$sigTable
regions_cohen_SA1v2 <- c(ttestSA_1v2$sigTable$CohenD)
names(regions_cohen_SA1v2) <- c(ttestSA_1v2$sigTable$FSregion)

bplotSA_cohen1v2 <- vis.region.values.on.subject(subjects_dir = template_subject_dir, 
                                                 subject_id = template_subject, 
                                                 atlas = params$atlas,
                                                 lh_region_value_list = regions_cohen_SA1v2, 
                                                 rh_region_value_list = emptyRegionList, 
                                                 surface = 'white', 
                                                 draw_colorbar = T, 
                                                 makecmap_options = mkc_cohenSA, 
                                                 views = c('sd_lateral_lh', 'sd_medial_lh'),
                                                 rgloptions = rgl_options,
                                                 bg = 'curv')

# separate legend
legendCohenCT <- coloredmesh.plot.colorbar.separate(bplotCT_cohen1v2, show = T, png_options = list(filename = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/cbar_CT.png', res = 100), image.plot_extra_options = list('legend.width' = 2, 'legend.mar' = 14.1, 'legend.only' = T, 'col' = colFn_raw('n' = 10L), 'legend.lab' = 'Cohen`s D', 'legend.cex' = 1.1, 'legend.line' = 4))
legendCohenSA <- coloredmesh.plot.colorbar.separate(bplotSA_cohen1v2, show = T, png_options = list(filename = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/cbar_SA.png', res = 100), image.plot_extra_options = list('legend.width' = 2, 'legend.mar' = 14.1, 'legend.only' = T, 'col' = colFn_raw('n' = 10L), 'legend.lab' = 'Cohen`s D', 'legend.cex' = 1.1, 'legend.line' = 4))

# Exporting brainplots as PNG
export(bplotCT_cohen1v2, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/CT_ASD1v2.png')

export(bplotSA_cohen1vC, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/SA_ASD1vC.png')
export(bplotSA_cohen2vC, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/SA_ASD2vC.png')
export(bplotSA_cohen1v2, draw_colorbar = F, view_angles = c('sd_lateral_lh', 'sd_medial_lh', 'sd_ventral'), grid_like = F, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/SA_ASD1v2.png')

# Comparison: ASD1/2 v C
paths1 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/SA_ASD1vC.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/SA_ASD2vC.png')
arrange.brainview.images.grid(num_per_row = 1, brainview_images = paths1, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/ASDvC.png')

# Comparison: ASD1 v ASD2
paths2 = c('~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/CT_ASD1v2.png', '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/SA_ASD1v2.png')
arrange.brainview.images.grid(num_per_row = 1, brainview_images = paths2, output_img = '~/Documents/Psychologie/MA/Results/Anatomical/PNG/CohensD/ASD1v2.png')


## Preliminary exploration of ASD3

# ASD3 v C
ttestCT_3vC <- ttest_brain_function(allCT, 0, 3, 1)
ttestCT_3vC$sigTable

ttestSA_3vC <- ttest_brain_function(allSA, 0, 3, nyholt.bonf_SA)
ttestSA_3vC$sigTable

# ASD3 v ASD1
ttestCT_3v1 <- ttest_brain_function(allCT, 1, 3, 1)
ttestCT_3v1$sigTable

ttestSA_3v1 <- ttest_brain_function(allSA, 1, 3, 1)
ttestCT_3v1$sigTable

# ASD3 v ASD2
ttestCT_3vC <- ttest_brain_function(allCT, 0, 3, 1)
ttestCT_3vC$sigTable

ttestSA_3v2 <- ttest_brain_function(allSA, 0, 3, 1)
ttestSA_3v2$sigTable

