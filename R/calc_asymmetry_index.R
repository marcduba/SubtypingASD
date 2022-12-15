###################################################
# Calculating Asymmetry Index and Visualization  #
##################################################

library(dplyr)
library(fsbrain)
source('~/Documents/Psychologie/MA/Code/R/GeneralParameters.R')

setwd(params$FreeSurferDir)
template_subject = params$template_sub
template_subject_dir = params$template_sub_dir

atlas_region_names = get.atlas.region.names(atlas = params$atlas, template_subject = template_subject, hemi = 'lh')

regionnames <- atlas_region_names[atlas_region_names != 'Unknown']
regionnames <- regionnames[regionnames != 'Medial_wall']

setwd(params$PathBrainData)

# Loading files
ASD_lh_area <- read.delim('./ASD/ASD_lh_area.txt', header = T, sep = '\t', dec = '.')
ASD_rh_area <- read.delim('./ASD/ASD_rh_area.txt', header = T, sep = '\t', dec = '.')
ASD_lh_thickness <- read.delim('./ASD/ASD_lh_thickness.txt', header = T, sep = '\t', dec = '.')
ASD_rh_thickness <- read.delim('./ASD/ASD_rh_thickness.txt', header = T, sep = '\t', dec = '.')

names(ASD_lh_area)[names(ASD_lh_area) == "lh.aparc.a2009s.area"] <- "ID"
names(ASD_rh_area)[names(ASD_rh_area) == "rh.aparc.a2009s.area"] <- "ID"
names(ASD_lh_thickness)[names(ASD_lh_thickness) == "lh.aparc.a2009s.thickness"] <- "ID"
names(ASD_rh_thickness)[names(ASD_rh_thickness) == "rh.aparc.a2009s.thickness"] <- "ID"

rownames(ASD_lh_area) <- c(ASD_lh_area$ID)
rownames(ASD_rh_area) <- c(ASD_rh_area$ID)
rownames(ASD_lh_thickness) <- c(ASD_lh_thickness$ID)
rownames(ASD_rh_thickness) <- c(ASD_rh_thickness$ID)

HC_lh_area <- read.delim('./HC/HC_lh_area.txt', header = T, sep = '\t', dec = '.')
HC_rh_area <- read.delim('./HC/HC_rh_area.txt', header = T, sep = '\t', dec = '.')
HC_lh_thickness <- read.delim('./HC/HC_lh_thickness.txt', header = T, sep = '\t', dec = '.')
HC_rh_thickness <- read.delim('./HC/HC_rh_thickness.txt', header = T, sep = '\t', dec = '.')

names(HC_lh_area)[names(HC_lh_area) == "lh.aparc.a2009s.area"] <- "ID"
names(HC_rh_area)[names(HC_rh_area) == "rh.aparc.a2009s.area"] <- "ID"
names(HC_lh_thickness)[names(HC_lh_thickness) == "lh.aparc.a2009s.thickness"] <- "ID"
names(HC_rh_thickness)[names(HC_rh_thickness) == "rh.aparc.a2009s.thickness"] <- "ID"

rownames(HC_lh_area) <- c(HC_lh_area$ID)
rownames(HC_rh_area) <- c(HC_rh_area$ID)
rownames(HC_lh_thickness) <- c(HC_lh_thickness$ID)
rownames(HC_rh_thickness) <- c(HC_rh_thickness$ID)

# Formula for Asymmetry Index
#AI_SA = (SA_LH - SA_RH)/(SA_LH + SA_RH)
#AI_CT = (CT_LH - CT_RH)/(CT_LH + CT_RH)

AI_function <- function(dflh, dfrh) {
  sub = data.frame((dflh[, 2:75]) - (dfrh[, 2:75]))
  sum = data.frame((dflh[, 2:75]) + (dfrh[, 2:75]))
  AI = data.frame(sub / sum)
  return(AI)
}

# Manually checking if function produces correct output
ASD_SA_sub = data.frame((ASD_lh_area[ , 2:75]) - (ASD_rh_area[ , 2:75]))
ASD_SA_sum = data.frame((ASD_lh_area[ , 2:75]) + (ASD_rh_area[ , 2:75]))
ASD_AI_SA = ASD_SA_sub / ASD_SA_sum

# ASD SA
AI_ASD_SA <- AI_function(ASD_lh_area, ASD_rh_area)

colnames(AI_ASD_SA) <- regionnames

# ASD CT
AI_ASD_CT <- AI_function(ASD_lh_thickness, ASD_rh_thickness)

colnames(AI_ASD_CT) <- regionnames

# HC SA
AI_HC_SA <- AI_function(HC_lh_area, HC_rh_area)

colnames(AI_HC_SA) <- regionnames

# HC CT
AI_HC_CT <- AI_function(HC_lh_thickness, HC_rh_thickness)

colnames(AI_HC_CT) <- regionnames

# Writing stats files
write.table(AI_ASD_SA, file = './ASD/AI_ASD_SA.txt', quote = F, sep = '\t', row.names = T, col.names = T)
write.table(AI_ASD_CT, file = './ASD/AI_ASD_CT.txt', quote = F, sep = '\t', row.names = T, col.names = T)
write.table(AI_HC_SA, file = './HC/AI_HC_SA.txt', quote = F, sep = '\t', row.names = T, col.names = T)
write.table(AI_HC_CT, file = './HC/AI_HC_CT.txt', quote = F, sep = '\t', row.names = T, col.names = T)
