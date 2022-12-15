## Combining stats files

source('~/Documents/Psychologie/MA/Code/R/GeneralParameters.R')

setwd(params$PathBrainData)

# ASD subjects
ASD_lh_area_cb <- read.delim("./ASD/lh_area_a2009s_stats_cbic.txt", header = T, sep = "\t", dec = ".")
ASD_lh_area_cb$site <- 'cb'
ASD_lh_area_cu <- read.delim("./ASD/lh_area_a2009s_stats_cuny.txt", header = T, sep = "\t", dec = ".")
ASD_lh_area_cu$site <- 'cu'
ASD_lh_area_ru <- read.delim("./ASD/lh_area_a2009s_stats_ru.txt", header = T, sep = "\t", dec = ".")
ASD_lh_area_ru$site <- 'ru'
ASD_lh_area <- rbind(ASD_lh_area_cb, ASD_lh_area_cu, ASD_lh_area_ru)
names(ASD_lh_area)[names(ASD_lh_area) == "lh.aparc.a2009s.area"] <- "ID"


ASD_rh_area_cb <- read.delim("./ASD/rh_area_a2009s_stats_cbic.txt", header = T, sep = "\t", dec = ".")
ASD_rh_area_cb$site <- 'cb'
ASD_rh_area_cu <- read.delim("./ASD/rh_area_a2009s_stats_cuny.txt", header = T, sep = "\t", dec = ".")
ASD_rh_area_cu$site <- 'cu'
ASD_rh_area_ru <- read.delim("./ASD/rh_area_a2009s_stats_ru.txt", header = T, sep = "\t", dec = ".")
ASD_rh_area_ru$site <- 'ru'
ASD_rh_area <- rbind(ASD_rh_area_cb, ASD_rh_area_cu, ASD_rh_area_ru)


ASD_lh_thickness_cb <- read.delim("./ASD/lh_thickness_a2009s_stats_cbic.txt", header = T, sep = "\t", dec = ".")
ASD_lh_thickness_cb$site <- 'cb'
ASD_lh_thickness_cu <- read.delim("./ASD/lh_thickness_a2009s_stats_cuny.txt", header = T, sep = "\t", dec = ".")
ASD_lh_thickness_cu$site <- 'cu'
ASD_lh_thickness_ru <- read.delim("./ASD/lh_thickness_a2009s_stats_ru.txt", header = T, sep = "\t", dec = ".")
ASD_lh_thickness_ru$site <- 'ru'
ASD_lh_thickness <- rbind(ASD_lh_thickness_cb, ASD_lh_thickness_cu, ASD_lh_thickness_ru)


ASD_rh_thickness_cb <- read.delim("./ASD/rh_thickness_a2009s_stats_cbic.txt", header = T, sep = "\t", dec = ".")
ASD_rh_thickness_cb$site <- 'cb'
ASD_rh_thickness_cu <- read.delim("./ASD/rh_thickness_a2009s_stats_cuny.txt", header = T, sep = "\t", dec = ".")
ASD_rh_thickness_cu$site <- 'cu'
ASD_rh_thickness_ru <- read.delim("./ASD/rh_thickness_a2009s_stats_ru.txt", header = T, sep = "\t", dec = ".")
ASD_rh_thickness_ru$site <- 'ru'
ASD_rh_thickness <- rbind(ASD_rh_thickness_cb, ASD_rh_thickness_cu, ASD_rh_thickness_ru)


# HC subjects
HC_lh_area_cb <- read.delim("./HC/lh_area_a2009s_stats_cbic.txt", header = T, sep = "\t", dec = ".")
HC_lh_area_cb$site <- 'cb'
HC_lh_area_cu <- read.delim("./HC/lh_area_a2009s_stats_cuny.txt", header = T, sep = "\t", dec = ".")
HC_lh_area_cu$site <- 'cu'
HC_lh_area_ru <- read.delim("./HC/lh_area_a2009s_stats_ru.txt", header = T, sep = "\t", dec = ".")
HC_lh_area_ru$site <- 'ru'
HC_lh_area <- rbind(HC_lh_area_cb, HC_lh_area_cu, HC_lh_area_ru)


HC_rh_area_cb <- read.delim("./HC/rh_area_a2009s_stats_cbic.txt", header = T, sep = "\t", dec = ".")
HC_rh_area_cb$site <- 'cb'
HC_rh_area_cu <- read.delim("./HC/rh_area_a2009s_stats_cuny.txt", header = T, sep = "\t", dec = ".")
HC_rh_area_cu$site <- 'cu'
HC_rh_area_ru <- read.delim("./HC/rh_area_a2009s_stats_ru.txt", header = T, sep = "\t", dec = ".")
HC_rh_area_ru$site <- 'ru'
HC_rh_area <- rbind(HC_rh_area_cb, HC_rh_area_cu, HC_rh_area_ru)


HC_lh_thickness_cb <- read.delim("./HC/lh_thickness_a2009s_stats_cbic.txt", header = T, sep = "\t", dec = ".")
HC_lh_thickness_cb$site <- 'cb'
HC_lh_thickness_cu <- read.delim("./HC/lh_thickness_a2009s_stats_cuny.txt", header = T, sep = "\t", dec = ".")
HC_lh_thickness_cu$site <- 'cu'
HC_lh_thickness_ru <- read.delim("./HC/lh_thickness_a2009s_stats_ru.txt", header = T, sep = "\t", dec = ".")
HC_lh_thickness_ru$site <- 'ru'
HC_lh_thickness <- rbind(HC_lh_thickness_cb, HC_lh_thickness_cu, HC_lh_thickness_ru)


HC_rh_thickness_cb <- read.delim("./HC/rh_thickness_a2009s_stats_cbic.txt", header = T, sep = "\t", dec = ".")
HC_rh_thickness_cb$site <- 'cb'
HC_rh_thickness_cu <- read.delim("./HC/rh_thickness_a2009s_stats_cuny.txt", header = T, sep = "\t", dec = ".")
HC_rh_thickness_cu$site <- 'cu'
HC_rh_thickness_ru <- read.delim("./HC/rh_thickness_a2009s_stats_ru.txt", header = T, sep = "\t", dec = ".")
HC_rh_thickness_ru$site <- 'ru'
HC_rh_thickness <- rbind(HC_rh_thickness_cb, HC_rh_thickness_cu, HC_rh_thickness_ru)

# Writing stats files
write.table(ASD_lh_area, file = './ASD/ASD_lh_area.txt', quote = F, sep = "\t", row.names = F, col.names = T)
write.table(ASD_rh_area, file = './ASD/ASD_rh_area.txt', quote = F, sep = "\t", row.names = F, col.names = T)
write.table(ASD_lh_thickness, file = './ASD/ASD_lh_thickness.txt', quote = F, sep = "\t", row.names = F, col.names = T)
write.table(ASD_rh_thickness, file = './ASD/ASD_rh_thickness.txt', quote = F, sep = "\t", row.names = F, col.names = T)

write.table(HC_lh_area, file = './HC/HC_lh_area.txt', quote = F, sep = '\t', row.names = F, col.names = T)
write.table(HC_rh_area, file = './HC/HC_rh_area.txt', quote = F, sep = '\t', row.names = F, col.names = T)
write.table(HC_lh_thickness, file = './HC/HC_lh_thickness.txt', quote = F, sep = '\t', row.names = F, col.names = T)
write.table(HC_rh_thickness, file = './HC/HC_rh_thickness.txt', quote = F, sep = '\t', row.names = F, col.names = T)

