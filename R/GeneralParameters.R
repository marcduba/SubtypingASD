######################
# General Parameters #
######################

# Definition of general parameters for data analysis

params <- list(
  # Paths to data and files
  WorkingDirectory = '~/Documents/Psychologie/MA/',
  BasicDataPath = '/Volumes/methlab_data/HBN/BD.csv',
  FreeSurferDir = '/Applications/freesurfer/7.2.0',
  PathBrainData = '~/Documents/Psychologie/MA/Data/',

  # Subjects
  ASDdiag = 'Autism Spectrum Disorder',
  HCdiag = 'No Diagnosis Given',
  EHQCutoff = 50,
  
  # Behavioral Measures
  BehavMeasures = c('ASSQ.ASSQ_Total', 'SRS.SRS_Total', 'RBS.RBS_Total', 'SCQ.SCQ_Total', 'Basic_Demos.Age', 'Basic_Demos.Sex'), 
  SRSsubscales = c('SRS.SRS_AWR', 'SRS.SRS_COG', 'SRS.SRS_COM', 'SRS.SRS_DSMRRB', 'SRS.SRS_MOT', 'SRS.SRS_RRB', 'SRS.SRS_SCI'),
  paletteASD = c('#00AFBB', '#E7B800', '#FC4E07'),
  
  # MRI parameters
  Protocols = list(HCP = 'HCP', VNav = 'VNav', VNavN = 'VNavNorm', Other = 'Unknown'),
  template_sub_dir = '/Users/marcdubacher/Documents/Psychologie/MA/Data',
  template_sub = 'fsaverage',
  atlas = 'aparc.a2009s',
  
  # Clustering parameters
  Distance = 'euclidean',
  ClustMethod = 'ward.D2', 
  
  
  # Function to split IDs
  # adapted from Andreas Papageorgiou
  IDsplitFunction <- function(RunID) {
    splt1 <- strsplit(RunID, "_")
    EID <- strsplit(splt1[[1]][1],"-")[[1]][2]
    Protocol <- strsplit(splt1[[1]][2],"-")[[1]][2]
    if (is.na(Protocol)) {
      Protocol = params$Protocols$Other
    }
    Run <- 1
    if (length(splt1[[1]]) == 4) {
      Run <- str_sub(splt1[[1]][3], start = -1)
    }
    
    return(c(EID, Protocol, Run))
  },
  
  # Nyholt measure for multiple comparisons
  # supplied by Sabine Dziemian/Andreas Papageorgiou
  # low Meff indicates high linkage disequilibrium, high Meff indicates low linkage disequilibrium
  getMeff <- function(D){
    DT <- na.omit(D)
    DT <- as.matrix(DT)
    
    eigen <- eigen(cor(DT))
    v.eigen <- var(eigen$values)
    
    M <- ncol(D)
    Meff <- 1 + (M - 1) * (1 - (v.eigen/M))
    
    return(Meff)
    
  },
  
  # Function to compute Nyholt & Bonferroni correction
  nyholt.bonf.p <- function(data, p) {
    Meff <- getMeff(data)
    return(p / Meff)
  }
)




