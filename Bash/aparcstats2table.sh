
## Setting Up FreeSurfer
export FREESURFER_HOME=/Applications/freesurfer/7.2.0
source $FREESURFER_HOME/SetUpFreeSurfer.sh


## Extracting stats files for lh/rh area and thickness

##ASD

# CBIC:
export SUBJECTS_DIR=/Volumes/methlab_data/HBN/MRI/Site-CBIC_Derivatives_UZH/

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_CB.txt --hemi rh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/rh_area_a2009s_stats_cbic.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_CB.txt --hemi lh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/lh_area_a2009s_stats_cbic.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_CB.txt --hemi rh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/rh_thickness_a2009s_stats_cbic.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_CB.txt --hemi lh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/lh_thickness_a2009s_stats_cbic.txt


# CUNY:
export SUBJECTS_DIR=/Volumes/methlab_data/HBN/MRI/Site-CUNY_Derivatives_UZH/

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_CU.txt --hemi lh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/lh_area_a2009s_stats_cuny.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_CU.txt --hemi rh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/rh_area_a2009s_stats_cuny.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_CU.txt --hemi rh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/rh_thickness_a2009s_stats_cuny.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_CU.txt --hemi lh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/lh_thickness_a2009s_stats_cuny.txt


# RU:
export SUBJECTS_DIR=/Volumes/methlab_data/HBN/MRI/Site-RU_Derivatives_UZH/

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_RU.txt --hemi lh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/lh_area_a2009s_stats_ru.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_RU.txt --hemi rh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/rh_area_a2009s_stats_ru.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_RU.txt --hemi rh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/rh_thickness_a2009s_stats_ru.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ASD_RU.txt --hemi lh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/ASD/lh_thickness_a2009s_stats_ru.txt


## HC


# CBIC:
export SUBJECTS_DIR=/Volumes/methlab_data/HBN/MRI/Site-CBIC_Derivatives_UZH/

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_CB.txt --hemi rh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/rh_area_a2009s_stats_cbic.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_CB.txt --hemi lh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/lh_area_a2009s_stats_cbic.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_CB.txt --hemi rh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/rh_thickness_a2009s_stats_cbic.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_CB.txt --hemi lh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/lh_thickness_a2009s_stats_cbic.txt


# CUNY:
export SUBJECTS_DIR=/Volumes/methlab_data/HBN/MRI/Site-CUNY_Derivatives_UZH/

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_CU.txt --hemi lh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/lh_area_a2009s_stats_cuny.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_CU.txt --hemi rh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/rh_area_a2009s_stats_cuny.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_CU.txt --hemi rh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/rh_thickness_a2009s_stats_cuny.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_CU.txt --hemi lh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/lh_thickness_a2009s_stats_cuny.txt


# RU:
export SUBJECTS_DIR=/Volumes/methlab_data/HBN/MRI/Site-RU_Derivatives_UZH/

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_RU.txt --hemi lh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/lh_area_a2009s_stats_ru.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_RU.txt --hemi rh --meas area --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/rh_area_a2009s_stats_ru.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_RU.txt --hemi rh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/rh_thickness_a2009s_stats_ru.txt

aparcstats2table --subjectsfile=/Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/HC_RU.txt --hemi lh --meas thickness --parc aparc.a2009s --skip --tablefile /Users/marcdubacher/Documents/Psychologie/MA/Data/HC/lh_thickness_a2009s_stats_ru.txt

