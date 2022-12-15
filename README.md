# SubtypingASD

This repository contains all the scripts and data to reproduce the Subtyping Approach. This analysis was part of the master's thesis in Psychology at the University of Zurich

## Procedure
The exactly reproduce the analysis, run the scripts in the following order:

1. subject_selection.R
2. MRI_ID_collector.sh
3. matching_IDs.R
4. aparcstats2table.sh
5. combining_statsfiles.R
6. calc_asymmetry_index.R
7. hierarchical_clustering.R
8. combining_datafiles.R
9. neuroanatomical_analysis.R
10. behavioral_analysis.R
11. prediction_analysis.R

## Folders

### Bash
Scripts to collect subject IDs from Server and to interact with FreeSurfer

### Data
.txt files from FreeSurfer processing and combination of neuroanatomical and behavioral data

### Output
Various output files created throughout the analysis

### R
Scripts for data processing and to run the analysis

## Note
GeneralParameter.R file contains definitions that are used throughout the analysis. Make sure to have the latest version of FreeSurfer (7.2.0) and R (4.2.1) installed. 
