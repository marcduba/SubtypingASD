#!/bin/bash

## Creating list of subject IDs for every site
# CBIC:
cd /Volumes/methlab_data/HBN/MRI/Site-CBIC_Derivatives_UZH/
ls -d * > /Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ID_CB.txt

# CUNY:
cd /Volumes/methlab_data/HBN/MRI/Site-CUNY_Derivatives_UZH/
ls -d * > /Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ID_CU.txt

# RU:
cd /Volumes/methlab_data/HBN/MRI/Site-RU_Derivatives_UZH/
ls -d * > /Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ID_RU.txt

# SI:
cd /Volumes/methlab_data/HBN/MRI/Site-SI_Derivatives_UZH/
ls -d * > /Users/marcdubacher/Documents/Psychologie/MA/Data/IDlists/ID_SI.txt


