#!/bin/bash

# parametersName="lab/Parameters"
parametersName="willy/Parameters"

path2ComparingProteins="/home/willy/RedoxChallenges/ComparingProteins/"



# first create all pts-files

./centerSelect.sh ./parameters/"$parametersName"20.txt
$path2ComparingProteins./CompareIsosurfaces.R ./parameters/willy/ComparingProteinsParameters20.txt



./centerSelect.sh ./parameters/"$parametersName"30.txt
$path2ComparingProteins./CompareIsosurfaces.R ./parameters/willy/ComparingProteinsParameters30.txt



./centerSelect.sh ./parameters/"$parametersName"40.txt
$path2ComparingProteins./CompareIsosurfaces.R ./parameters/willy/ComparingProteinsParameters40.txt



##!/bin/bash
#
## first create all pts-files
#./centerSelect.sh labParameters.txt
#
#
## second compare the pts-files
#/home/sysgen/Documents/LWB/ComparingProteins/./CompareIsosurfaces.R /home/sysgen/Documents/LWB/ComparingProteins/Parameter
