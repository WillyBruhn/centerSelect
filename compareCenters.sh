#!/bin/bash

# first create all pts-files
./centerSelect.sh labParameters.txt


# second compare the pts-files
/home/sysgen/Documents/LWB/ComparingProteins/./CompareIsosurfaces.R /home/sysgen/Documents/LWB/ComparingProteins/Parameter
