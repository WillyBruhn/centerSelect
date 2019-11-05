#!/bin/bash

output="/home/willy/RedoxChallenges/centerSelect/Redox_old/Output/"
output=$1
#echo "-----------------------------------------------------------"
#echo "WhichMissing.sh called on $output"
#echo "-----------------------------------------------------------"


cd $output

find . -maxdepth 1 -type d  | while read dir
do
count=`ls -1 "$dir"/*_active_center.pts 2>/dev/null | wc -l`
    if [ $count = 0 ]
    then
        echo $dir
    fi
done
