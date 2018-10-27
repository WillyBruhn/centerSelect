#!/bin/bash

output="/home/willy/RedoxChallenges/centerSelect/Redox_old/Output/"

cd $output

find . -maxdepth 1 -type d  | while read dir
do
count=`ls -1 "$dir"/active_center.png 2>/dev/null | wc -l`
    if [ $count = 0 ]
    then
        echo $dir
    fi
done
