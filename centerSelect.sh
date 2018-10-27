#!/bin/bash

path2centerSelect="/home/willy/RedoxChallenges/centerSelect/"
path2files="/home/willy/RedoxChallenges/centerSelect/Redox_old/Output/"
boxSize=30

# point[3] > center[3] - depth; positive depth means more points
depth=10
eps=0.3


missingList=$(./WhichMissing.sh)

cd $path2files

#for d in */ ; do
for d in $missingList; do
    #echo "$d"
    
    filename="${d%/}"

	f="$path2files/$filename/$filename"

	rm "$f""HeadOff.pqr"

	sed '1,1d' "$f.pqr"  >> $f."tmp.pqr"
	sed '$d' $f."tmp.pqr"  >> "$f""HeadOff.pqr"

	rm $f."tmp.pqr"
 
	
	$path2centerSelect/./centerSelecter.R $path2centerSelect $path2files $filename $boxSize $depth $eps

done


echo "still missing active centers of ..."
cd $path2centerSelect
./WhichMissing.sh
