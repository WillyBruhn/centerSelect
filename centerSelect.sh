#!/bin/bash

parametersFile="./parameters.txt"
if [ $# -eq 0 ]
  then
    echo "No arguments supplied ..."
else
	parametersFile=$1
fi

echo "reading in parameters from $parametersFile"

echo "found the following parameters..."
echo "paramameter          value"
exit
NLS="-----------------------------------------------------------"
echo $NLS
#-----------------------------------------------------------------------------------
IFS=";"
while read parameter value
do
    echo "$parameter              $value"

	if [ "$parameter" == "path2centerSelect" ]; then
		path2centerSelect=$value

	elif [ "$parameter" == "path2files" ]; then
		path2files=$value

	elif [ "$parameter" == "boxSize" ]; then
		boxSize=$value


	elif [ "$parameter" == "depth" ]; then
		depth=$value


	elif [ "$parameter" == "eps" ]; then
		eps=$value


	elif [ "$parameter" == "eps" ]; then
		eps=$value


	elif [ "$parameter" == "onlyDoMissing" ]; then
		onlyDoMissing=$value

	else 
		echo "corrupted format in $parametersfile/parameters.txt"
		echo "aborting ..."
		exit
	fi

done < $parametersFile
echo $NLS
#-----------------------------------------------------------------------------------

#path2centerSelect="/home/willy/RedoxChallenges/centerSelect/"
#path2files="/home/willy/RedoxChallenges/centerSelect/Redox_old/Output/"
#boxSize=30
#
## point[3] > center[3] - depth; positive depth means more points
#depth=10
#eps=0.3

function selectCenter {
	d=$1

	echo $d
	echo $filename

    filename="${d%/}"

	f="$path2files/$filename/$filename"

	rm "$f""HeadOff.pqr"

	sed '1,1d' "$f.pqr"  >> $f."tmp.pqr"
	sed '$d' $f."tmp.pqr"  >> "$f""HeadOff.pqr"

	rm $f."tmp.pqr"
 
	
	$path2centerSelect/./centerSelecter.R $path2centerSelect $path2files $filename $boxSize $depth $eps
}


#missingList=$($path2centerSelect./WhichMissing.sh $path2files)


#while read -r line; do
#    echo "... $line ..."
#done <<< "$missingList"

#exit

cd $path2files

#for d in */ ; do

echo "onlyDoMissing == $onlyDoMissing ..."

if [ $onlyDoMissing == true ]; then
	
	echo "... performing only for missing folders"

	missingList=$($path2centerSelect./WhichMissing.sh $path2files)

	while read -r d; do
		selectCenter $d
	done <<< "$missingList"
else 

	echo "... performing for all folders"

	for d in */ ; do
		selectCenter $d
	done
fi


echo $NLS
#-----------------------------------------------------------------------------------
echo "still missing active centers of ..."
$path2centerSelect./WhichMissing.sh $path2files




