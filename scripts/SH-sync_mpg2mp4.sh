#!/bin/bash

source=${1:-"./"}
target=${2:-"./"}
pdb=${3:-"aaaa"}

current=`pwd`

echo $source
echo $target
echo $pdb

cd $source

for dir in $pdb
do 
    echo "--- working on "$dir
    cd $dir
    cd $dir

    pwd

    for mpg in *-combi.mpg
    do echo $mpg
	if [ ! -f $target/$dir/$dir/$dir-`basename $mpg .mpg`.mp4 ]
	then
	    ffmpeg -hide_banner -loglevel panic -i $mpg -vcodec libx264 -crf 28 $target/$dir/$dir/$dir-`basename $mpg .mpg`.mp4
	else
	    echo $target/$dir/$dir/$dir-`basename $mpg .mpg`.mp4 "already exists --- skipping!"
	fi
    done
    cd ../..
done

cd $currentdir

