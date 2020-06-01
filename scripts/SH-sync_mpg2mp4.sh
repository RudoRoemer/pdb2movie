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
    do 
	echo $mpg "--> .mp4"
	if [ ! -f $target/$dir/$dir/$dir-`basename $mpg .mpg`.mp4 ]
	then
#	    ffmpeg -hide_banner -loglevel panic drawtext="text="`basename $mpg`": fontcolor=white: fontsize=24: box=1: boxcolor=black@0.5: boxborderw=5: x=(w-text_w)/4: y=(h-text_h)/4" -i $mpg -vcodec libx264 -crf 28 $target/$dir/$dir/$dir-`basename $mpg .mpg`.mp4
	    ffmpeg -hide_banner -loglevel panic -i $mpg -vcodec libx264 -crf 28 $target/$dir/$dir/$dir-`basename $mpg .mpg`.mp4
	else
	    echo $target/$dir/$dir/$dir-`basename $mpg .mpg`.mp4 "already exists --- skipping!"
	fi
    done
    cd ../../
done

cd $current

