#!/bin/bash

source=${1:-"./"}
target=${2:-"./"}
pdb=${3:-"aaaa"}
#overwrite=${4:-"n"}

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
    cd Runs

    echo `pwd`

    allecuts=`ls -1d *`
    echo $allecuts

    for ecut in $allecuts
    do 
        newestfile=`ls -tr 3.0 | tail -1`
	echo $newestfile
        targetzip=$target/$dir/$dir/$dir-$ecut.zip 
	echo $targetzip

	if [[ ! -f $targetzip ]] || [[ $newestfile -nt $targetzip ]]
	then
	    cd $ecut
	   echo $ecut "--> .zip"
           zip -uqR $targetzip "tmp_froda*.pdb"
	   cd ..
	else
	    echo "NEW" $targetzip "already exists --- skipping!"
	fi
    done
    cd ../../../
done

cd $current

