#!/bin/bash

# settings from input

pdb=${1:-1bpi}
ecuts=${2:-"2.0"}
modes=${3:-"07"}

echo "Making movie plan for "$pdb "with Ecuts=" $ecuts "and modes=" $modes

if [ -d $pdb/Runs_single-chain/ ]
then
    cd $pdb/Runs_single-chain/
elif [ -d $pdb/Runs/ ]
then
    cd $pdb/Runs/
elif [ -d $pdb/out_tidy/Runs/ ] # Flex-SARS-CoV-2-mutations
then
    cd $pdb/out_tidy/Runs/
else
    echo "MMP: did not find Runs directory --- aborting!"
fi

echo "MMP: working in "`pwd`

for ecut in $ecuts
do
    echo $ecut

    for mode in $modes
    do
	echo $mode

	pwd
	moviedir=Ecut`basename $ecut`Mode`basename $mode`
	mkdir -p $moviedir

	cd $moviedir
	#cd $ecut/Mode$mode-neg/
	for froda in ../$ecut/Mode$mode-neg/tmp_froda*.pdb
	do 
	    cp -s $froda N-`basename $froda`
	done

	for froda in ../$ecut/Mode$mode-pos/tmp_froda*.pdb
	do 
	    cp -s $froda P-`basename $froda`
	done
	cd ../

#	cd $ecut/Mode$mode-neg/
#	for froda in tmp_froda*.pdb
#	do 
#	    cp -s $froda ../../$moviedir/N-$froda
#	done
#	cd ../Mode$mode-pos/
#	for froda in tmp_froda*.pdb
#	do 
#	    cp -s $froda ../../$moviedir/P-$froda
#	done
#	cd ../../

	cd $moviedir

	ls N-*froda*.pdb | sort -rn >movie.plan-$ecut-$mode.plan
	echo "" >>movie.plan-$ecut-$mode.plan
	ls P-*froda*.pdb | sort -n >>movie.plan-$ecut-$mode.plan
 
	cd ..

    done

done

cd ../../


