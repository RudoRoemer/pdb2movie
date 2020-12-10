#!/bin/bash

modes=${1:-"7"}
ecuts=${2:-"1.0"}

for ecut in $ecuts
do

for mode in $modes
do

echo "CovID: mode=" $mode ", ecut=" $ecut ", direction=" $direct

if [ -f Run-$ecut-mode$mode-combi.mpg ]
then 
    echo "WRNG:" Run-$ecut-mode$mode-combi.mpg " exists --- skipped!"
else
    if [ ! -f Run-$ecut-mode$mode-pos.mpg ]
    then 
	echo "WRNG:" Run-$ecut-mode$mode-pos.mpg " does NOT exist --- skipped COMBI creation!"
    else
	if [ ! -f Run-$ecut-mode$mode-neg.mpg ]
	then 
	    echo "WRNG:" Run-$ecut-mode$mode-neg.mpg " does NOT exist --- skipped COMBI creation!"
	else
	    echo "--- making " Run-$ecut-mode$mode-combi.mpg 
	    cat Run-$ecut-mode$mode-pos.mpg Run-$ecut-mode$mode-neg.mpg > Run-$ecut-mode$mode-combi.mpg 
	fi
    fi
fi

done
done

