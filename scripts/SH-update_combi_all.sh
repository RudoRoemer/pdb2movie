#!/bin/bash

allecuts=`ls -1 *.mpg | cut -d "-" -f 2 | sort | uniq`
allmodes=`ls -1 *.mpg | cut -d "-" -f 3 | cut -b 5-6 | sort | uniq`

echo "ecuts=" $allecuts ", allmodes=" $allmodes

for ecut in $allecuts
do

for mode in $allmodes
do

echo "CovID: mode=" $mode ", ecut=" $ecut ", direction=" $direct

if [ ! -f Run-$ecut-mode$mode-pos.mpg ] 
then
    echo "WRNG:" Run-$ecut-mode$mode-pos.mpg " does NOT exist --- skipping COMBI creation!"
    break
fi

if [ ! -f Run-$ecut-mode$mode-neg.mpg ] 
then
    echo "WRNG:" Run-$ecut-mode$mode-neg.mpg " does NOT exist --- skipping COMBI creation!"
    exit
fi

if [[ -f Run-$ecut-mode$mode-combi.mpg ]]
then 
    if [[ Run-$ecut-mode$mode-combi.mpg -ot Run-$ecut-mode$mode-pos.mpg ]] || [[ Run-$ecut-mode$mode-combi.mpg -ot Run-$ecut-mode$mode-neg.mpg ]] 
	then #make new combi
	    echo "--- UPDATING" Run-$ecut-mode$mode-combi.mpg 
	    cat Run-$ecut-mode$mode-pos.mpg Run-$ecut-mode$mode-neg.mpg > Run-$ecut-mode$mode-combi.mpg 
    else
	echo "WRNG: newer" Run-$ecut-mode$mode-combi.mpg " exists --- skipped!"
    fi
else
    echo "--- making FIRST" Run-$ecut-mode$mode-combi.mpg 
    cat Run-$ecut-mode$mode-pos.mpg Run-$ecut-mode$mode-neg.mpg > Run-$ecut-mode$mode-combi.mpg 
fi

done
done

