#!/bin/bash

for ecut in 1.0 2.0 3.0
do

for mode in 07 08 09 10 11 12
do

echo "CovID: mode=" $mode ", ecut=" $ecut ", direction=" $direct

if [ -f Run-$ecut-mode$mode-pos.mpg ] 
then
    echo "WRNG:" Run-$ecut-mode$mode-pos.mpg " does NOT exist --- skipped COMBI creation!"
    exit
fi

if [ -f Run-$ecut-mode$mode-neg.mpg ] 
then
    echo "WRNG:" Run-$ecut-mode$mode-neg.mpg " does NOT exist --- skipped COMBI creation!"
    exit
fi

if [[ -f Run-$ecut-mode$mode-combi.mpg ]]
then 
    if [[ Run-$ecut-mode$mode-combi.mpg -nt Run-$ecut-mode$mode-pos.mpg ]] && [[ Run-$ecut-mode$mode-combi.mpg -nt Run-$ecut-mode$mode-new.mpg ]] 
	then #make new combi
	    echo "--- making NEW" Run-$ecut-mode$mode-combi.mpg 
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

