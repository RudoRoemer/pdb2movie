#!/bin/bash 
pdb=$1 
for chain in $(grep "^ATOM" $pdb | cut -b 22 | sort -u) 
do 
    sed -n "/^.\{21\}$chain/p" $pdb > ${pdb%.pdb}_$chain.pdb 
done 
