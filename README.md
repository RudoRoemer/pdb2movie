# pdb2movie
Code for generating movies of the most relevant movement modes of proteins from their PDB files, using rigidity clustering analysis


Requirements:     
PyMOL     
FreeMOL mpeg_encode component     
 
Instructions:      
Recompiling FIRST is probably necessary to make sure it references residue tables correctly in your filesystem. Go to FIRST-190916-SAW/src and run 'make'.    

Usage:     
python pdb2movie.py FILE [options]    
FILE is your desired PDB file    

Options:     

--keep MOLECULE LIST - stops the cleaning routine from removing the listed molecules from the PDB file     
--output PATH - sets the output location for the simulations and videos    
--waters - keeps water molecules in the PDB file (equivalent to --keep HOH)     
