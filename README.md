# pdb2movie
Code for generating movies of the most relevant movement modes of proteins from their PDB files, using rigidity clustering analysis


Requirements:     
PyMOL     
FreeMOL mpeg_encode component     
Python libraries: glob, argparse        

Instructions:      
Recompiling FIRST is necessary!
1) Go to FIRST-190916-SAW/src and run 'make clean'.
2) run 'make'.
3) You should be ready to go!    

Usage:     
python pdb2movie.py FILE [options]    
FILE is your desired PDB file    

Options:     

--keep MOLECULE LIST - Stops the cleaning routine from removing the listed molecules from the PDB file     
--output PATH - Sets the output location for the simulations and videos    
--waters - Keeps water molecules in the PDB file (equivalent to --keep HOH)     
--confs CONFS -        Total number of configurations to be calculated     
--freq FREQ   -        Frequency of saving intermediate configurations        
--step STEP    -       Size of random step       
--dstep DSTEP    -     Size of directed step        
--res WID HEI    -     Resolution of generated videos (width, height)      
--modes MODE LIST -    Movement modes to be investigated           
--ecuts ECUTS LIST -   Energy cutoff values        
--video FILE        -  Python file with PyMOL commands to be run before generating video            
--threed    -     Generates anaglyph stereo videos     
--combi     -          Creates videos combining positive and negative directions for each mode/cutoff energy      
--multiple      -     Keeps multiple chains from the original PDB file (default: uses only chain A)   
