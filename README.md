# pdb2movie
Code for generating movies of the most relevant movement modes of proteins from their PDB files, using rigidity clustering analysis

Requirements:     
PyMOL or VMD     
Python libraries: argparse     
ffmpeg   

## Instructions:      
Recompiling FIRST and diagstd is necessary!
1) Go to FIRST-190916-SAW/src and run `make`.
2) You should be ready to go! 

(Note that if your C++/Fortran compilers are not the ones we use, you might need to edit some Makefiles.)   

## Usage:     
`python pdb2movie.py FILE [options]`    
`FILE` is the PBD file of the protein for which you want to generate videos.    

## Options:     

--output PATH          - Set the output location for the videos, defaulting to `FILE` without the `.pdb`  
--overwrite            - Delete the contents of the output folder before beginning  
--forceoverwrite       - Overwrite, and silence all are-you-sure messages  
--single               - Only use chain A from the pdb file (default: use all chains)  
--keep MOLECULE LIST   - Stops the cleaning routine from removing the listed molecules from the PDB file  
--waters               - Keeps water molecules in the PDB file (equivalent to --keep HOH)  
--modes MODE LIST      - Movement modes to be investigated  
--ecuts ECUTS LIST     - Energy cutoff values  
--dstep DSTEP          - Size of directed step  
--step STEP            - Size of random step  
--confs CONFS          - Total number of configurations to be calculated  
--freq FREQ            - Frequency of saving intermediate configurations  
--nomovie              - Do not generate any videos  
--drawingengine ENGINE - Use 'vmd' or 'pymol' to render pdb files to frames (defaults to pymol for now)  
--video FILE LIST      - File(s) with PyMOL or VMD commands to be run before generating video  
--nocombi              - Prevent creation of videos combining both pos and neg directions  
--videocodec CODEC     - Use 'mp4' or 'hevc' to enode the videos, resulting in .mp4 or .mov files (default: mp4)  
--res WID HEI          - Video resolution (width, height), range [16, 8192]  
--fps FPS              - Frames per second of the videos, range [1, 240]  

See also the web-server-based implementation at https://pdb2movie.warwick.ac.uk.
