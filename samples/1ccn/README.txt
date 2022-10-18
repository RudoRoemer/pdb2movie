###############################################################################
# README
###############################################################################

Sample of a PDB file and output produced by pdb2movie

###############################################################################
# Instructions
###############################################################################

> cd pdb2movie
> cp samples/1ccn.pdb .

Then run as follows

> python3 pdb2movie.py 1ccn.pdb --combi --confs 10 --freq 50 --modes 7 8 --ecuts 2.0 --res 1920 1080 --multiple --videocodec mp4 --drawingengine pymol

This produces the .mp4 movies given in samples/1ccn. For .mov movies, do

> python3 pdb2movie.py 1ccn.pdb --combi --confs 10 --freq 50 --modes 7 8 --ecuts 2.0 --res 1920 1080 --multiple --videocodec hevc --drawingengine pymol

For the VMD produced movies in samples/1ccn/vmd_movies/, use "vmd" instead of "pymol" in the last command line argument or simply remove the "--drawingengine" argument.
