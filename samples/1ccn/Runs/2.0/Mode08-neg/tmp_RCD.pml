# Rigid cluster decomposition coloring script for PyMol
#
# Created by Dan Farrell, Brandon Hespenheide.
# Department of Physics and Astronomy
# Biophysics Theory Group
# Arizona State University
############################################################

from pymol import cmd
from pymol.cgo import *

bg_color white
# script to make a movie of the FRODA generated structures.

from glob import glob
filelist = glob ("tmp_froda_*.pdb") 
filelist.sort()
cmd.load("tmp_RCD.pdb", "tmp_froda")
for file in filelist: cmd.load( file, "tmp_froda") 
show lines, tmp_froda
color black
color 0x0000b2, ( b > 0.99 and b < 1.01)
color 0x57eb0f, ( b > 1.99 and b < 2.01)
color 0xe75ab7, ( b > 2.99 and b < 3.01)
# Draw hbonds, hydrophoic tethers and stacked rings as distance objects
set dash_gap, 0.1
distance hbonds = id 00291 , id 00233
distance hbonds = id 00431 , id 00369
distance hbonds = id 00037 , id 00455
distance hbonds = id 00393 , id 00340
distance hbonds = id 00428 , id 00369
distance hbonds = id 00609 , id 00567
distance hbonds = id 00181 , id 00114
distance hbonds = id 00211 , id 00148
distance hbonds = id 00381 , id 00325
distance hbonds = id 00116 , id 00070
distance hbonds = id 00151 , id 00079
distance hbonds = id 00198 , id 00124
distance hbonds = id 00498 , id 00004
distance hbonds = id 00448 , id 00379
distance hbonds = id 00132 , id 00068
distance hbonds = id 00010 , id 00493
distance hbonds = id 00013 , id 00526
distance hbonds = id 00164 , id 00098
distance hbonds = id 00460 , id 00034
distance hbonds = id 00226 , id 00159
distance hbonds = id 00372 , id 00311
distance hbonds = id 00412 , id 00350
color red, hbonds
hide labels, hbonds
disable hbonds
distance hydrophobics = id 00023 , id 00476
distance hydrophobics = id 00035 , id 00601
distance hydrophobics = id 00036 , id 00495
distance hydrophobics = id 00045 , id 00125
distance hydrophobics = id 00046 , id 00177
distance hydrophobics = id 00177 , id 00447
distance hydrophobics = id 00179 , id 00370
distance hydrophobics = id 00312 , id 00351
distance hydrophobics = id 00354 , id 00408
distance hydrophobics = id 00380 , id 00478
distance hydrophobics = id 00457 , id 00604
distance hydrophobics = id 00458 , id 00606
distance hydrophobics = id 00459 , id 00496
color green, hydrophobics
hide labels, hydrophobics
disable hydrophobics
# Rigid Cluster 1 has 80 atoms.
create RC1, ( b > 0.99 and b < 1.01)
show sticks, RC1
set line_width = 3, RC1
color 0x0000b2, RC1

# Rigid Cluster 2 has 65 atoms.
create RC2, ( b > 1.99 and b < 2.01)
show sticks, RC2
set line_width = 3, RC2
color 0x57eb0f, RC2

# Rigid Cluster 3 has 24 atoms.
create RC3, ( b > 2.99 and b < 3.01)
show sticks, RC3
set line_width = 3, RC3
color 0xe75ab7, RC3

# Rigid Cluster BIN2
create BIN2, ( b > 3.99 and b < 8.01)
show sticks, BIN2
set line_width = 3, BIN2
color gray, BIN2
disable BIN2

# Rigid Cluster BIN1
create BIN1, ( b > 8.99 and b < 147.01)
show sticks, BIN1
set line_width = 3, BIN1
color gray, BIN1
disable BIN1

