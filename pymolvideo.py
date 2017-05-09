#run this script by typing pymol view.py at the command line in the directory that contains the pdb files you want to view as a movie.
#edit name appropriately - this is for the standard input_froda_* files
from pymol import cmd
#this is to view all of the pdb files in your current working directory in number order, using pymol.
from glob import glob
# print sys.argv
folder=sys.argv[1]
prot=folder.rsplit("/",1)[1]
lst = glob(folder+"/*.pdb")
lst.sort()
for fil in lst: cmd.load(fil,"mov")
filename=tst.mpg

#cmd.select( 'flexchain', 'resi 22-30' )
#cmd.hide( 'cartoon', 'flexchain' )
#cmd.select( 'a_domain', 'resi 31-140' )
#cmd.color( 'blue', 'a_domain' )
#cmd.select( 'b_domain', 'resi 141-240' )
#cmd.color( 'green', 'b_domain' )
#cmd.select( 'bprime_domain', 'resi 241-359' )
#cmd.color( 'yellow', 'bprime_domain' )
#cmd.select( 'x_linker', 'resi 360-374' )
#cmd.color( 'black', 'x_linker' )
#cmd.select( 'aprime_domain', 'resi 375-484' )
#cmd.color( 'orange', 'aprime_domain' )
#cmd.select( 'c_domain', 'resi 485-504' )
#cmd.color( 'red', 'c_domain' )
#cmd.select( 'residue61', 'resi 61 and name ca')
#cmd.color( 'yellow', 'residue61')
#cmd.show( 'sphere', 'residue61' )
#cmd.select( 'residue406', 'resi 406 and name ca')
#cmd.color( 'yellow', 'residue406')
#cmd.show( 'sphere', 'residue406' )


#cmd.set(full_screen='on')

### cut below here and paste into script ###
# cmd.set_view('0.979684353,   -0.145519406,   -0.137993217,\
#      0.135169536,   -0.029161714,    0.990392923,\
#     -0.148145691,   -0.988925874,   -0.008899692,\
#      0.000000000,    0.000000000, -348.124450684,\
#    -14.108730316,  -18.215091705,   66.387222290,\
#    274.463958740,  421.784942627,  -20.000000000' )

### cut above here and paste into script ###
cmd.hide()
cmd.show('cartoon')
cmd.bg_color('white')
cmd.movie.add_state_sweep(2,1,'start=1')
cmd.movie.produce(filename,mode='ray',quality=100,quiet=1)
cmd.save(folder+'cartoon.pse')

print 'CARTOON rendering finished'

cmd.quit()
