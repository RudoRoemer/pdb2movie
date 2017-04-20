#run this script by typing pymol view.py at the command line in the directory that contains the pdb files you want to view as a movie.
#edit name appropriately - this is for the standard input_froda_* files

#this is to view all of the pdb files in your current working directory in number order, using pymol.
from glob import glob
import sys
import os

def gen_video(args,folder):
    cutlist=[1.0, 2.0]
    modelist=range(7,12)
    modelist=[format(i, '02d') for i in modelist]
    signals=['pos','neg']
    for cut in cutlist:
        for mode in modelist:
            for sign in signals:

                # Desired pymol commands here to produce and save figures

                currfolder=folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/"
                os.system("pymol -cq pymolvideo.py -- "+currfolder)





#cmd.set(full_screen='on')

### cut below here and paste into script ###
# cmd.set_view('0.979684353,   -0.145519406,   -0.137993217,\
#      0.135169536,   -0.029161714,    0.990392923,\
#     -0.148145691,   -0.988925874,   -0.008899692,\
#      0.000000000,    0.000000000, -348.124450684,\
#    -14.108730316,  -18.215091705,   66.387222290,\
#    274.463958740,  421.784942627,  -20.000000000' )




if __name__ == "__main__":#
    cutlist=[1.0, 2.0]
    modelist=range(7,12)
    # print sys.argv
    folder=sys.argv[1]
    # pymol_test()
    gen_video(folder,cutlist,modelist)
