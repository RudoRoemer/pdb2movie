#run this script by typing pymol view.py at the command line in the directory that contains the pdb files you want to view as a movie.
#edit name appropriately - this is for the standard input_froda_* files
import pymol
#this is to view all of the pdb files in your current working directory in number order, using pymol.
from glob import glob
import sys
import __main__


def gen_video(folder,cutlist,modelist):
    __main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
    pymol.finish_launching()
    modelist=[format(i, '02d') for i in modelist]
    signals=['pos','neg']
    for cut in cutlist:
        for mode in modelist:
            for sign in signals:
                pymol.cmd.reinitialize()
                # Desired pymol commands here to produce and save figures
                sleep(0.5) # (in seconds)
                currfolder=folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/"
                lst = glob(currfolder+"input*.pdb")
                lst.sort()
                for fil in lst: pymol.cmd.load(fil,"mov")
                pymol.cmd.hide()
                pymol.cmd.show('cartoon')
                pymol.cmd.bg_color('white')
                pymol.cmd.movie.add_state_sweep(2,1,'start=1')
                pymol.cmd.movie.produce(currfolder+'cartoon.mpg',mode='ray',quality=100,quiet=1)
                pymol.cmd.save('cartoon.pse','','pse',quiet=0)




#cmd.set(full_screen='on')

### cut below here and paste into script ###
# cmd.set_view('0.979684353,   -0.145519406,   -0.137993217,\
#      0.135169536,   -0.029161714,    0.990392923,\
#     -0.148145691,   -0.988925874,   -0.008899692,\
#      0.000000000,    0.000000000, -348.124450684,\
#    -14.108730316,  -18.215091705,   66.387222290,\
#    274.463958740,  421.784942627,  -20.000000000' )

### cut above here and paste into script ###
def pymol_test():
    # __main__.pymol_argv = ['./usr/bin/pymol','-qc'] # Pymol: quiet and no GUI
    pymol.finish_launching()


if __name__ == "__main__":#
    cutlist=[1.0, 2.0]
    modelist=range(7,12)
    # print sys.argv
    folder=sys.argv[1]
    pymol_test()
    # gen_video(folder,cutlist,modelist)
