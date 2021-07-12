'''
pdb2movie.py - this is where the main program lives! 
it binds together all other functions defined in the other python scripts.
'''

import sys
import os

import helpers
import cleanpdb
import runfirst
import runelnemo
import runfroda
import generate_video_serial


# no specific function is defined here, since this is our main program!

if __name__ == "__main__":#

    # parse commmand-line arguments
    args=helpers.parsing_args(sys.argv)

    # print "test"
    # if (args.threed):
    #     userinput=raw_input("WARNING: PyMOL windows will open during generation and they won't close by themselves. You have been warned. Are you sure you want to continue? [y/n]   ")
    #     if (userinput!="y"):
    #         quit()

    # set exec_folder to the full path of this script
    exec_folder=os.path.dirname(os.path.abspath(sys.argv[0]))

    # set the output folder if not specified
    if (not args.output):
        args.output = [args.pdbfile[0][:-4]]

    # change directory to the output folder
    helpers.go_to_output_folder(args)

    # now we will just call the functions defined in other files in sequence
    # and that's all we need to do!

    # first we clean up the PDB file...
    print ("---------------------------------------------------------------")
    print ("pdb2movie: cleaning the PDB file")
    print ("----------------------------------------------------------------")
    clean_file=cleanpdb.cleanPDB(args)

    # then we run FIRST on it
    print ("---------------------------------------------------------------")
    print ("pdb2movie: starting rigidity analysis with FIRST")
    print ("----------------------------------------------------------------")
    hydro_file=runfirst.firstsim(exec_folder,args,clean_file)

    # we get the folder path because we forgot to save it when creating/emptying it beforehand
    folder=os.path.abspath(hydro_file.rsplit("/",1)[0])

    # these don't even have outputs that need to be saved:
    # first we run ElNemo...
    print ("---------------------------------------------------------------")
    print ("pdb2movie: starting the eleastic network modelling")
    print ("----------------------------------------------------------------")
    runelnemo.elnemosim(exec_folder,args,hydro_file)

    # then FRODA...
    print ("---------------------------------------------------------------")
    print ("pdb2movie: starting the FRODA analysis")
    print ("----------------------------------------------------------------")
    runfroda.frodasim(exec_folder,args,hydro_file)

    if args.nomovie:
        # --nomovie setting
        print ("---------------------------------------------------------------")
        print ("pdb2movie: NOT generating any videos")
        print ("----------------------------------------------------------------")
    else:
        # We generate the videos
        print ("---------------------------------------------------------------")
        print ("pdb2movie: generating the videos")
        print ("----------------------------------------------------------------")
        generate_video_serial.gen_video(exec_folder,args,folder)
        
    # finally, we are done
    print ("---------------------------------------------------------------")
    print ("pdb2movie: all should be done, enjoy!")
    print ("----------------------------------------------------------------")
