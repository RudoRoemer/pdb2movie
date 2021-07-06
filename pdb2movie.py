'''
pdb2movie.py - this is where the main program lives! 
it binds together all other functions defined in the other python scripts.
'''

import sys
import os
import cleanpdb
import runfirst
import runelnemo
import runfroda
import generate_video_serial
import argparse

'''
parsing_args: takes all command-line arguments and parse them into a structure with argument fields

Inputs:
list sys_args: list of system arguments as received in the command line

Outputs:
structure args: structured object with fields corresponding to the possible parameters from command line
'''

def parsing_args(sys_args):

    # the argparse library takes care of all the parsing from a list of
    # command-line arguments to a structure
    parser = argparse.ArgumentParser(description='Runs simulations and generates videos for the most likely movement modes given a PDB file.',usage='%(prog)s pdbfile [options]')

    parser.add_argument('--keep',  nargs="+",
                        help='List of molecules to be kept')
    parser.add_argument('--output',  nargs=1,
                        help='Output directory')
    parser.add_argument('--res',  nargs=2, type=int, default=[640, 480], 
                        help='Video resolution (width, height), range [16, 8192]')
    parser.add_argument('--waters',  action='store_true',
                        help='Flag for keeping water molecules')
    parser.add_argument('--multiple',  action='store_true',
                        help='Keep multiple chains (default: uses only chain A)')
    parser.add_argument('--nomovie',  action='store_true',
                        help='Stop calculation after FRODA and before any attempt to generate images and movies')
    parser.add_argument('--combi',  action='store_true',
                        help='Combine both positive and negative directions into a single movie')
    parser.add_argument('--threed',  action='store_true',
                        help='Flag for generating anaglyph stereo movies')
    parser.add_argument('--confs',  nargs=1,
                        help='Total number of configurations to be calculated')
    parser.add_argument('--freq',  nargs=1,
                        help='Frequency of saving intermediate configurations')
    parser.add_argument('--step',  nargs=1,
                        help='Size of random step')
    parser.add_argument('--dstep',  nargs=1,
                        help='Size of directed step')
    parser.add_argument('--modes',  nargs="+",
                        help='Movement modes to be investigated')
    parser.add_argument('--ecuts',  nargs="+",
                        help='Energy cutoff values')
    parser.add_argument('--video',  nargs=1,
                        help='File with PyMOL or VMD commands to be run before generating video')
    parser.add_argument('pdbfile', metavar='PDB', type=str, nargs=1,
                        help='Initial PDB file')
    parser.add_argument('--videocodec', type=str, default="mp4", 
                        help="Use 'mp4' or 'hevc' to enode the videos, resulting in .mp4 or .mov files (defaults to mp4)")
    parser.add_argument('--drawingengine', type=str, default="pymol", 
                        help="Use 'vmd' or 'pymol' to render pdb files to frames (defaults to pymol for now)")
    parser.add_argument('--fps',  nargs=1, type=int, default=30,
                        help='Frames per second of the videos, range [1, 240]')

    # actually do the parsing for all system args other than 0
    # (which is the python script name) and return the structure generated
    args = parser.parse_args(sys_args[1:])

    #ensure pdbfile and output are full paths, not relative ones (so they are correct from here on, regardless of pwd)
    args.pdbfile[0] = os.path.abspath(args.pdbfile[0])
    args.output[0] = os.path.abspath(args.output[0])
    
    return args

# no specific function is defined here, since this is our main program!

if __name__ == "__main__":#

    # first things first: we need to parse command-line arguments with the function we have defined
    args=parsing_args(sys.argv)

    # print "test"
    # if (args.threed):
    #     userinput=raw_input("WARNING: PyMOL windows will open during generation and they won't close by themselves. You have been warned. Are you sure you want to continue? [y/n]   ")
    #     if (userinput!="y"):
    #         quit()

    # set exec_folder to the full path of this script
    exec_folder=os.path.dirname(os.path.abspath(sys.argv[0]))

    #print(exec_folder)

    # if an output folder was defined, we need to make that directory if it doesn't exist,
    # and warn the user that the folder will be emptied if it exists
    # and then make that the present working directory
    if (args.output):

        try:
            os.mkdir(args.output[0])
            os.chdir(args.output[0])
        except Exception:
            userinput=raw_input("WARNING: everything in output folder will be deleted! Are you sure you want to continue? [y/N]   ")
            if (userinput=="y"):
                os.system("rm -r "+args.output[0]+"/*")
                os.chdir(args.output[0])
            else:
                quit()
            pass

    # if we haven't specified an output folder, we make a subfolder at current directory with the protein name (or empty an existing subfolder)
    else:

        try:
            os.mkdir(args.pdbfile[0][:-4])
            os.chdir(args.pdbfile[0][:-4])
        except Exception:
            os.system("rm -r "+args.pdbfile[0][:-4]+"/*")
            os.chdir(args.pdbfile[0][:-4])
            pass

    # now we will just call the functions defined in other files in sequence
    # and that's all we need to do!

    # first we clean up the PDB file...
    print ("---------------------------------------------------------------")
    print ("pdb2movie: cleaning the PDB file")
    print ("----------------------------------------------------------------")
    clean_file=cleanpdb.cleanPDB(args,exec_folder)

    # then we run FIRST on it
    print ("---------------------------------------------------------------")
    print ("pdb2movie: starting rigidity analysis with FIRST")
    print ("----------------------------------------------------------------")
    hydro_file=runfirst.firstsim(exec_folder,args,clean_file)

    # we get the folder path because we forgot to save it when creating/emptying it beforehand
    folder=hydro_file.rsplit("/",1)[0]
    
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
