'''
generate_video.py - takes a series of PDB files and generates a video

'''


from __future__ import print_function
import sys
import os
#import multiprocessing
import argparse


'''
parsing_video_args: takes all command-line arguments and parse them into a structure with argument fields
This function is only used when generate_video is called as a separate script, not as part of PDB2movie!

Inputs:
list sys_args: list of system arguments as received in the command line

Outputs:
structure args: structured object with fields corresponding to the possible parameters from command line

'''
def parsing_video_args(sys_args):

    # the argparse library takes care of all the parsing from a list of command-line arguments to a structure
    parser = argparse.ArgumentParser(description='Generates videos for the most likely movement modes given a folder where the runs are stored.', usage='%(prog)s folder [options]')

    # you just need to give the parser all arguments you want to parse for and everything just works!
    parser.add_argument('--threed',  action='store_true',
                        help='Flag for generating anaglyph stereo movies')
    parser.add_argument('--res',  nargs=2, type=int, default=[640, 480], 
                        help='Video resolution (width, height), range [16, 8192]')
    parser.add_argument('--combi',  action='store_true',
                        help='Combine both positive and negative directions into a single movie')
    parser.add_argument('--modes',  nargs="+",
                        help='Movement modes to be investigated')
    parser.add_argument('--ecuts',  nargs="+",
                        help='Energy cutoff values')
    parser.add_argument('--video',  nargs=1,
                        help='File with PyMOL or VMD commands to be run before generating video')
    parser.add_argument('--videocodec', type=str, default="mp4", 
                        help="Use 'mp4' or 'hevc' to enode the videos, resulting in .mp4 or .mov files (defaults to mp4)")
    parser.add_argument('--drawingengine', type=str, default="pymol", 
                        help="Use 'vmd' or 'pymol' to render pdb files to frames (defaults to pymol for now)")
    parser.add_argument('folder', metavar='folder', type=str, nargs=1,
                        help='folder')
    parser.add_argument('--fps',  nargs=1, type=int, default=30,
                        help='Frames per second of the videos, range [1, 240]')

    # actually do the parsing for all system args other than 0 (which is the python script name) and return the structure generated
    args = parser.parse_args(sys_args[1:])
    return args



'''
call_pymol: simple wrapper for calling a Linux command

'''
'''
def call_pymol(command):
    name = multiprocessing.current_process().name
    print('--- starting:', name)
    os.system(command)
    print('--- exiting:', name)
'''



'''
gen_video: takes a structured folder full of PDB files for conformers and generate videos out of them

Inputs:
string exec_folder: folder where the python scripts are located (full path)
struct args: structure containing all arguments already parsed

'''
def gen_video(exec_folder, args):
    
    print ("---------------------------------------------------------------")
    print ("gen_video:")
    print ("----------------------------------------------------------------")

    # print(folder)

    # check for a list of modes and cutoff energies in the arguments, fills it in with defaults if no specification
    if args.modes:
        modelist = [int(x) for x in args.modes]
    else:
        modelist = range(7, 9)

    if args.ecuts:
        cutlist = [float(x) for x in args.ecuts]
    else:
        cutlist = [2.0] #[1.0, 2.0]

    modelist = [format(i, '02d') for i in modelist]
    signals = ['pos', 'neg']

    
    # ensure the width and height inputs are within allowed ranges
    if (args.res[0] < 16 or args.res[0] > 8192) or (args.res[1] < 16 or args.res[1] > 8192):
        print("width or heigth out of range [16, 8192]")
        return

    # ensure the fps input is within allowed ranges
    if (args.fps < 1 or args.fps > 240):
        print("fps out of range [1, 240]")
        return


    # ensure the drawingengine input is an allowed option
    if args.drawingengine != "pymol" and args.drawingengine != "vmd":
        print("drawingengine invalid")
        return


    # set commandfile to inputted path, including none by default
    commandfile = ""

    if args.video:
        commandfile = args.video

    # ensure the videocodec input is an allowed option
    if args.videocodec != "mp4" and args.videocodec != "hevc":
        print("videocodec invalid")
        return

    # set the fileextension according to the videocodec chosen
    fileextension = ".mp4"
    if args.videocodec == "hevc":
        fileextension = ".mov"

    # set present working directory
    folder = os.getcwd()

    print ("---------------------------------------------------------------")
    print ("gen_video: converting from pdb files to videos")
    print ("----------------------------------------------------------------")

    #for all combinations of cutoffs, modes, and signs
    for cut in cutlist:
        for mode in modelist:
            for sign in signals:

                #this is where the relevant pdbs are located
                currfolder = folder + "/Runs/"+str(cut)+"/Mode"+mode+"-"+sign

                #this is where the video will be put and what it will be called
                videoname = folder + "/Run-" + str(cut) + "-mode" + mode + "-" + sign + fileextension

                #create a link to tmp_RCD.pdb so it can be accessed as tmp_froda_00000000.pdb
                #os.system('ln -s '+currfolder+'/tmp_RCD.pdb '+currfolder+'/tmp_froda_00000000.pdb')
                os.system('cp '+currfolder+'/tmp_RCD.pdb '+currfolder+'/tmp_froda_00000000.pdb')

                #make the videos from the pdbs, using the make_video_pymol.sh or make_video_vmd.sh bash script
                print('gen_video: '+exec_folder + '/make_video_' + 
					args.drawingengine + '.sh ' + str(args.res[0]) + ' ' + str(args.res[1]) + ' ' + 
					str(args.fps) + ' ' + currfolder+' '+videoname + ' ' + args.videocodec + ' ' + 
					''.join([str(i) for i in commandfile]))
                os.system(exec_folder + '/make_video_' + 
					args.drawingengine + '.sh ' + str(args.res[0]) + ' ' + str(args.res[1]) + ' ' + 
					str(args.fps) + ' ' + currfolder+' '+videoname + ' ' + args.videocodec + ' ' + 
					''.join([str(i) for i in commandfile]))
    
            # combine the pos and neg videos into if desired
            if args.combi:

                #create a videolist file (specifying the videos to combine) to give ffmpeg
                filename = folder+"/Run-"+str(cut)+"-mode"+mode+"-"
                videolist = folder + "/Run-" + str(cut) + "-mode" + mode + "-list" 
                outF = open(videolist, "w")
                print("file " + filename+'pos'+fileextension, end="\n", file=outF)
                print("file " + filename+'neg'+fileextension, end="", file=outF)
                outF.close()

                #run ffmpeg to create the combi
                os.system('ffmpeg -hide_banner -loglevel warning -safe 0 -f concat -i ' + videolist + ' -c copy -movflags faststart -y '+filename+'combi' + fileextension)

                #fix the permissions on the combi video
                os.system('rm ' + videolist)
                os.system('chmod 744 '+filename+'combi'+fileextension)

    return



# cmd.set(full_screen='on')

# cut below here and paste into script ###
# cmd.set_view('0.979684353,   -0.145519406,   -0.137993217,\
#      0.135169536,   -0.029161714,    0.990392923,\
#     -0.148145691,   -0.988925874,   -0.008899692,\
#      0.000000000,    0.000000000, -348.124450684,\
#    -14.108730316,  -18.215091705,   66.387222290,\
#    274.463958740,  421.784942627,  -20.000000000' )



# calling this script by itself works but is inadvisable

if __name__ == "__main__":
    args = parsing_video_args(sys.argv)

    # set exec_folder to the full path of the location of this script
    exec_folder=os.path.dirname(os.path.abspath(sys.argv[0]))

    # go into the output folder
    os.chdir(args.folder[0])

    # pymol_test()
    # prepare_script(sys.argv,"t1t.mpg")
    gen_video(exec_folder, args)
