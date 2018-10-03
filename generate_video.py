'''
generate_video.py - takes a series of PDB files and generates a video

'''




import sys
import os
import multiprocessing
import argparse
import subprocess


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
    parser.add_argument('--res',  nargs=2,
                        help='Video resolution (width, height)')
    parser.add_argument('--combi',  action='store_true',
                        help='Combine both positive and negative directions into a single movie')
    parser.add_argument('--modes',  nargs="+",
                        help='Movement modes to be investigated')
    parser.add_argument('--ecuts',  nargs="+",
                        help='Energy cutoff values')
    parser.add_argument('--video',  nargs=1,
                        help='Python file with PyMOL commands to be run before generating video')
    parser.add_argument('folder', metavar='PDB', type=str, nargs=1,
                        help='Folder where the runs can be found')

    # actually do the parsing for all system args other than 0 (which is the python script name) and return the structure generated
    args = parser.parse_args(sys_args[1:])
    return args



'''
call_pymol: simple wrapper for calling a Linux command
'''

def call_pymol(command):
    os.system(command)



'''
gen_video: takes a folder full of PDB files for conformers and generate a video out of them

Inputs:
string exec_folder: folder where the python scripts are located (full path)
struct args: structure containing all arguments already parsed
string folder: folder where the PDB files to be turned into a video are located (full path)


'''


def gen_video(exec_folder, args, folder):
    
    print ("---------------------------------------------------------------")
    print ("gen_video:")
    print ("----------------------------------------------------------------")

    # print(folder)

    # check for a list of modes and cutoff energies in the arguments, fills it in with defaults if no specification
    if args.modes:
        modelist = [int(x) for x in args.modes]
    else:
        modelist = range(7, 12)

    if args.ecuts:
        cutlist = [float(x) for x in args.ecuts]
    else:
        cutlist = [1.0, 2.0]

    modelist = [format(i, '02d') for i in modelist]
    signals = ['pos', 'neg']


    # initialise a list of jobs for multiprocessing 
    jobs = []


    # loop over all combinations of cutoffs, modes, signs
    for cut in cutlist:
        for mode in modelist:
            for sign in signals:

                # generates a filename with the correct cutoff, mode and sign
                filename = folder + "/Run-"+str(cut)+"-mode"+mode+"-"+sign+".mpg"
                print (filename)


                # generates a separate script for this combination of cutoff, mode and sign
                prepare_script(exec_folder, args, filename, cut, mode, sign, folder)
                # Desired pymol commands here to produce and save figures

                currfolder = folder + "/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/"

                # calls pymol with the generated script, with a new process for each combination of cutoffs, modes and signs
                if args.threed:
                    command = 'pymol -q '+folder+'/pymolvideo'+str(cut)+mode+sign+'.py -- '+currfolder
                else:
                    command = 'pymol -cq '+folder+'/pymolvideo'+str(cut)+mode+sign+'.py -- '+currfolder
                p = multiprocessing.Process(target=call_pymol, args=(command,))
                jobs.append(p)
                p.start()


    # this is here so that we wait for every process to finish before proceeding
    for job in jobs:
        job.join()

    # checks whether freemol is present in the system, if it is not we will have some extra work
    if (os.system('grep \'FREEMOL\' $(which pymol)')):

        # if freemol is not there, we will loop over all combinations once more, generate the correct filename and so on
        for cut in cutlist:
            for mode in modelist:
                for sign in signals:
                    filename = folder+"/Run-"+str(cut)+"-mode"+mode+"-"+sign+".mp4"
                    tmpfolder = filename.rsplit("/", 1)[1][:-3]

                    currfolder = folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/"
                    # command = ['convert', '-quality', ' 100', folder+'/'+tmpfolder+'tmp/*.ppm', filename]


                    # finally, we will generate a video using ffmpeg instead of freemol, based on the ppm screenshots pymol has generated before
                    command = 'ffmpeg -pattern_type glob -i '+'\"'+folder+'/'+tmpfolder+'tmp/*.ppm'+'\" -c:v copy -pix_fmt yuv420p -me_method epzs -threads 4 -r 30.000030 -g 45 -bf 2 -trellis 2 -y -b:v 6000k '+filename
                    
                    # command = ['ffmpeg', '-pattern_type', 'glob', '-i', '\"'+folder+'/'+tmpfolder+'tmp/*.ppm'+'\"', '-c:v', 'mpeg2video', '-pix_fmt', 'yuv420p', '-me_method', 'epzs', '-threads', '4', '-r', '30.000030', '-g', '45', '-bf', '2', '-trellis', '2', '-y', '-b', '6000k', filename]
                    # print(command)
                    os.system("bash -c '{0}'".format(command))
                    # subprocess.call(command)




    # now we loop over cutoffs and modes, and if we want combined movies we do that purely by concatenating two videos
    for cut in cutlist:
        for mode in modelist:
            if args.combi:
                filename = folder+"/Run-"+str(cut)+"-mode"+mode+"-"
                os.system('cat '+filename+'pos.mpg '+filename+'neg.mpg > '+filename+'combi.mpg')
                os.system('chmod 744 '+filename+'combi.mpg')

            # we also need to fix permissions for the all the videos 
            for sign in signals:
                # os.system('rm '+folder+'/pymolvideo'+str(cut)+mode+sign+'.py')
                filename = folder+"/Run-"+str(cut)+"-mode"+mode+"-"+sign+".mpg"
                os.system('chmod 744 '+filename)
                tmpfolder = filename.rsplit("/", 1)[1][:-3]
                os.system('rm -r '+folder+'/'+tmpfolder+'tmp/')
    return

# cmd.set(full_screen='on')

# cut below here and paste into script ###
# cmd.set_view('0.979684353,   -0.145519406,   -0.137993217,\
#      0.135169536,   -0.029161714,    0.990392923,\
#     -0.148145691,   -0.988925874,   -0.008899692,\
#      0.000000000,    0.000000000, -348.124450684,\
#    -14.108730316,  -18.215091705,   66.387222290,\
#    274.463958740,  421.784942627,  -20.000000000' )





'''
prepare_script: generates a python script to be used in pymol for generating a specific video

Inputs:
string exec_folder: folder where the python scripts are located (full path)
struct args: structure containing all arguments already parsed
string filename: path to the resulting movie at the end of the video generating process
float cut: cutoff energy for this video
int mode: movement mode for this video
string sign: direction of movement for this video (positive or negative)
string folder: folder where the PDB files to be turned into a video are located (full path)

'''


def prepare_script(exec_folder, args, filename, cut, mode, sign, folder):

    # string="cat video_template.py <(echo filename=\'"+filename+"\') video_minimal.py >pymolvideo.py"
    # string='cat video_template.py <(echo \"stereo anaglyph\") <(echo filename=\\"'+filename+'\\") ' +args.video[0]+' video_minimal.py > pymolvideo.py'


    # the way this works is we generate a pretty big string with a command that will concatenate a bunch of stuff into a py file.
    # the start of this py file is video_template.py

    string = 'cat '+exec_folder+'/video_template.py '

    # if optional resolution arguments were passed, we concatenate that into our final script
    if args.res:
        string = string+'<(echo \"cmd.viewport('+str(args.res[0])+','+str(args.res[1])+')\") '

    # if the videos are to be generated in 3d, we need to set stereo options in pymol as well
    if args.threed:
        string = string+'<(echo \"cmd.set(\\"stereo_mode\\",10)\") <(echo \"cmd.stereo(\\"on\\")\") '

    # this concatenates a line specifying the filename for the generated movie
    string = string+'<(echo filename=\\"'+filename+'\\") '

    # if we are given an external file with commands to run before generating a video, we concatenate that as well
    if args.video:
        string = string+args.video[0]+' '

    # the end of the file is on video_minimal.py, where the actual video generation lives, so we concatenate it as well into the final python file
    string = string+exec_folder+'/video_minimal.py > '+folder+'/pymolvideo'+str(cut)+mode+sign+'.py'

    # then, we run the big string command!
    os.system("bash -c '{0}'".format(string))




# if we're running this as a separate script, we need to parse arguments and then call gen_video, that's pretty much it!
if __name__ == "__main__":
    args = parsing_video_args(sys.argv)
    folder = args.folder[0]
    exec_folder = sys.argv[0].rsplit("/", 1)[0]
    if (exec_folder.endswith(".py")):
        exec_folder = "."
    # pymol_test()
    # prepare_script(sys.argv,"t1t.mpg")
    gen_video(exec_folder, args, folder)
