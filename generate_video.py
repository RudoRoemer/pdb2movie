# run this script by typing pymol view.py at the command line in the directory that contains the pdb files you want to view as a movie.
# edit name appropriately - this is for the standard input_froda_* files

# this is to view all of the pdb files in your current working directory in number order, using pymol.
import sys
import os
import multiprocessing
import argparse
import subprocess


def parsing_video_args(sys_args):
    parser = argparse.ArgumentParser(description='Generates videos for the most likely movement modes given a folder where the runs are stored.', usage='%(prog)s folder [options]')

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

    # args = parser.parse_args(sys_args)
    args = parser.parse_args(sys_args[1:])
    return args


def call_pymol(command):
    os.system(command)


def gen_video(exec_folder, args, folder):
    print(folder)
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
    jobs = []

    # for cut in cutlist:
    #     for mode in modelist:
    #         for sign in signals:
    #             filename = folder + "/Run-"+str(cut)+"-mode"+mode+"-"+sign+".mpg"
    #             print (filename)
    #             prepare_script(exec_folder, args, filename, cut, mode, sign, folder)
    #             # Desired pymol commands here to produce and save figures

    #             currfolder = folder + "/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/"
    #             if args.threed:
    #                 command = 'pymol -q '+folder+'/pymolvideo'+str(cut)+mode+sign+'.py -- '+currfolder
    #             else:
    #                 command = 'pymol -cq '+folder+'/pymolvideo'+str(cut)+mode+sign+'.py -- '+currfolder
    #             p = multiprocessing.Process(target=call_pymol, args=(command,))
    #             jobs.append(p)
    #             p.start()

    # for job in jobs:
    #     job.join()

    if (os.system('grep \'FREEMOL\' $(which pymol)')):
        for cut in cutlist:
            for mode in modelist:
                for sign in signals:
                    filename = folder+"/Run-"+str(cut)+"-mode"+mode+"-"+sign+".mpg"
                    tmpfolder = filename.rsplit("/", 1)[1][:-3]

                    currfolder = folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/"
                    # command = ['convert', '-quality', ' 100', folder+'/'+tmpfolder+'tmp/*.ppm', filename]
                    command = 'ffmpeg -pattern_type glob -i '+'\"'+folder+'/'+tmpfolder+'tmp/*.ppm'+'\" -c:v mpeg2video -pix_fmt yuv420p -me_method epzs -threads 4 -r 30.000030 -g 45 -bf 2 -trellis 2 -y -b 6000k '+filename
                    # command = ['ffmpeg', '-pattern_type', 'glob', '-i', '\"'+folder+'/'+tmpfolder+'tmp/*.ppm'+'\"', '-c:v', 'mpeg2video', '-pix_fmt', 'yuv420p', '-me_method', 'epzs', '-threads', '4', '-r', '30.000030', '-g', '45', '-bf', '2', '-trellis', '2', '-y', '-b', '6000k', filename]
                    print(command)
                    os.system("bash -c '{0}'".format(command))
                    # subprocess.call(command)

    for cut in cutlist:
        for mode in modelist:
            if args.combi:
                filename = folder+"/Run-"+str(cut)+"-mode"+mode+"-"
                os.system('cat '+filename+'pos.mpg '+filename+'neg.mpg > '+filename+'combi.mpg')
                os.system('chmod 744 '+filename+'combi.mpg')
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


def prepare_script(exec_folder, args, filename, cut, mode, sign, folder):

    # string="cat video_template.py <(echo filename=\'"+filename+"\') video_minimal.py >pymolvideo.py"
    # string='cat video_template.py <(echo \"stereo anaglyph\") <(echo filename=\\"'+filename+'\\") ' +args.video[0]+' video_minimal.py > pymolvideo.py'

    string = 'cat '+exec_folder+'/video_template.py '
    if args.res:
        string = string+'<(echo \"cmd.viewport('+str(args.res[0])+','+str(args.res[1])+')\") '
    if args.threed:
        string = string+'<(echo \"cmd.set(\\"stereo_mode\\",10)\") <(echo \"cmd.stereo(\\"on\\")\") '
    string = string+'<(echo filename=\\"'+filename+'\\") '
    if args.video:
        string = string+args.video[0]+' '
    string = string+exec_folder+'/video_minimal.py > '+folder+'/pymolvideo'+str(cut)+mode+sign+'.py'
    os.system("bash -c '{0}'".format(string))


if __name__ == "__main__":
    args = parsing_video_args(sys.argv)
    folder = args.folder[0]
    exec_folder = sys.argv[0].rsplit("/", 1)[0]
    if (exec_folder.endswith(".py")):
        exec_folder = "."
    # pymol_test()
    # prepare_script(sys.argv,"t1t.mpg")
    gen_video(exec_folder, args, folder)
