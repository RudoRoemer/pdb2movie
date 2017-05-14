import sys
import os
import cleanpdb
import runfirst
import runelnemo
import runfroda
import generate_video
import argparse








def parsing_args(sys_args):
    parser = argparse.ArgumentParser(description='Runs simulations and generates videos for the most likely movement modes given a PDB file.',usage='%(prog)s PDB [options]')

    parser.add_argument('--keep',  nargs="+",
                        help='List of molecules to be kept')
    parser.add_argument('--output',  nargs=1,
                        help='Output directory')
    parser.add_argument('--waters',  action='store_true',
                        help='Flag for keeping water molecules')
    parser.add_argument('--multiple',  action='store_true',
                        help='Keep multiple chains (default: uses only chain A)')
    parser.add_argument('--transp',  action='store_true',
                        help='Use transparent background for video')
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
                        help='Python file with PyMOL commands to be run before generating video')
    parser.add_argument('pdbfile', metavar='PDB', type=str, nargs=1,
                        help='Initial PDB file')

    # args = parser.parse_args(sys_args)
    args = parser.parse_args(sys_args[1:])
    return args








if __name__ == "__main__":#
    args=parsing_args(sys.argv)
    # print "test"
    # if (args.threed):
    #     userinput=raw_input("WARNING: PyMOL windows will open during generation and they won't close by themselves. You have been warned. Are you sure you want to continue? [y/n]   ")
    #     if (userinput!="y"):
    #         quit()
    exec_folder=sys.argv[0].rsplit("/",1)[0]
    if (exec_folder.endswith(".py")):
        exec_folder="."
    # print(exec_folder)
    if (args.output):

        try:
            os.mkdir(args.output[0])
        except Exception:
            userinput=raw_input("WARNING: everything in output folder will be deleted! Are you sure you want to continue? [y/n]   ")
            if (userinput=="y"):
                os.system("rm -r "+args.output[0]+"/*")
            else:
                quit()
            pass
    else:

        try:
            os.mkdir(args.pdbfile[0][:-4])
        except Exception:
            os.system("rm -r "+args.pdbfile[0][:-4]+"/*")
            pass

    clean_file=cleanpdb.cleanPDB(args)
    hydro_file=runfirst.firstsim(exec_folder,args,clean_file)


    folder=hydro_file.rsplit("/",1)[0]
    # print("\n\n\n\nHERE IS THE RESULT!!!!\n\n")
    # print(hydro_file,folder)
    runelnemo.elnemosim(exec_folder,args,hydro_file)
    runfroda.frodasim(exec_folder,args,hydro_file)
    generate_video.gen_video(exec_folder,args,folder)
