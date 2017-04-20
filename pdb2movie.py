import sys
import argparse
import os
import cleanpdb
import runfirst
import runelnemo
import runfroda
import generate_video
import argparse

if __name__ == "__main__":#

    parser = argparse.ArgumentParser(description='Runs simulations and generates videos for the most likely movement modes given a PDB file.',usage='%(prog)s PDB [options]')

    parser.add_argument('--keep',  nargs="+",
                        help='List of molecules to be kept')
    parser.add_argument('--output',  nargs=1,
                        help='Output directory')
    parser.add_argument('--waters',  action='store_true',
                        help='Flag for keeping water molecules')
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
    parser.add_argument('pdbfile', metavar='PDB', type=str, nargs=1,
                        help='Initial PDB file')

    # args = parser.parse_args(sys_args)
    args = parser.parse_args(sys.argv[1:])
    if (args.output):
        try:
            os.mkdir(args.output[0])
        except Exception:
            os.system("rm -r "+args.output[0]+"/*")
            pass
    else:

        try:
            os.mkdir(args.pdbfile[0][:-4])
        except Exception:
            os.system("rm -r "+args.pdbfile[0][:-4]+"/*")
            pass

    clean_file=cleanpdb.cleanPDB(args)
    hydro_file=runfirst.firstsim(args,clean_file)
    folder=hydro_file.rsplit("/",1)[0]
    runelnemo.elnemosim(args,hydro_file)
    runfroda.frodasim(args,hydro_file)
    generate_video.gen_video(args,folder)
