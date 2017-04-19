import sys
import argparse
import os
import cleanpdb
import runfirst
import runelnemo
import runfroda
import generate_video

if __name__ == "__main__":#
    try:
        os.mkdir(sys.argv[1][:-4])
    except Exception:
        os.system("rm "+sys.argv[1][:-4]+"/*")
        pass
    cutlist=[1.0, 2.0]
    modelist=range(7,12)
    clean_file=cleanpdb.cleanPDB(sys.argv)
    hydro_file=runfirst.firstsim(sys.argv,clean_file)
    folder=hydro_file.rsplit("/",1)[0]
    runelnemo.elnemosim(sys.argv,hydro_file)
    runfroda.frodasim(sys.argv,hydro_file,cutlist,modelist)
    generate_video.gen_video(folder,cutlist,modelist)
