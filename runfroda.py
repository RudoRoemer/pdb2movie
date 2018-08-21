'''
runfroda.py - functions to run FRODA_simulations, generating conformers for the proteins
'''

import os
import multiprocessing


'''
call_froda: simple wrapper for calling sheel commands

Inputs: 
- string command: shell command to be run
'''


def call_froda(command):
    os.system(command)



'''
frodasim: main function for running FRODA simulations on a protein

Inputs:
- string exec_folder: folder where the python scripts are located (full path)
- struct args: structure containing all arguments already parsed
- string hydro_file: full path to PDB file after addition of hydrogens (by extension, also includes path to where all outputs are)

'''



def frodasim(exec_folder,args,hydro_file):


    # first things first: we get the output folder from the path to the PDB file with hydrogens
    folder=hydro_file.rsplit("/",1)[0]






    # now we need to care about a bunch of command-line arguments: if they were passed, we set them, otherwise we use default values
    # that is true for confs, frew, step, dstep, modes, ecuts
    if args.confs:
        totconf=int(args.confs[0])
    else:
        totconf=5000

    if args.freq:
        freq=int(args.freq[0])
    else:
        freq=100

    if args.step:
        step=float(args.step[0])
    else:
        step=0.1

    if args.dstep:
        dstep=float(args.dstep[0])
    else:
        dstep=0.01

    if args.modes:
        modelist=[int(x) for x in args.modes]
    else:
        modelist=range(7,12)

    if args.ecuts:
        cutlist=[float(x) for x in args.ecuts]
    else:
        cutlist=[1.0, 2.0]


    # we isolate the name of the protein, which is in the PDB_path somewhere!
    prot=hydro_file.rsplit("/",1)[1][:-10]
    # print modelist, cutlist

    # now we make modelist into a list of strings instead of ints and define the signs positive and negative
    modelist=[format(i, '02d') for i in modelist]
    signals=["pos","neg"]

    # initialising job list for multiprocessing
    jobs=[]

    # now we need to create the folder where all conformers are going to live (or empty and existing one)
    try:
        os.mkdir(folder+"/Runs/")
    except Exception:
        os.system("rm -r "+folder+"/Runs/*")
        pass


    # now for every cutoff energy we will create a subfolder before anything
    for cut in cutlist:
        try:
            os.mkdir(folder+"/Runs/"+str(cut))
        except Exception:
            os.system("rm -r "+folder+"/Runs/"+str(cut)+"/*")
            pass

        # and for every mode and sign combination we will, before anything, create a subfolder for that
        for mode in modelist:

            for sign in signals:
                try:
                    os.mkdir(folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign)
                except Exception:
                    os.system("rm -r "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/*")
                    pass


                # now there's some setup before starting FRODA: we need to copy the correct mode file from Modes...
                os.system("cp "+folder+"/Modes/"+"mode"+mode+".in "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/mode.in")
                # ... the relevant PDB_file, and cov.out, and hbonds.out, and hphobes.out as input files (also we need to create an empty stacked.in)
                os.system("cp "+hydro_file+" "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/tmp.pdb")
                os.system("touch "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/stacked.in")
                os.system("cp "+folder+"/cov.out "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/cov.in")
                os.system("cp "+folder+"/hbonds.out "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/hbonds.in")
                os.system("cp "+folder+"/hphobes.out "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/hphobes.in")


                # now we actually run FRODA - there's a slight difference between positive and negative directions
                if (sign=="neg"):
                    command=exec_folder+"/./FIRST-190916-SAW/src/FIRST "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/tmp.pdb"+" -non -E -"+str(cut)+" -FRODA -mobRC1 -freq "+str(freq)+" -totconf "+str(totconf)+" -modei -step "+str(step)+" -dstep -"+str(dstep)+" -covin -hbin -phin -srin"
                else:
                    command=exec_folder+"/./FIRST-190916-SAW/src/FIRST "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/tmp.pdb"+" -non -E -"+str(cut)+" -FRODA -mobRC1 -freq "+str(freq)+" -totconf "+str(totconf)+" -modei -step "+str(step)+" -dstep "+str(dstep)+" -covin -hbin -phin -srin"

                # we assign a new process to this call of froda (the reason why call_froda exists at all), append it to job list and start it
                p = multiprocessing.Process(target=call_froda,args=(command,))
                jobs.append(p)
                p.start()

    # this is here so that we wait for all processes to finish before proceeding
    for job in jobs:
        job.join()


    # now some housekeeping: we remove all temp files we created at each subfolder
    for cut in cutlist:
        for mode in modelist:
            for sign in signals:
                os.system('rm '+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+'/hphobes.in')
                os.system('rm '+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+'/hbonds.in')
                os.system('rm '+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+'/cov.in')
                os.system('rm '+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+'/stacked.in')
                os.system('rm '+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+'/tmp.pdb')

    return
