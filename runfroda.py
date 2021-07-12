'''
runfroda.py - functions to run FRODA_simulations, generating conformers for the proteins

'''


import os
from multiprocessing import Pool


'''
call_command: simple wrapper for calling shell commands

Inputs: 
- string command: shell command to be run

'''
def call_command(command):
    print ("--- calling:" + command)
    os.system(command)
    print ("--- finished calling:" + command)
    return



'''
frodasim: main function for running FRODA simulations on a protein

Inputs:
- string exec_folder: folder where the python scripts are located (full path)
- struct args: structure containing all arguments already parsed
- string hydro_file: full path to PDB file after addition of hydrogens (by extension, also includes path to where all outputs are)

'''
def frodasim(exec_folder,args,hydro_file):

    print ("---------------------------------------------------------------")
    print ("frodasim:")
    print ("----------------------------------------------------------------")

    # first things first: we get the output folder from the path to the PDB file with hydrogens
    folder=hydro_file.rsplit("/",1)[0]

    # now we need to care about a bunch of command-line arguments: if they were passed, we set them, otherwise we use default values
    # that is true for confs, freq, step, dstep, modes, ecuts
    if args.confs:
        totconf=int(args.confs[0])
    else:
        totconf=1000

    if args.freq:
        freq=int(args.freq[0])
    else:
        freq=50

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
        modelist=range(7, 9)

    if args.ecuts:
        cutlist=[float(x) for x in args.ecuts]
    else:
        cutlist=[2.0] #[1.0, 2.0]

    print ("---------------------------------------------------------------")
    print ("runfroda: preparing folders for FRODA")
    print ("----------------------------------------------------------------")

    # now we make modelist into a list of strings instead of ints and define the signs positive and negative
    modelist=[format(i, '02d') for i in modelist]
    signals=["pos","neg"]

    # create a folder (or empty an existing folder) for all the conformers and other outputs
    try:
        os.mkdir(folder+"/Runs/")
    except Exception:
        os.system("rm -r "+folder+"/Runs/*")
        pass

    # create a list of commands, that will later be given to a pool of processes to run
    commands = []

    # create and organise folders, and add commands (that will generate the conformers) to 'commands'
    for cut in cutlist:

        # folder for each cut
        try:
            os.mkdir(folder+"/Runs/"+str(cut))
        except Exception:
            os.system("rm -r "+folder+"/Runs/"+str(cut)+"/*")
            pass

        for mode in modelist:

            # subfolder for each mode
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

                # now we generate commands to run FRODA
                # note that there's a slight difference (in the dstep) between positive and negative directions
                if (sign=="neg"):
                    command=exec_folder+"/./FIRST-190916-SAW/src/FIRST "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/tmp.pdb"+" -non -E -"+str(cut)+" -FRODA -resrmsd -mobRC1 -freq "+str(freq)+" -totconf "+str(totconf)+" -modei -step "+str(step)+" -dstep -"+str(dstep)+" -covin -hbin -phin -srin -L "+exec_folder+"/FIRST-190916-SAW"
                else:
                    command=exec_folder+"/./FIRST-190916-SAW/src/FIRST "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/tmp.pdb"+" -non -E -"+str(cut)+" -FRODA -resrmsd -mobRC1 -freq "+str(freq)+" -totconf "+str(totconf)+" -modei -step "+str(step)+" -dstep "+str(dstep)+" -covin -hbin -phin -srin -L "+exec_folder+"/FIRST-190916-SAW"

                #add each command to the list of commands
                commands.append(command)

    print ("---------------------------------------------------------------")
    print ("runfroda: running FRODA to generate conformer PDBs")
    print ("----------------------------------------------------------------")
    
    # give the commands to a pool of processes (one for each core) to run
    # the 'chunksize' is 1, meaning the commands are given to processes individually as the processes become free
    Pool().map(call_command, commands, 1)

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
