'''
runfroda.py - functions to run FRODA_simulations, generating conformers for the proteins

'''


import os
from multiprocessing import Pool, cpu_count


'''
call_command: simple wrapper for calling shell commands

Inputs: 
- string command: shell command to be run

'''
def call_command(command):
    print ("--- calling:" + command)
    os.system(command)
    print ("--- finished calling:" + command)
    print(os.getcwd())
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

    # now we need to care about a bunch of command-line arguments: if they were passed, we set them, otherwise we use default values
    # that is true for confs, freq, step, dstep, modes, ecuts
    if args.confs:
        totconf = int(args.confs[0])
    else:
        totconf = 1000

    if args.freq:
        freq = int(args.freq[0])
    else:
        freq = 50

    if args.step:
        step = float(args.step[0])
    else:
        step = 0.1

    if args.dstep:
        dstep = float(args.dstep[0])
    else:
        dstep = 0.01

    if args.modes:
        modelist = [int(x) for x in args.modes]
    else:
        modelist = range(7, 9)

    if args.ecuts:
        cutlist = [float(x) for x in args.ecuts]
    else:
        cutlist = [2.0] #[1.0, 2.0]

    print ("---------------------------------------------------------------")
    print ("runfroda: preparing folders and commands for FRODA")
    print ("----------------------------------------------------------------")

    # now we make modelist into a list of strings instead of ints and define the signs positive and negative
    modelist = [format(i, '02d') for i in modelist]
    signals = ["pos", "neg"]

    base = os.getcwd()

    # folder for all the conformers (and other outputs) with this step and dstep combination
    runs_dir = "Runs"
    if (args.single):
        runs_dir += "-singlechain"
    if (step != 0.1):
        runs_dir += "-step" + str(step)
    if (dstep != 0.01):
        runs_dir += "-dstep" + str(dstep)
    os.system("mkdir -p " + runs_dir)
    os.chdir(runs_dir)

    # create a list of commands, that will later be given to a pool of processes to run
    commands = []

    # create and organise folders, and add commands (that will generate the conformers) to 'commands'
    for cut in cutlist:

        # folder for each cut
        os.system("mkdir -p " + str(cut))
        os.chdir(str(cut))

        for mode in modelist:
            for sign in signals:

                # subfolder for each mode-sign combination
                os.system("mkdir -p Mode" + mode + "-" + sign)
                os.chdir("Mode" + mode + "-" + sign)

                # determine if we need to generate conformers
                # if confs is not more than previously, and freq is a multiple of previously, then we do not so we skip
                files = os.listdir(".")
                confs_file = [x for x in files if (x[:5] == "CONFS")]
                freq_file = [x for x in files if (x[:4] == "FREQ")]
                if (confs_file and freq_file):
                    if (int(confs_file[0][5:]) >= totconf):
                        if (freq % int(freq_file[0][4:]) == 0):
                            print("   conformers already generated for " + str(cut) + "-mode" + mode + "-" + sign)
                            os.chdir("..")
                            continue

                # remove previous output files
                os.system("rm -f *")
                #os.system("rm -f CONFS* FREQ* problem.log tmp_bond.txt tmp_data.txt tmp_results.txt tmp_froda_00000000.pdb tmp_RCD.*")

                # we need to some files in this folder for FRODA
                os.system("cp " + base + "/mode" + mode + ".in ./mode.in")
                os.system("cp " + base + "/" + hydro_file + " ./tmp.pdb")
                os.system("touch stacked.in")
                os.system("cp " + base + "/cov.out ./cov.in")
                os.system("cp " + base + "/hbonds.out ./hbonds.in")
                os.system("cp " + base + "/hphobes.out ./hphobes.in")

                # now we generate commands to run FRODA
                # note that there's a slight difference (in the dstep) between positive and negative directions
                command = "cd " + os.path.abspath(".") + "\n"
                if (sign == "neg"):
                    command += exec_folder + "/./FIRST-190916-SAW/src/FIRST tmp.pdb -non -E -" + str(cut) + " -FRODA -resrmsd -mobRC1 -freq " + str(freq) + " -totconf " + str(totconf) + " -modei -step " + str(step) + " -dstep -" + str(dstep) + " -covin -hbin -phin -srin -L " + exec_folder + "/FIRST-190916-SAW"
                else:
                    command += exec_folder + "/./FIRST-190916-SAW/src/FIRST tmp.pdb -non -E -" + str(cut) + " -FRODA -resrmsd -mobRC1 -freq " + str(freq) + " -totconf " + str(totconf) + " -modei -step " + str(step) + " -dstep " + str(dstep) + " -covin -hbin -phin -srin -L " + exec_folder + "/FIRST-190916-SAW"

                # record confs and freq when finished so we know we have already done them if we run this again
                command += "\n touch CONFS" + str(totconf)
                command += "\n touch FREQ" + str(freq)

                # add each command to the list of commands
                commands.append(command)

                os.chdir("..")

        os.chdir("..")

    print ("---------------------------------------------------------------")
    print ("runfroda: running FRODA commands to generate conformer PDBs")
    print ("----------------------------------------------------------------")

    if (commands == []):
        print("   no commands to run; conformers all already generated")
    else:
        # set the number of threads to be the number of cores if not already specified
        if args.frodathreads:
            frodathreads = int(args.frodathreads[0])
        else:
            frodathreads = cpu_count()

        # give the commands to a pool of processes to run
        # the 'chunksize' is 1, meaning the commands are given to processes individually as the processes become free
        Pool(processes=frodathreads).map(call_command, commands, 1)

    # now some housekeeping: we remove all temp files we created at each subfolder
    for cut in cutlist:
        for mode in modelist:
            for sign in signals:
                os.system('rm -f ' + str(cut) + "/Mode" + mode + "-" + sign + '/hphobes.in')
                os.system('rm -f ' + str(cut) + "/Mode" + mode + "-" + sign + '/hbonds.in')
                os.system('rm -f ' + str(cut) + "/Mode" + mode + "-" + sign + '/cov.in')
                os.system('rm -f ' + str(cut) + "/Mode" + mode + "-" + sign + '/stacked.in')
                os.system('rm -f ' + str(cut) + "/Mode" + mode + "-" + sign + '/tmp.pdb')

    os.chdir("..")

    return
