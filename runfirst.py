# -*- coding: utf-8 -*-
'''
runfirst.py - functions to run FIRST simulations, analyzing rigidity in proteins

'''


import os


'''
firstsim: main function for running FIRST analysis on a protein

Inputs:
- string exec_folder: folder where the python scripts are located (full path)
- string cleanpdb: full path to PDB file after cleaning (by extension, also includes path to where all outputs are)

Outputs:
- string hydropdb: full path to PDB file after addition of hydrogens 

'''
def firstsim(exec_folder, cleanpdb):

    print ("---------------------------------------------------------------")
    print ("firstsim:")
    print ("----------------------------------------------------------------")

    # now, we need to run reduce to add hydrogens to protein residues - this generates a PDB file with hydrogens
    print ("---------------------------------------------------------------")
    print ("firstsim: adding hydrogens")
    print ("----------------------------------------------------------------")

    basename = cleanpdb[:-10]

    # if the output file is already here we don't need to do anything
    if (os.path.isfile(basename + "_hydro.pdb")):
        print("   hydro file already generated: " + os.path.basename(basename + "_hydro.pdb"))
    else:
        # add hydrogens
        os.system(exec_folder + "/reduce.3.23.130521 -DB " + exec_folder + "/reduce_het_dict.txt -build " + cleanpdb + " > " + basename + "_hydro_temp.pdb")

        # now, we call an external function to renumber the atoms taking the hydrogens into account
        renum_atoms(basename + "_hydro_temp.pdb")

        # rename the file to indicate it is complete (remove "temp")
        os.system("mv " + basename + "_hydro_temp.pdb " + basename + "_hydro.pdb")

    # finally we run FIRST with the PDB after hydrogen addition!
    print ("---------------------------------------------------------------")
    print ("firstsim: running FIRST with new PDB file after hydrogens added")
    print ("----------------------------------------------------------------")

    # if the output files are not here, or if FIRST was previously stopped midway through, we need to run FIRST
    if (os.path.isfile("FIRST_in_progress") or not os.path.isfile(basename + "_hydro_results.txt")):
        os.system("rm -f *hydro_* *.out *list *map*")
        os.system("touch FIRST_in_progress")
        os.system(exec_folder + "/FIRST-190916-SAW/src/FIRST " + basename + "_hydro.pdb -non -dil 1 -E -0 -covout -hbout -phout -srout -L " + exec_folder + "/FIRST-190916-SAW")
        os.system("rm FIRST_in_progress")
    else:
        print("   FIRST already run")

    # finally, return the path of the PBD with added hydrogens
    return basename + "_hydro.pdb"



'''
renum_atoms: renumber the atoms in a PDB file to make sure they're in order and without gaps

Inputs:
string filename: path to origin PDB file

'''
def renum_atoms(filename):

    print ("---------------------------------------------------------------")
    print ("renum_atoms:")
    print ("----------------------------------------------------------------")

    # open the input file and a temp file
    inputfile = open(filename, 'r')
    tempfile = open("tmp.pdb", 'w')

    counter = 1

    # looping over every line in the input file
    for line in inputfile:

        # if line is an atom...
        if (line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM'):

            # we chop up the line, cutting out the current atom number and putting the value of counter as atom number
            line = line[:6] + format(counter, '05d') + line[11:]

            # write it to temp file and increase counter
            tempfile.write(line)
            counter = counter + 1

    # now we can close them both
    inputfile.close()
    tempfile.close()

    # finally, we replace the original file with the temporary one we just made
    os.system("rm " + filename)
    os.system("mv tmp.pdb " + filename)
