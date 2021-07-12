# -*- coding: utf-8 -*-

'''
cleanpdb.py - removes rotamers and non-protein molecules

'''


import sys
import os

import helpers
try:
    from exceptions import RuntimeError
except ImportError:
    from builtins import RuntimeError

'''
remove_rotamers: uses pymol to keep a single rotamer from the pdb file

Inputs:
string output_filename: a filename (full path) where the pdb file with rotamers is located
string exec_folder: location where the python scripts are located

Outputs:
structure args: structured object with fields corresponding to the possible parameters from command line
'''
def remove_rotamers(output_filename,exec_folder):
    # calls pymol using remove_rotamers.py - for details on this, see that file!
    os.system('pymol -qc '+exec_folder+'/remove_rotamers.py '+output_filename)


'''
cleanPDB: 

Inputs: 
argument list args: object containing all command-line arguments as parsed by pdb2movie

Outputs:
string output_filename: a filename (full path) where the clean pdb file will be located

'''
def cleanPDB(args):

    # opens the input file received as one of the arguments for reading
    inputfile=open(args.pdbfile[0],'r')

    # print which molecules are being protected from being removed by cleanPDB
    if (args.keep!=[]):
        print("Keeping the following molecules: ")
        for i in args.keep:
            print(i)

    # set output filename based on arguments received
    output_filename="./"+args.pdbfile[0].rsplit("/",1)[1][:-4]+"_clean.pdb"

    # if the cleaned file is already here we don't need to do anything
    if (os.path.isfile(output_filename)):
        print("clean file already generated: " + os.path.basename(output_filename))
        return output_filename

    # open a file to write to
    os.system("rm -f " + output_filename + "_incomplete")
    output=open(output_filename + "_incomplete",'w')

    # initialise a list of residues
    residues=[]

    # goes over every single line from tge 
    for line in inputfile:
        #print(line)
        # only looks at atoms in the PDB file
        if (line[0:6].strip()=='ATOM'):


            # if there are multiple chains and the "multiple" flag has not been set, only use the first chain
            if ((line[21]!='A' and line[21]!=' ') and not(args.multiple)):

                continue

            # add the residue number to the list of residues
            residues.append(int(line[23:26].strip()))


        # this section looks at the HETATM lines in the file
        if (line[0:6]=='HETATM'):

            # if any extraneous molecules need to be kept, check whether this HETATM is part of one of them and keep it or not
            if args.keep:

                if line[17:20].strip() in args.keep:

                    output.write(line)
                    #print("printing line" +line)
            continue

        # if you didn't hit a continue, you get here and that line is kept
        #print("printing line" +line)
        output.write(line)

    # takes only unique values of residues list
    residues=list(set(residues))

    # check whether there are missing residues by comparing the list with a range
    try:
	    for i in range(residues[0],residues[-1]):
	        if (not(i in residues)):
	            print("WARNING: residue "+str(i)+" is missing!")
    except:
	    raise RuntimeError("cleanpdb: chain A seems to be empty/non-existent! --- aborting")

    # close files
    inputfile.close()
    output.close()

    # rename file to indicate it is complete
    os.system("mv " + output_filename + "_incomplete " + output_filename)

    # calls pymol to remove remaining rotamers
    #remove_rotamers(output_filename,exec_folder)

    #done!
    return output_filename


# calling this script by itself works as if it were pdb2movie works but is inadvisable

if __name__ == "__main__":
    # parse commmand-line arguments
    args=helpers.parsing_args(sys.argv)

    # set the output folder if not specified
    if (not args.output):
        args.output = [args.pdbfile[0][:-4]]

    # change directory to the output folder
    helpers.go_to_output_folder(args)

    cleanPDB(args)
