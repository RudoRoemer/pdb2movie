# -*- coding: utf-8 -*-

'''
runfirst.py - functions to run FIRST simulations, analyzing rigidity in proteins
'''


import os


'''
firstsim - main function for running FIRST analysis on a protein

Inputs:
- string exec_folder: folder where the python scripts are located (full path)
- struct args: structure containing all arguments already parsed
- string cleanpdb: full path to PDB file after cleaning (by extension, also includes path to where all outputs are)

Outputs:
- string hydropdb: full path to PDB file after addition of hydrogens 

'''


def firstsim(exec_folder,args,cleanpdb):
    # print "./reduce.3.23.130521 -DB reduce_het_dict.txt -build "+cleanpdb+" > "+cleanpdb[:-9]+"hydro.pdb"

    # first we isolate the name of the protein, which is in the clean PDB path somewhere!
    prot=cleanpdb[:-10].rsplit("/",1)[1]
    # print "prot is "+prot

    # now, we need to run reduce to add hydrogens to protein residues - this generates a PDB file with hydrogens
    os.system(exec_folder+"/./reduce.3.23.130521 -DB "+exec_folder+"/reduce_het_dict.txt -build "+cleanpdb+" > "+cleanpdb[:-9]+"hydro.pdb")

    # we also isolate the folder where outputs are being saved
    folder=cleanpdb.rsplit("/",1)[0]

    # now, we call an external function to renumber the atoms taking the hydrogens into account
    renum_atoms(cleanpdb[:-9]+"hydro.pdb",folder)


    # finally we run FIRST with the PDB after hydrogen addition!
    os.system(exec_folder+"/./FIRST-190916-SAW/src/FIRST "+cleanpdb[:-9]+"hydro.pdb -non -dil 1 -E -0 -covout -hbout -phout -srout")

    # housekeeping: move the output files to the output folder
    os.system("mv "+prot+"_hydro_hbdilute.txt "+folder+"/"+prot+"_hydro_hbdilute.txt")
    os.system("mv "+prot+"_hydro_hbdpath.txt "+folder+"/"+prot+"_hydro_hbdpath.txt")
    os.system("mv "+prot+"_hydro_hbd_lrc.txt "+folder+"/"+prot+"_hydro_hbd_lrc.txt")
    os.system("mv residue_map.txt "+folder+"/residue_map.txt")
    # print("\n\n\n\n"+cleanpdb[:-9]+"hydro.pdb")

    # finally, return the hydro-added PDB path
    return cleanpdb[:-9]+"hydro.pdb"
    # print folder,prot


'''
renum_atoms: renumber the atoms in a PDB file to make sure they're in order and without gaps

Inputs:
string filename: path to origin PDB file
string folder: path to output folder

'''


def renum_atoms(filename,folder):

    # open the input file and a temp file at the output folder
    inputfile=open(filename,'r')
    tempfile=open(folder+"/tmp.pdb",'w')

    # let's learn how to count!
    counter=1
    # test=format(1, '05d')
    # print test

    # looping over every line in the input file
    for line in inputfile:

        # if line is an atom...
        if (line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM'):

            # we chop up the line, cutting out the current atom number and putting the value of counter as atom number
            line=line[:6]+format(counter, '05d')+line[11:]
            # write it to temp file and increase counter
            tempfile.write(line)
            counter=counter+1

    # now we can close them both
    inputfile.close()
    tempfile.close()

    # finally, we replace the original file with the temporary one
    os.system("rm "+filename)
    os.system("mv "+folder+"/tmp.pdb "+filename)
