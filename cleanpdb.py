'''
cleanpdb.py - removes rotamers and non-protein molecules

'''


import sys
import os


'''
remove_rotamers: uses pymol to keep a single rotamer from the pdb file

Inputs:
string output_filename: a filename (full path) where the pdb file with rotamers is located
string exec_folder: location where the python scripts are located


'''
def remove_rotamers(output_filename,exec_folder):
    # calls pymol using remove_rotamers.py - for details on this, see that file!
    os.system('pymol -qc '+exec_folder+'/remove_rotamers.py '+output_filename)




'''
cleanPDB: 

Inputs: 
- argument list args: object containing all command-line arguments as parsed by pdb2movie
- string exec_folder: location where the python scripts are located

Outputs:
- string output_filename: a filename (full path) where the clean pdb file will be located

'''
def cleanPDB(args,exec_folder):
# def main():

    # set args.keep as a list of molecules that need to be kept, based on the arguments received
    if (args.keep==None):
        args.keep=[]
    if args.waters:
        args.keep.append('HOH')
    if (args.keep!=[]):
        print("Keeping the following molecules: ")
        for i in args.keep:
            print(i)

    # opens the input file received as one of the arguments for reading
    inputfile=open(args.pdbfile[0],'r')


    # initialise a list of residues
    residues=[]

    # set output filename based on arguments received and open that file for writing
    if (args.output):
        output_filename=args.output[0]+"/"+args.pdbfile[0].rsplit("/",1)[1][:-4]+"_clean.pdb"

    else:
        folder=args.pdbfile[0].split("/")
        output_filename=args.pdbfile[0][:-4]+"/"+folder[-1][:-4]+"_clean.pdb"
    output=open(output_filename,'w')



    # goes over every single line from tge 
    for line in inputfile:

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
            continue

        # if you didn't hit a continue, you get here and that line is kept
        output.write(line)

    # takes only unique values of residues list
    residues=list(set(residues))


    #check whether there are missing residues by comparing the list with a range
    for i in range(residues[0],residues[-1]):
        if (not(i in residues)):
            print("WARNING: residue "+str(i)+" is missing!")

    # close files
    inputfile.close()
    output.close()

    # calls pymol to remove remaining rotamers
    remove_rotamers(output_filename,exec_folder)

    #done!
    return output_filename


# calling this as a single script will probably not work, I think?
#if __name__ == "__main__":#
#    cleanPDB(sys.argv)
