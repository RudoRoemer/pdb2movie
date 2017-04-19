import sys
import argparse

def cleanPDB(sys_args):
# def main():

    parser = argparse.ArgumentParser(description='Clean PDB files for ElNemo.',usage='%(prog)s PDB [options]')

    parser.add_argument('--keep',  nargs="*",
                        help='List of molecules to be kept')
    parser.add_argument('--output',  nargs=1,
                        help='Output file')
    parser.add_argument('--waters',  action='store_true',
                        help='Flag for keeping water molecules')
    parser.add_argument('pdbfile', metavar='PDB', type=str, nargs=1,
                        help='PDB file to be cleaned')

    # args = parser.parse_args(sys_args)
    args = parser.parse_args(sys_args[1:])
    if (args.keep==None):
        args.keep=[]
    if args.waters:
        args.keep.append('HOH')
    inputfile=open(args.pdbfile[0],'r')
    if (args.keep!=[]):
        print("Keeping the following molecules: ")
        for i in args.keep:
            print(i)
        # print("\n")
    residues=[]
    # lines=open(args.pdbfile[0],'r')
    if (args.output):
        output_filename=args.output[0]
    else:
        folder=args.pdbfile[0].split("/")
        # print folder
        # print folder[-1][:-4]
        output_filename=args.pdbfile[0][:-4]+"/"+folder[-1][:-4]+"_clean.pdb"
    output=open(output_filename,'w')

    for line in inputfile:


        if (line[0:6].strip()=='ATOM'):
            if (line[21]!='A'):
                # print line[21]
                continue
            residues.append(int(line[23:26].strip()))
        # print line[0:6]
        if (line[0:6]=='HETATM'):
            if args.keep:
                if line[17:20].strip() in args.keep:
                    output.write(line)
            continue
        output.write(line)
    residues=list(set(residues))
    for i in range(residues[0],residues[-1]):
        if (not(i in residues)):
            print("WARNING: residue "+str(i)+" is missing!")
    inputfile.close()
    output.close()
    return output_filename
#
if __name__ == "__main__":#
    print sys.argv
    cleanPDB(sys.argv)
