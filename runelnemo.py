'''
runelnemo.py - functions related to running Elastic Network model ElNemo

'''


import os
import sys

import helpers


'''
elnemosim: main function for running ElNemo simulations

Inputs:
- string exec_folder: folder where the python scripts are located (full path)
- struct args: structure containing all arguments already parsed
- string hydropdb: full path to PDB file after addition of hydrogens (by extension, also includes path to where all outputs are)

'''
def elnemosim(exec_folder, args, hydropdb):

    print ("---------------------------------------------------------------")
    print ("elnemosim:")
    print ("----------------------------------------------------------------")

     # check for a list of modes in the arguments, fills it in with defaults if no specification
    if args.modes:
        modelist=[int(x) for x in args.modes]
    else:
        modelist=range(7, 9)

    # call external function to generate a structure file - more details in that function
    print ("---------------------------------------------------------------")
    print ("elnemosim: generate structure file")
    print ("----------------------------------------------------------------")

    if (os.path.isfile(hydropdb[:-10] + ".structure")):
        print("   structure file already generated: " + os.path.basename(hydropdb)[:-10] + ".structure")
        struct_file=hydropdb[:-10] + ".structure"
    else:
        struct_file=generate_structure(hydropdb)

    # now we have a series of shell commands to set up input files for the ElNemo tools

    # now, we need to generate a pdbmat options file using an external function as well - more details there!
    print ("---------------------------------------------------------------")
    print ("elnemosim: generate pdbmat options file")
    print ("----------------------------------------------------------------")

    if (os.path.isfile("pdbmat.dat")):
        print("   pdbmat options file already generated: pdbmat.dat")
    else:
        generate_pdbmat(struct_file)

    # finally, we run pdbmat
    print ("---------------------------------------------------------------")
    print ("elnemosim: running local pdbmat()")
    print ("----------------------------------------------------------------")

    if (os.path.isfile("pdbmat_in_progress") or not os.path.isfile("pdbmat.dat_run")):
        os.system("touch pdbmat_in_progress")
        os.system(exec_folder+"/pdbmat")
        os.system("rm pdbmat_in_progress")
    else:
        print("   pdbmat already run")

    # now, we run diagstd
    print ("---------------------------------------------------------------")
    print ("elnemosim: running local diagstd()")
    print ("----------------------------------------------------------------")

    if (os.path.isfile("diagstd_in_progress") or not os.path.isfile("pdbmat.eigenfacs")):
        os.system("touch diagstd_in_progress")
        os.system(exec_folder + "/FIRST-190916-SAW/src/diagstd")
        os.system("rm diagstd_in_progress")
    else:
        print("   diagstd already run")


    # these variables store the lowest and highest mode we must run modesplit with
    # the code here makes this range as tight as possible, accounting for mode[num].in files already generated
    min_mode = modelist[0]
    max_mode = modelist[-1]

    for _, _, files in os.walk("."):
        numbers=[]
        for file in files:
            if (file[:4]=="mode" and file[-3:]==".in"):
                numbers.append(file[4:-3])
        numbers.sort()
        up_down_numbers = numbers.copy()
        numbers.reverse()
        up_down_numbers += numbers
        for number in up_down_numbers:
            if (int(number)==min_mode):
                min_mode = min_mode + 1
            if (int(number)==max_mode):
                max_mode = max_mode - 1
            if (min_mode > max_mode):
                return
        break

    # now, we run modesplit with the generated structure file, the output of diagstd and the list of modes we generated (actually, a range from lowest mode to highest mode with everything in between)
    os.system(exec_folder+"/modesplit "+struct_file+" "+"pdbmat.eigenfacs "+str(min_mode)+" "+str(max_mode))




'''
generate_structure: generate a .structure file from a PDB

Inputs: 
string filename: PDB file to be transformed into a .structure

Outputs: 
string outputfile: path to the generated .structure file

'''
def generate_structure(filename):

    # let's open up the PDB file and generate a filename for the output file based on its path
    inputfile=open(filename,'r')

    # now we open the output file as well
    outputfile=open(filename[:-10] + ".structure", 'w')

    # now we can loop over all lines from the PDB file...
    for line in inputfile:

        # ...and write only the alpha carbons to the output file
        if (line[0:6].strip()=='ATOM' and line[12:15].strip()=='CA'):
            outputfile.write(line)

    # finished with the files
    inputfile.close()
    outputfile.close()

    # rename the file to indicate it is complete (remove "temp")
    os.system("mv " + filename[:-10] + ".structure_temp " + filename[:-10] + ".structure")

    return filename[:-10] + ".structure"



'''
generate_pdbmat: a simple wrapper to a bunch of write calls to generate a very specific pdbmat.dat file with the correct options for diagstd

Inputs:
string struct_file: path to the .structure file that was already generated at this point

'''
def generate_pdbmat(struct_file):

    # we start by creating a pdbmat.dat file and opening it
    datfile=open('pdbmat_temp.dat','w')

    # the only variable in this is the name of the struct file
    datfile.write("Coordinate FILENAME = " + os.path.basename(struct_file) + "\n")
    datfile.write("INTERACtion DISTance CUTOF =     12.000\n")
    datfile.write("INTERACtion FORCE CONStant =      1.000\n")
    datfile.write("Origin of MASS values      =       CONS ! CONstant, or from COOrdinate file.\n")
    datfile.write("Output PRINTing level      =          0 ! =1: more detailled. =2: debug level.\n")
    datfile.write("Bond DEFINITION            =       NONE ! NONe, ALL, or between CONsecutive atoms.\n")
    datfile.write("Maximum bond LENGTH        =      0.000\n")
    datfile.write("BOND FORCE CONStant        =      0.000\n")
    datfile.write("ANGLE FORCE CONStant       =      0.000\n")
    datfile.write("LevelSHIFT                 =    1.0E-09 ! non-zero value often required (numerical reasons).\n")
    datfile.write("Matrix FORMAT              =       FREE ! Free, or Binary, matrix saved.\n")
    datfile.close()

    # rename the file to indicate it is complete (remove "temp")
    os.system("mv pdbmat_temp.dat pdbmat.dat")



# calling this script by itself works as if it were pdb2movie but is inadvisable

if __name__ == "__main__":
    # parse commmand-line arguments
    args=helpers.parsing_args(sys.argv)

    hydro="./" + os.path.basename(sys.argv[1])[:-4] + "_hydro" + os.path.basename(sys.argv[1])[-4:]

    # set exec_folder to the full path of this script
    exec_folder=os.path.dirname(os.path.abspath(sys.argv[0]))

    # set the output folder if not specified
    if (not args.output):
        args.output = [args.pdbfile[0][:-4]]

    # change directory to the output folder
    helpers.go_to_output_folder(args)

    elnemosim(exec_folder, args, hydro)
