'''
runelnemo.py - functions related to running Elastic Network model ElNemo

'''



import os
import sys

'''
elnemosim: main function for running ElNemo simulations

Inputs:
- string exec_folder: folder where the python scripts are located (full path)
- struct args: structure containing all arguments already parsed
- string hydropdb: full path to PDB file after addition of hydrogens (by extension, also includes path to where all outputs are)

'''

def elnemosim(exec_folder,args,hydropdb):

    print ("---------------------------------------------------------------")
    print ("elnemosim:")
    print ("----------------------------------------------------------------")

     # check for a list of modes in the arguments, fills it in with defaults if no specification
    if args.modes:
        modelist=[int(x) for x in args.modes]
    else:
        modelist=range(7,12)

    # call external function to generate a structure file - more details in that function
    print ("---------------------------------------------------------------")
    print ("elnemosim: generate structure file")
    print ("----------------------------------------------------------------")

    struct_file=generate_structure(hydropdb)
    # from that filename, get the folder where all the outputs are
    folder=struct_file.rsplit("/",1)[0]
    # print(struct_file)


    # now we have a series of shell commands to set up input files for the ElNemo tools


    # now, we need to generate a pdbmat options file using an external function as well - more details there!
    print ("---------------------------------------------------------------")
    print ("elnemosim: generate pdbmat options file")
    print ("----------------------------------------------------------------")

    generate_pdbmat(struct_file,folder)
    print("calling "+exec_folder+"/pdbmat")

    # finally, we run pdbmat
    print ("---------------------------------------------------------------")
    print ("elnemosim: running local pdbmat()")
    print ("----------------------------------------------------------------")

    os.system(exec_folder+"/pdbmat")

    # now, we run diagstd
    print ("---------------------------------------------------------------")
    print ("elnemosim: running local diagstd()")
    print ("----------------------------------------------------------------")

    os.system(exec_folder + "/FIRST-190916-SAW/src/diagstd")

    # now, we run modesplit with the generated structure file, the output of diagstd and the list of modes we generated (actually, a range from lowest mode to highest mode with everything in between)
    os.system(exec_folder+"/modesplit "+struct_file+" "+folder+"/pdbmat.eigenfacs "+str(modelist[0])+" "+str(modelist[-1]))

    # now we have some housekeeping to do - putting all movement modes into a "Modes" folder
    try:
        os.mkdir(folder+"/Modes/")
    except Exception:
        os.system("rm -r "+folder+"/Modes/*")
        pass
    os.system("mv "+folder+"/mode*.in "+folder+"/Modes/")



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
    outputfile=filename[:-10]+".structure"
    #print outputfile

    # now we open the output file as well
    tempfile=open(outputfile,'w')

    # now we can loop over all lines from the PDB file...
    for line in inputfile:

        # ...and write only the alpha carbons to the output file
        if (line[0:6].strip()=='ATOM' and line[12:15].strip()=='CA'):
            # line=line[:7]+format(counter, '04d')+line[11:]
            tempfile.write(line)
            # counter=counter+1

    # we close both files and return the path to the output file
    inputfile.close()
    tempfile.close()
    return outputfile
    # os.system("rm "+filename)
    # os.system("mv tmp.pdb "+filename)


'''
generate_pdbmat: a simple wrapper to a bunch of write calls to generate a very specific pdbmat.dat file with the correct options for diagstd

Inputs:
string struct_file: path to the .structure file that was already generated at this point
string folder: path to the folder where outputs are being written

'''


def generate_pdbmat(struct_file,folder):
    # print("pdbmat.dat at "+folder+"/pdbmat.dat")
    filename=folder + '/pdbmat.dat'
    # we start by creating a pdbmat.dat file and opening it
    datfile=open(filename,'w')

    # the only variable in this is the name of the struct file
    datfile.write("Coordinate FILENAME = " + struct_file + "\n")
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





# calling this script by itself is probably something you shouldn't do! But you can pass a filename for a PDB file with hydrogens and it should work, I guess?

if __name__ == "__main__":#
    args=[]
    hydro=sys.argv[1]
    exec_folder=sys.argv[0].rsplit("/",1)[0]
    if (exec_folder.endswith(".py")):
        exec_folder="."
    # pymol_test()
    # prepare_script(sys.argv,"t1t.mpg")
    elnemosim(exec_folder,args,hydro)
