import os
import sys

def elnemosim(exec_folder,args,hydropdb):
    # print(hydropdb)
    if args.modes:
        modelist=[int(x) for x in args.modes]
    else:
        modelist=range(7,12)

    struct_file=generate_structure(hydropdb)
    folder=struct_file.rsplit("/",1)[0]
    # print(struct_file)
    os.system("cp "+struct_file+" tmp_"+str(os.getpid())+".structure")
    os.system("cp "+exec_folder+"/pdbmat "+folder+"/pdbmat")
    generate_pdbmat(struct_file,folder)
    print("calling "+folder+"/pdbmat")
    os.system(folder+"/pdbmat")
    # os.system("rm "+folder+"/tmp_"+str(os.getpid())+".structure")
    # os.system("rm "+folder+"/pdbmat")
    os.system("cp "+exec_folder+"/./FIRST-190916-SAW/src/diagstd "+folder+"/diagstd")
    os.system(folder+"/diagstd")
    # os.system("rm "+folder+"/diagstd")
    os.system("cp "+exec_folder+"/modesplit "+folder+"/modesplit")
    os.system("mv pdbmat.eigenfacs "+folder+"/pdbmat.eigenfacs")
    os.system(folder+"/modesplit "+struct_file+" "+folder+"/pdbmat.eigenfacs "+str(modelist[0])+" "+str(modelist[-1]))
    # os.system("rm "+folder+"/modesplit")

    os.system("mv pdbmat.* mode*.in "+folder+"/")

    try:
        os.mkdir(folder+"/Modes/")
    except Exception:
        os.system("rm -r "+folder+"/Modes/*")
        pass
    os.system("mv "+folder+"/mode*.in "+folder+"/Modes/")
    os.system("rm tmp_"+str(os.getpid())+".structure")
    os.system("rm "+folder+"/pdbmat")
    os.system("rm "+folder+"/diagstd")
    os.system("rm "+folder+"/modesplit")
    # print folder,prot


def generate_structure(filename):

    inputfile=open(filename,'r')
    outputfile=filename[:-10]+".structure"
    print outputfile
    tempfile=open(outputfile,'w')

    for line in inputfile:

        if (line[0:6].strip()=='ATOM' and line[12:15].strip()=='CA'):
            # line=line[:7]+format(counter, '04d')+line[11:]
            tempfile.write(line)
            # counter=counter+1
    inputfile.close()
    tempfile.close()
    return outputfile
    # os.system("rm "+filename)
    # os.system("mv tmp.pdb "+filename)


def generate_pdbmat(struct_file,folder):
    # print("pdbmat.dat at "+folder+"/pdbmat.dat")
    filename='pdbmat.dat'
    # print(filename)
    datfile=open(filename,'w')
    datfile.write("Coordinate FILENAME = tmp_"+str(os.getpid())+".structure\n")
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





if __name__ == "__main__":#
    args=[]
    hydro=sys.argv[1]
    exec_folder=sys.argv[0].rsplit("/",1)[0]
    if (exec_folder.endswith(".py")):
        exec_folder="."
    # pymol_test()
    # prepare_script(sys.argv,"t1t.mpg")
    elnemosim(exec_folder,args,hydro)
