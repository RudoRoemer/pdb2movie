import os

def firstsim(exec_folder,args,cleanpdb):
    # print "./reduce.3.23.130521 -DB reduce_het_dict.txt -build "+cleanpdb+" > "+cleanpdb[:-9]+"hydro.pdb"
    prot=cleanpdb[:-10].rsplit("/",1)[1]
    # print "prot is "+prot
    os.system(exec_folder+"/./reduce.3.23.130521 -DB "+exec_folder+"/reduce_het_dict.txt -build "+cleanpdb+" > "+cleanpdb[:-9]+"hydro.pdb")
    folder=cleanpdb.rsplit("/",1)[0]


    renum_atoms(cleanpdb[:-9]+"hydro.pdb",folder)

    os.system(exec_folder+"/./FIRST-190916-SAW/src/FIRST "+cleanpdb[:-9]+"hydro.pdb -non -dil 1 -E -0 -covout -hbout -phout -srout")


    os.system("mv "+prot+"_hydro_hbdilute.txt "+folder+"/"+prot+"_hydro_hbdilute.txt")
    os.system("mv "+prot+"_hydro_hbdpath.txt "+folder+"/"+prot+"_hydro_hbdpath.txt")
    os.system("mv "+prot+"_hydro_hbd_lrc.txt "+folder+"/"+prot+"_hydro_hbd_lrc.txt")
    os.system("mv residue_map.txt "+folder+"/residue_map.txt")
    # print("\n\n\n\n"+cleanpdb[:-9]+"hydro.pdb")
    return cleanpdb[:-9]+"hydro.pdb"
    # print folder,prot


def renum_atoms(filename,folder):

    inputfile=open(filename,'r')
    tempfile=open(folder+"/tmp.pdb",'w')
    counter=1
    # test=format(1, '05d')
    # print test
    for line in inputfile:

        if (line[0:6].strip()=='ATOM'):
            line=line[:6]+format(counter, '05d')+line[11:]
            tempfile.write(line)
            counter=counter+1
    inputfile.close()
    tempfile.close()
    os.system("rm "+filename)
    os.system("mv "+folder+"/tmp.pdb "+filename)
