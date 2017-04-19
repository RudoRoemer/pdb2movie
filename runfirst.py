import os

def firstsim(sys_args,cleanpdb):
    # print "./reduce.3.23.130521 -DB reduce_het_dict.txt -build "+cleanpdb+" > "+cleanpdb[:-9]+"hydro.pdb"
    os.system("./reduce.3.23.130521 -DB reduce_het_dict.txt -build "+cleanpdb+" > "+cleanpdb[:-9]+"hydro.pdb")
    renum_atoms(cleanpdb[:-9]+"hydro.pdb")
    os.system("./FIRST-190916-SAW/src/FIRST "+cleanpdb[:-9]+"hydro.pdb -non -dil 1 -E -0 -covout -hbout -phout -srout")
    folder=cleanpdb.rsplit("/",1)[0]
    prot=folder.rsplit("/",1)[1]
    os.system("mv "+prot+"_hydro_hbdilute.txt "+folder+"/"+prot+"_hydro_hbdilute.txt")
    os.system("mv "+prot+"_hydro_hbdpath.txt "+folder+"/"+prot+"_hydro_hbdpath.txt")
    os.system("mv "+prot+"_hydro_hbd_lrc.txt "+folder+"/"+prot+"_hydro_hbd_lrc.txt")
    return cleanpdb[:-9]+"hydro.pdb"
    # print folder,prot


def renum_atoms(filename):

    inputfile=open(filename,'r')
    tempfile=open("tmp.pdb",'w')
    counter=1
    # test=format(1, '05d')
    # print test
    for line in inputfile:

        if (line[0:6].strip()=='ATOM'):
            line=line[:7]+format(counter, '04d')+line[11:]
            tempfile.write(line)
            counter=counter+1
    inputfile.close()
    tempfile.close()
    os.system("rm "+filename)
    os.system("mv tmp.pdb "+filename)
