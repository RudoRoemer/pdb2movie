import os
import multiprocessing



def call_froda(command):
    os.system(command)




def frodasim(args,hydro_file):
    folder=hydro_file.rsplit("/",1)[0]







    if args.confs:
        totconf=int(args.confs[0])
    else:
        totconf=5000

    if args.freq:
        freq=int(args.freq[0])
    else:
        freq=100

    if args.step:
        step=float(args.step[0])
    else:
        step=0.1

    if args.dstep:
        dstep=float(args.dstep[0])
    else:
        dstep=0.01

    if args.modes:
        modelist=[int(x) for x in args.modes]
    else:
        modelist=range(7,12)

    if args.ecuts:
        cutlist=[float(x) for x in args.ecuts]
    else:
        cutlist=[1.0, 2.0]

    prot=hydro_file.rsplit("/",1)[1][:-10]
    # print modelist, cutlist
    modelist=[format(i, '02d') for i in modelist]
    # print modelist
    signals=["pos","neg"]
    jobs=[]
    try:
        os.mkdir(folder+"/Runs/")
    except Exception:
        os.system("rm -r "+folder+"/Runs/*")
        pass
    for cut in cutlist:
        try:
            os.mkdir(folder+"/Runs/"+str(cut))
        except Exception:
            os.system("rm -r "+folder+"/Runs/"+str(cut)+"/*")
            pass
        for mode in modelist:

            for sign in signals:
                try:
                    os.mkdir(folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign)
                except Exception:
                    os.system("rm -r "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/*")
                    pass
                os.system("cp "+folder+"/Modes/"+"mode"+mode+".in "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/mode.in")
                os.system("cp "+hydro_file+" "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/tmp.pdb")
                os.system("touch "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/stacked.in")
                os.system("cp "+folder+"/cov.out "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/cov.in")
                os.system("cp "+folder+"/hbonds.out "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/hbonds.in")
                os.system("cp "+folder+"/hphobes.out "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/hphobes.in")

                if (sign=="neg"):
                    command="./FIRST-190916-SAW/src/FIRST "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/tmp.pdb"+" -non -E -"+str(cut)+" -FRODA -mobRC1 -freq "+str(freq)+" -totconf "+str(totconf)+" -modei -step "+str(step)+" -dstep -"+str(dstep)+" -covin -hbin -phin -srin"
                else:
                    command="./FIRST-190916-SAW/src/FIRST "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/tmp.pdb"+" -non -E -"+str(cut)+" -FRODA -mobRC1 -freq "+str(freq)+" -totconf "+str(totconf)+" -modei -step "+str(step)+" -dstep "+str(dstep)+" -covin -hbin -phin -srin"
                p = multiprocessing.Process(target=call_froda,args=(command,))
                jobs.append(p)
                p.start()
    for job in jobs:
        job.join()
    return
