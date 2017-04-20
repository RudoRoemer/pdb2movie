import os

def frodasim(args,hydro_file):
    folder=hydro_file.rsplit("/",1)[0]
    os.system("touch "+folder+"/stacked.in")
    os.system("cp "+folder+"/cov.out "+folder+"/cov.in")
    os.system("cp "+folder+"/hbonds.out "+folder+"/hbonds.in")
    os.system("cp "+folder+"/hphobes.out "+folder+"/hphobes.in")


    step=0.1
    dstep=0.01
    cutlist=[1.0, 2.0]
    modelist=range(7,12)

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

    

    modelist=[format(i, '02d') for i in modelist]
    # print modelist
    signals=["pos","neg"]
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
            os.system("cp "+folder+"/Modes/"+"mode"+mode+".in "+folder+"/mode.in")
            for sign in signals:
                try:
                    os.mkdir(folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign)
                except Exception:
                    os.system("rm -r "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/*")
                    pass

                if (sign=="neg"):
                    os.system("./FIRST-190916-SAW/src/FIRST "+hydro_file+" -non -E -"+str(cut)+" -FRODA -mobRC1 -freq "+str(freq)+" -totconf "+str(totconf)+" -modei -step "+str(step)+" -dstep -"+str(dstep)+" -covin -hbin -phin -srin")
                else:
                    os.system("./FIRST-190916-SAW/src/FIRST "+hydro_file+" -non -E -"+str(cut)+" -FRODA -mobRC1 -freq "+str(freq)+" -totconf "+str(totconf)+" -modei -step "+str(step)+" -dstep "+str(dstep)+" -covin -hbin -phin -srin")
                os.system("mv "+folder+"/R504C_PaPBP3_hydro_froda* "+folder+"/Runs/"+str(cut)+"/Mode"+mode+"-"+sign+"/")
