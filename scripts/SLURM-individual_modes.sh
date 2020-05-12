#!/bin/bash

ecut=${1:-1.0}
mode=${2:-7}   
direct=${3:-pos}

echo "CovID: mode=" $mode ", ecut=" $ecut ", direction=" $direct

if [ -f Run-$ecut-mode$mode-$direct.mpg ]
then 
    echo "WRNG:" Run-$ecut-mode$mode-$direct.mpg " exists --- please delete before proceeding!"
    exit
fi

jobfile=`printf "$ecut-$mode-$direct.sh"`
echo $jobfile

cat > ${jobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
##SBATCH --mem-per-cpu=2012
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=03:00:00

pwd
echo "--- working on $ecut, $mode, $direct"

pymol -cq ./pymolvideo$ecut$mode$direct.py -- ./Runs/$ecut/Mode$mode-$direct/

echo "--- finished with $ecut, $mode, $direct"
EOD

#cat ${jobfile}

chmod 755 ${jobfile}
#(sbatch -q devel ${jobfile})
(sbatch -q taskfarm ${jobfile})
#(sbatch ${jobfile})
#source ${jobfile}

