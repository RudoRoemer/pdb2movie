#!/bin/bash

dir=${1:-../Covid-19}   
pdbname=${2:-test}
confs=${3:-10}
modes=${4:-"7"}
ecuts=${5:-"1.0"}
options=${6:""}

echo "CodID: dir=" $dir ", confs=" $confs ", modes=" $modes ", Ecuts=" $ecuts

for pdb in $pdbname #\
#5r80 5r81 5r82 5r83 5r84 \
#6lu7 6lvn 6lxt \
#6m03 6m17 6m3m \
#6vsb 6vw1 6vww 6vyo \
#6w01 6w02 6w4b \
#6y84 \
#wuhan_spike wuhan_E_protein_5x29 wuhan_ace2_min \
#wuhan_ace2_MD
do
#echo $pdb

jobfile=`printf "$pdb.sh"`
echo $jobfile

cat > ${jobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
##SBATCH --mem-per-cpu=2012
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=48:00:00

pwd
echo "--- working on $pdb"

cd $dir
pwd

echo $pdb-`echo $modes | sed "s/ /-/g"`_`echo $ecuts | sed "s/ /-/g"`

mkdir -p $pdb-`echo $modes | sed "s/ /-/g"`_`echo $ecuts | sed "s/ /-/g"`
cd $pdb-`echo $modes | sed "s/ /-/g"`_`echo $ecuts | sed "s/ /-/g"`
pwd

cp ../$pdb.pdb .

git clone https://github.com/RudoRoemer/pdb2movie.git
cd pdb2movie
git checkout devel
git pull

cd FIRST-190916-SAW/src
make
cd ../../

pwd

echo "python pdb2movie.py ../$pdb.pdb --combi --confs $confs --freq 50 --modes $modes --ecuts $ecuts --res 1920 1080 --multiple $options >& ../$pdb.log"
python pdb2movie.py ../$pdb.pdb --combi --confs $confs --freq 50 --modes $modes --ecuts $ecuts --res 1920 1080 --multiple $options >& ../$pdb.log 

echo "--- finished with $pdb"
EOD

#cat ${jobfile}

chmod 755 ${jobfile}
#(sbatch -q devel ${jobfile})
(sbatch -q taskfarm ${jobfile})
#(sbatch ${jobfile})
#source ${jobfile}

done
