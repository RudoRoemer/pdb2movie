#!/bin/bash

dir=$1
pdbs=$2
options=$3
threads=${4:-1}

for pdb in $pdbs
do

jobfile=`printf "$pdb.sh"`

cat > ${jobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$threads
#SBATCH --mem-per-cpu=2012
#SBATCH --time=48:00:00

python3 ../pdb2movie.py $dir$pdb $options

EOD

(sbatch -q taskfarm ${jobfile})

done
