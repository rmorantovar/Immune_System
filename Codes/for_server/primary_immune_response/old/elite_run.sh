#! /bin/sh

#SBATCH --job-name=primary_immune_response
#SBATCH --output=/home/rmorantovar/out/job_out_%A_%a.out
#SBATCH --error=/home/rmorantovar/out/job_err_%A_%a.out
#SBATCH --mem=25GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --array=1-99

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/rmorantovar/gsl/lib

env="env_roberto"
pythondir=/home/rmorantovar/anaconda3/envs/${env}/bin
srcdir=/home/rmorantovar/exe
datadir=/home/rmorantovar/out

n_ens=$((400+$SLURM_ARRAY_TASK_ID))

echo "Start"

sleep $(($SLURM_ARRAY_TASK_ID*4))

cd $srcdir
echo "./immune_response_elite_binary_server.x -a 6 -b 3 -k 12 -t 0 -T 7.0 -E TCRen -C 100000000 -B 100000000 -s TACNSEYPNTTRAKCGRWYC -q 3.0 --ensemble -N $n_ens"
./immune_response_elite_binary_server.x -a 6 -b 3 -k 12 -t 0 -T 7.0 -E TCRen -C 100000000 -B 100000000 -s TACNSEYPNTTRAKCGRWYC -q 3.0 --ensemble -N $n_ens
cd ~

echo "Done"

