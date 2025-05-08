#! /bin/sh

#SBATCH --job-name=primary_immune_response
#SBATCH --output=/home/rmorantovar/out/job_out_%A_%a.out
#SBATCH --error=/home/rmorantovar/out/job_err_%A_%a.out
#SBATCH --mem=25GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --array=1,10,100,1000,10000

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/rmorantovar/gsl/lib

env="env_roberto"
pythondir=/home/rmorantovar/anaconda3/envs/${env}/bin
srcdir=/home/rmorantovar/exe
datadir=/home/rmorantovar/out

B=$((($SLURM_ARRAY_TASK_ID)*100000))
C=$((100000000000/($SLURM_ARRAY_TASK_ID)))

echo "Start"

sleep 100

cd $srcdir
echo "./immune_response_binary_server.x -a 6 -b 3 -k 120 -t 0 -T 7.5 -E TCRen -C $C -B $B -s TACNSEYPNTTRAKCGRWYC -q 2 --ensemble -N 400"
./immune_response_binary_server.x -a 6 -b 3 -k 120 -t 0 -T 7.5 -E TCRen -C $C -B $B -s TACNSEYPNTTRAKCGRWYC -q 2 --ensemble -N 400
cd ~

echo "Done"
