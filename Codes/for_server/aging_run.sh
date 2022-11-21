#! /bin/sh

#SBATCH --job-name=primary_immune_response
#SBATCH --output=/home/rmorantovar/out/job_out_%A_%a.out
#SBATCH --error=/home/rmorantovar/out/job_err_%A_%a.out
#SBATCH --mem=25GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --array=1,2,5,10

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/rmorantovar/gsl/lib

env="env_roberto"
pythondir=/home/rmorantovar/anaconda3/envs/${env}/bin
srcdir=/home/rmorantovar/exe
datadir=/home/rmorantovar/out

B=$((10000000*$SLURM_ARRAY_TASK_ID))

echo "Start"

cd $srcdir
echo "./primary_immune_response.x -a 6 -b 3 -k 12 -t 0 -T 7.5 -E TCRen -C 100000000 -B $B -s TACNSEYPNTTRAKCGRWYC -q 3.0 --ensemble -N 500"
./primary_immune_response.x -a 6 -b 3 -k 12 -t 0 -T 7.5 -E TCRen -C 100000000 -B $B -s TACNSEYPNTTRAKCGRWYC -q 3.0 --ensemble -N 500
cd ~

echo "Done"

