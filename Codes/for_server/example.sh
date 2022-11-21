#! /bin/sh

#SBATCH --job-name=job
#SBATCH --output=job_out_%j.out
#SBATCH --error=job_err_%j.out
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/rmorantovar/gsl/lib

env="env_roberto"
pythondir=/home/rmorantovar/anaconda3/envs/${env}/bin
srcdir=/home/rmorantovar/exe
datadir=/home/rmorantovar/out

echo "Start"

${srcdir}/primary_immune_response.x -a 6 -b 0.5 -k 1 -t 0 -T 7.5 -E TCRen -C 100000 -B 200000000 -s EYTACNSEYPNTTKCGRWYCGRYPN -q 2.5 --ensemble -N 2

echo "Done"
