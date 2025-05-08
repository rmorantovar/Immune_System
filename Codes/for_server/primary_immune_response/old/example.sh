#! /bin/sh

#SBATCH --job-name=primary_immune_response
#SBATCH --output=/home/rmorantovar/out/job_out_%j.out
#SBATCH --error=/home/rmorantovar/out/job_err_%j.out
#SBATCH --mem=25GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/rmorantovar/gsl/lib

env="env_roberto"
pythondir=/home/rmorantovar/anaconda3/envs/${env}/bin
srcdir=/home/rmorantovar/exe
datadir=/home/rmorantovar/out

echo "Start"

echo "/immune_response_binary_server.x -a 6 -b 3 -k 12 -t 0 -T 7.0 -E TCRen -C 100000000 -B 100000000 -s TACNSEYPNTTRAKCGRWYC -q 3.5 --ensemble -N 200"
${srcdir}/immune_response_binary_server.x -a 6 -b 3 -k 12 -t 0 -T 7.0 -E TCRen -C 100000000 -B 100000000 -s TACNSEYPNTTRAKCGRWYC -q 3.5 --ensemble -N 200

echo "Done"
