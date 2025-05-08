#! /bin/bash

#SBATCH --job-name=primary_immune_response
#SBATCH --output=/home/rmorantovar/out/out_%A_%a.out
#SBATCH --error=/home/rmorantovar/out/err_%A_%a.out
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --array=1,10


env="env_roberto"
pythondir=/home/rmorantovar/anaconda3/envs/${env}/bin
srcdir=/home/rmorantovar/python
datadir=/home/rmorantovar/out

export JOBLIB_START_METHOD="forkserver"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

L0=$((10000000*${SLURM_ARRAY_TASK_ID}))
cs=$((L0 / 10))
echo "Start"

sleep $(($SLURM_ARRAY_TASK_ID*10))

cd $srcdir
echo "primary_immune_response.py --N_ens 400 --L0 $L0 --t_lim 8 --E_lim -12.0 --E_m -24 --chunk_size $cs --p 3.0 --k_step 720 --n_jobs 50"
python primary_immune_response.py --N_ens 400 --L0 $L0 --t_lim 8 --E_lim -12.0 --E_m -24 --chunk_size $cs --p 3.0 --k_step 720 --n_jobs 50
cd ~

echo "Done"
