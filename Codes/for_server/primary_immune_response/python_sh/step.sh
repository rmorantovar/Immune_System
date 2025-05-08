#! /bin/bash

#SBATCH --job-name=primary_immune_response
#SBATCH --output=/home/rmorantovar/out/out_%A_%a.out
#SBATCH --error=/home/rmorantovar/out/err_%A_%a.out
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --array=0


env="env_roberto"
pythondir=/home/rmorantovar/anaconda3/envs/${env}/bin
srcdir=/home/rmorantovar/python
datadir=/home/rmorantovar/out

export JOBLIB_START_METHOD="forkserver"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Define your array with values
k_steps=(1. 3. 6. 12. 30. 60. 120. 300. 600. 1200.)
t_lims=(8.5 8.0 7.5 7.0 6.5 6.0 5.5 5.0 4.5 4.0)

k_step_h=${k_steps[${SLURM_ARRAY_TASK_ID}]}
k_step=$(echo "$k_step_h * 24" | bc -l)
t_lim=${t_lims[${SLURM_ARRAY_TASK_ID}]}

echo "Start"

sleep $(($SLURM_ARRAY_TASK_ID*10))

cd $srcdir
echo "primary_immune_response.py --N_ens 400 --L0 1000000000 --t_lim $t_lim --E_lim -12.0 --E_m -24 --chunk_size 100000000 --p 3.0 --k_step $k_step --n_jobs 50"
python primary_immune_response.py --N_ens 400 --L0 1000000000 --t_lim $t_lim --E_lim -12.0 --E_m -24 --chunk_size 100000000 --p 3.0 --k_step $k_step --n_jobs 50
cd ~

echo "Done"
