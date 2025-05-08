#! /bin/sh

#SBATCH --job-name=primary_immune_response
#SBATCH --output=/home/rmorantovar/out/out_%A_%a.out
#SBATCH --error=/home/rmorantovar/out/err_%A_%a.out
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --array=0-5


env="env_roberto"
pythondir=/home/rmorantovar/anaconda3/envs/${env}/bin
srcdir=/home/rmorantovar/python/primary_immune_response
datadir=/home/rmorantovar/out

export JOBLIB_START_METHOD="forkserver"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Define your array with values
L0s=(1000 10000 100000 1000000 10000000 100000000)
css=(100 1000 10000 100000 1000000 10000000)
t_lims=(8.8 8.0 7.2 6.4 5.6 4.8 4.0)
E_lims=(-7.5 -8 -10 -11 -11 -11 -11)

L0=${L0s[${SLURM_ARRAY_TASK_ID}]}
cs=${css[${SLURM_ARRAY_TASK_ID}]}
t_lim=${t_lims[${SLURM_ARRAY_TASK_ID}]}
E_lim=${E_lims[${SLURM_ARRAY_TASK_ID}]}

echo "Start"

sleep $(($SLURM_ARRAY_TASK_ID*10))

cd $srcdir
echo "recognition_for_server.py --N_ens 2 --L0 $L0 --t_lim $t_lim --E_lim $E_lim --E_m -24 --chunk_size $cs --p 3.0 --k_step 720 --n_jobs 50 --seqs 1"
python recognition_for_server.py --N_ens 2 --L0 $L0 --t_lim $t_lim --E_lim $E_lim --E_m -24 --chunk_size $cs --p 3.0 --k_step 720 --n_jobs 50 --seqs 1
cd ~

echo "Done"
