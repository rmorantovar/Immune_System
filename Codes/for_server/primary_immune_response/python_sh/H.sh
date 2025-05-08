#! /bin/sh

#SBATCH --job-name=primary_immune_response
#SBATCH --output=/home/rmorantovar/out/out_%A_%a.out
#SBATCH --error=/home/rmorantovar/out/err_%A_%a.out
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --array=2,3,4,5,6,7,8,9,10


env="env_roberto"
pythondir=/home/rmorantovar/anaconda3/envs/${env}/bin
srcdir=/home/rmorantovar/python/primary_immune_response
datadir=/home/rmorantovar/out

export JOBLIB_START_METHOD="forkserver"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

p=$(echo "scale=1; ${SLURM_ARRAY_TASK_ID} / 2" | bc)
t_lim=$(echo "scale=1; ${SLURM_ARRAY_TASK_ID} / 2 + 3.0" | bc)
E_lim=$(echo "scale=1;  -1 * ${SLURM_ARRAY_TASK_ID} / 2 - 4" | bc)

echo "Start"
echo "p=$p"
echo "t_lim=$t_lim"

sleep $(($SLURM_ARRAY_TASK_ID*10))

cd $srcdir
echo "recognition_for_server.py --N_ens 1 --L0 1000000000 --t_lim $t_lim --E_lim $E_lim --E_m -24 --chunk_size 100000000 --p $p --k_step 720 --n_jobs 50"
python recognition_for_server.py --N_ens 1 --L0 1000000000 --t_lim $t_lim --E_lim $E_lim --E_m -24 --chunk_size 100000000 --p $p --k_step 720 --n_jobs 50
cd ~

echo "Done"
