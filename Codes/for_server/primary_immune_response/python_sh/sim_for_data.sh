#!/bin/bash

# Parameter values
#ps=(1 1.5 2 2.5 3 3.5 4 4.5 5)
ps=(1 1.5 2 2.5)
# E_lims=(-8.5 -8.5 -12 -12 -12 -12 -12 -12 -12)
E_lims=(-8.5 -8.5 -12 -12)
# t_lims=(4.0 4.5 5.2 5.7 6.2 6.7 7.2 7.7 8.2)
# t_lims=(5.0 5.5 6.2 6.7 7.2 7.7 8.2 8.7 9.2)
t_lims=(5.0 5.5 6.2 6.7)
# lamAs=(5.5 7.5 9.0)
lamAs=(7.5)

# Output directory
output_dir=/Users/robertomorantovar/Dropbox/Research/Immune_system/primary_immune_response/simulations_lambda_A/log
srcdir=/Users/robertomorantovar/Dropbox/My_Documents/SCIENCE/PhD-LaessigGroup/Immune_System/_Repository/Codes/python/primary_immune_response

cd $srcdir

for lamA in "${lamAs[@]}"; do
	for index in "${!ps[@]}"; do
	    # Get parameters from both arrays using the same index
	    p="${ps[$index]}"
	    E_lim="${E_lims[$index]}"
	    t_lim="${t_lims[$index]}"

	    # Create a unique identifier based on parameters
	    identifier="${lamA}__${p}__${t_lim}__${E_lim}"
	    echo $identifier
	    # Run Python script using nohup
	    python recognition.py --N_ens 100 --L0 100000000 --t_lim "$t_lim" --E_lim "$E_lim" --E_m -24 --chunk_size 10000000 --p "$p" --k_step 720 --lamA "$lamA" > "$output_dir/output_$identifier.log"

	    # Optionally, add a delay between jobs to avoid overwhelming the system
	    sleep 20
	done
done