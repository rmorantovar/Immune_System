#!/bin/bash

# Parameter values
ps=(1.5 1.8 2.2 2.6)
E_lims=(-8.0 -8.0 -8.0 -8.0)
t_lims=(5.0 5.0 5.0 5.0)
L0s=(10000 100000 1000000 10000000)
k_steps=(102811 28586 9901 4031)
lamAs=(6.0)

# Output directory
output_dir=/Users/robertomorantovar/Dropbox/Research/Immune_system/memory/exploration/log
srcdir=/Users/robertomorantovar/Dropbox/My_Documents/Science/PhD-LaessigGroup/Immune_System/_Repository/Codes/python/memory

cd $srcdir

for lamA in "${lamAs[@]}"; do
	for index in "${!ps[@]}"; do
	    # Get parameters from both arrays using the same index
	    p="${ps[$index]}"
	    E_lim="${E_lims[$index]}"
	    t_lim="${t_lims[$index]}"
	    L0="${L0s[$index]}"
	    k_step="${k_steps[$index]}"

	    # Create a unique identifier based on parameters
	    identifier="${lamA}__${p}__${t_lim}__${E_lim}"
	    echo $identifier
	    # Run Python script using nohup
	    python recognition.py --N_ens 1 --L0 "$L0" --t_lim "$t_lim" --E_lim "$E_lim" --E_m -22 --chunk_size "$L0" --p "$p" --k_step "$k_step" --lamA "$lamA" > "$output_dir/output_$identifier.log"

	    # Optionally, add a delay between jobs to avoid overwhelming the system
	    sleep 10
	done
done