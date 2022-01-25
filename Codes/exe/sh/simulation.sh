#!/bin/sh

#  simulation_Kitas.sh
#  
#
#  Created by Roberto Moran Tovar on 25.01.22.
L=15
T0=0
alphas=(0.2 0.5 1.0 2.0 5.0)
betas=(5.0)
N=2000
Ne=2000
GMs_names=('exponential'  'linear')
GMs=(0, 1)
EMs=('MJ')
antigen='FMLFMAVFVMTSWYC'

c++ ../../Cpp/Dynamics_ensemble.cpp -lgsl -lgslcblas -o ../dynamics_ensemble_prueba.x
chmod +x ../dynamics_ensemble_prueba.x
for beta in "${betas[@]}"
    do
    for alpha in "${alphas[@]}"
    	do
    	for i in "${!GMs[@]}"
    		do
    		for j in "${!EMs[@]}"
    			do
    			T=$(echo "scale=1; 25.0 / $alpha" | bc -l)
    			printf "_____________________ \n"
    			printf "Parameters: L:$L ; N:$N ; T:$T  ; T0:$T0 ; alpha:$alpha ; beta:$beta ; Ne:$Ne ; growth_model:${GMs_names[i]} ; energy_model:${EMs[j]}\n"
    			./../dynamics_ensemble_prueba.x $L $N $T $T0 $alpha $beta 0 $Ne ${GMs[i]} ${EMs[j]} $antigen
    			done
    		done
    	done
    done