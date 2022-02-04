#!/bin/sh

#  simulation_Kitas.sh
#  
#
#  Created by Roberto Moran Tovar on 25.01.22.
L=25
T0=0
alpha=1.0
beta=0.5
Ns=(200 2000 20000 200000)
Ne=2000
GMs_names=('exponential')
GMs=(0)
EMs=('MJ')
antigen='TACNSFVMTSSEYPNPNEFYMAWYC'
#antigen='TACNSFVMTSATNLFSEYPN'
#antigen='FMLFMAVFVMTSWYC'
#antigen='TACNSEYPNTTK'
#antigen='NTKTAATNLF'

c++ ../../Cpp/Dynamics_ensemble.cpp -lgsl -lgslcblas -o ../dynamics_ensemble.x
chmod +x ../dynamics_ensemble.x
for N in "${Ns[@]}"
    do
	for i in "${!GMs[@]}"
		do
		for j in "${!EMs[@]}"
			do
			T=$(echo "scale=1; 25.0 / $alpha" | bc -l)
			printf "_____________________ \n"
			printf "Parameters: L:$L ; N:$N ; T:$T  ; T0:$T0 ; alpha:$alpha ; beta:$beta ; Ne:$Ne ; growth_model:${GMs_names[i]} ; energy_model:${EMs[j]} ; antigen:$antigen\n"
			./../dynamics_ensemble.x $L $N $T $T0 $alpha $beta 0 $Ne ${GMs[i]} ${EMs[j]} $antigen
			done
		done
	done
