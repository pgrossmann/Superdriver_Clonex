#!/bin/bash

# mkdir -p "/links/shared/patrick_clonex/simulation_out/finalRuns"

SIMDIR="/links/grid/scratch/grossman/simulation_out/finalRuns" #"/links/shared/patrick_clonex/simulation_out/finalRuns" #"/home/grossman/simulation_out/finalRuns"

mkdir -p $SIMDIR
# mkdir ${SIMDIR}/nullModel

GENS=4500
INITCELLS=1000000
# UUU=0.0000001 #1e-07

# mkdir ${SIMDIR}/${GENS}gens_N10e9_n10e6ratenormal1e-8nullModel
# ./clonex -f ${SIMDIR}/${GENS}gens_N10e9_n10e6ratenormal1e-8nullModel -g ${GENS} -n ${INITCELLS} -N 1000000000 -c 1 -M 1000 -s 0 -R 50 -u 0.00000001
./perform_simulations.sh ${SIMDIR}/${GENS}gens_N10e9_n10e6ratenormal1e-8 ${GENS} ${INITCELLS} 1000000000 0.00000001 # done before

# ./perform_simulations.sh ${SIMDIR}/${GENS}gens_N10e7_n10e6 ${GENS} ${INITCELLS} 10000000
# ./perform_simulations.sh ${SIMDIR}/${GENS}gens_N10e8_n10e6 ${GENS} ${INITCELLS} 100000000
#  ./perform_simulations.sh ${SIMDIR}/${GENS}gens_N10e9_n10e6 ${GENS} ${INITCELLS} 1000000000 # done before
# ./perform_simulations.sh ${SIMDIR}/${GENS}gens_N10e10_n10e6 ${GENS} ${INITCELLS} 10000000000

# ./perform_simulations.sh ${SIMDIR}/4400gens_N10e7_n10e0 4400 1 10000000
# ./perform_simulations.sh ${SIMDIR}/4400gens_N10e8_n10e0 4400 1 100000000
# ./perform_simulations.sh ${SIMDIR}/4400gens_N10e9_n10e0 4400 1 1000000000
# ./perform_simulations.sh ${SIMDIR}/400gens_N10e10_n10e0 4400 1 10000000000

# ./perform_simulations.sh ${SIMDIR}/test_1200gens_N10e9 1200 1 1000000000
# ./perform_simulations.sh ${SIMDIR}/test_1200gens_N10e9_n10e6 1200 1000000 1000000000
# ./perform_simulations.sh ${SIMDIR}/test_1800gens_N10e9_n10e6 1200 1000000 10000000000
# 
# ./perform_simulations.sh ${SIMDIR}/test_1800gens_N10e9 1800 1 1000000000
# ./perform_simulations.sh ${SIMDIR}/test_1500gens_N10e9 1500 1 1000000000
# 
# ./perform_simulations.sh ${SIMDIR}/test_1500gens_N10e8 1500 1 100000000
# ./perform_simulations.sh ${SIMDIR}/test_1800gens_N10e8 1800 1 100000000

## vary mutation rate

# mkdir ${SIMDIR}/${GENS}gens_N10e9_n10e6normalrate1e-7nullModel
# 
# ./clonex -f ${SIMDIR}/${GENS}gens_N10e9_n10e6ratenormal1e-7nullModel -g ${GENS} -n ${INITCELLS} -N 1000000000 -c 1 -M 1000 -s 0 -R 50 -u 0.0000001
# ./perform_simulations.sh ${SIMDIR}/${GENS}gens_N10e9_n10e6ratenormal1e-7 ${GENS} ${INITCELLS} 1000000000 0.0000001
