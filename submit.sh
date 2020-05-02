#!/bin/bash
#$ -q mpi04-ht.q
#$ -cwd
#$ -N addSuperdriver
#$ -V
# #$ -pe 25

#/home/grossman/Dropbox/current_projects/master-thesis_cancer-progression/code/simona_new/C/clonex -h

/home/grossman/paper/code/C/perform_several_simulations.sh
