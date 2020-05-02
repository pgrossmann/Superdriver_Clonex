#!/bin/bash

SIM_DIR=$1
GENERATIONS=$2
NINIT=$3
NFINAL=$4
UU=$5 #0.0000001 #1e-07

mkdir $SIM_DIR

### function definition ##
perform(){
#   MM=$1
#   ll=$2
#   dd=$3
  
  local WHICH_CASE=$1
  local MM=$2[@]
#   local tt=$3[@]
  local dd=$3[@]
  local cc=$4[@]
  #echo ${!ref}
  for M in ${!MM}; do
#     for t in ${!tt}; do
      for d in ${!dd}; do
	for c in ${!cc}; do
# 	  TMP_OUTDIR=$SIM_DIR/${WHICH_CASE}_mutator${M}superdriver${t}driver${d}doubling${c}
	  TMP_OUTDIR=$SIM_DIR/${WHICH_CASE}_mutator${M}driver${d}factor${c}ratenormal${UU}
	  if [[ -d "$TMP_OUTDIR" ]]
	    then 
              echo "directory ${TMP_OUTDIR} exists, clonex is not run"
            else
#               echo "test ${TMP_OUTDIR}"
	      mkdir $TMP_OUTDIR && ./clonex -s $d -M $M -c $c -R 50 -g $GENERATIONS -n $NINIT -N $NFINAL -f $TMP_OUTDIR -u $UU &
# 	      wait
          fi
	done
      done
#       wait
#     done# 7 folders
  done
#  wait
}

###########
### runs ##
###########

##################################################################################
## mutator phenotype = False, superdriver = False, driver = True, passenger = True
## case A

# MM=1000
# ll=0
# 
# for d in $(seq 0.005 0.01 0.1); do 
#   TMP_OUTDIR=mutator${MM}superdriver${ll}driver$d
#   mkdir $TMP_OUTDIR
#   ./clonex -s $i -l $ll -M $MM -R 50 -f $TMP_OUTDIR
# done

mutators=( 1000 )
# superdrivers=( 0 )
# drivers=( 0.006 0.008 0.01 0.012 0.14 )
# drivers=( 0.001 0.005 0.01 0.05 0.1 0.5 1 )
drivers=( 0.001 0.005 0.01 0.05 0.1 )
multiplier=( 1 )

# perform A mutators drivers multiplier # 5 folders

#################################################################################
## mutator phenotype = False, superdriver = True, driver = True, passenger = True
## case B

mutators=( 1000 )
# superdrivers=( 1 )
# drivers=( 0.006 0.008 0.01 0.012 0.14 0.18 0.2 )
# drivers=( 0.001 0.005 0.01 0.05 0.1 0.5 1 )
drivers=( 0.005 0.01 0.02 0.025 0.03 0.04 0.05 )
# multiplier=( 1.5 1.75 2 2.25 2.5 2.75 3 )
multiplier=( 2.1 2.2 2.3 2.4 2.5 )

perform B mutators drivers multiplier # 24 folders

#################################################################################
## mutator phenotype = True, superdriver = False, driver = True, passenger = True
## case D

mutators=( 5 4 3 2 1 )
# superdrivers=( 0 )
# drivers=(  0.006 0.008 0.01 0.012 0.14  )
# drivers=( 1 0.5 0.1 0.05 0.01 0.005 0.001 )
drivers=( 0.1 0.05 0.01 0.005 0.001 )
multiplier=( 2 )

# perform D mutators drivers multiplier  # 25 folders

#################################################################################
## mutator phenotype = True, superdriver = True, driver = True, passenger = True
## case C

mutators=( 5 4 3 2 1 )
# superdrivers=( 1 )
# drivers=( 0.001 0.005 0.01 0.05 0.1 0.5 1 )
# drivers=( 1 0.5 0.1 0.05 0.01 0.005 0.001 )
drivers=( 0.1 0.05 0.01 0.005 0.001 ) # 5
multiplier=( 1.5 2 2.5 3 ) # 4
# multiplier=( 1.5 1.75 2 2.25 2.5 2.75 3 )

# perform C mutators drivers multiplier  # 100 folders
