#!/bin/bash 
source activate acpc
export OMP_NUM_THREADS=1
python Tracking.py $1 $2 $3 $4


