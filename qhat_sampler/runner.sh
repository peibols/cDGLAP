#!/bin/bash

N=1000
cent=("0-5" "5-10" "10-20" "20-30" "30-40" "40-50" "50-60" "60-70")

for ii in `seq 1 1`;
do

  for jj in `seq 0 7`;
  do

    echo "Doing cent ${cent[$jj]}"
    ./main_PbPb $N $ii ${cent[$jj]}

  done

done
