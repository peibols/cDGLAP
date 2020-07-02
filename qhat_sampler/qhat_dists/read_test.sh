#!/bin/bash

filename="cent_0-5.dat"

index=0
while IFS= read -r -a line;
do 
  index=$((index+1))
  arr=( $line )
  echo "$index qhat=${arr[0]} length=${arr[1]}"; 

done < $filename
