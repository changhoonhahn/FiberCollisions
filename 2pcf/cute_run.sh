#!/bin/bash

cd /home/users/kww231/src/CUTE-1.2/CUTE
echo "n_mocks="; read n

mock="Nseries"
scale="large"

for i in $(seq 21 $n); do 
    for correction in "true" "upweighted" "collrm"; do 
        param_file="/mount/riachuelo1/hahn/2pcf/param_files/param_"$correction"_"$mock$i"_"$scale".txt"
        ./CUTE $param_file 
    done
done
