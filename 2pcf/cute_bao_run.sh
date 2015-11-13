#!/bin/bash

cd /home/users/kww231/src/CUTE-1.2/CUTE
echo "n_mocks="; read n

mock="Nseries"

for i in $(seq 1 $n); do 
    for correction in "true"; do 
        param_file="/mount/riachuelo1/hahn/2pcf/param_files/param_"$correction"_"$mock$i"_bao.txt"
        ./CUTE $param_file 
    done
done
