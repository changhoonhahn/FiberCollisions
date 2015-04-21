#!/bin/bash
mock_prefix="sdssmock_gamma_lrgFull_zm_oriana"
mock_suffix="_no.rdcz.fibcoll.dat"
rand_file="sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat"

P0=20000; sscale=3600.0; Rbox=1800.0; box="3600"; FFT="FFT_"; power="power_"; grid="360"

datadir="/mount/riachuelo1/hahn/data/las_damas/fibcoll/"
randdir="/mount/chichipio2/rs123/MOCKS/randoms/"
FFTdir="/mount/riachuelo1/hahn/FFT/las_damas/fiber_collision/"
powerdir="/mount/riachuelo1/hahn/power/las_damas/fiber_collision/"

echo "n="; read n 
ifort -O3 -o FFT-fkp-delta-corr-lasdamas-mock-360grid.exe FFT-fkp-delta-corr-lasdamas-mock-360grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o power-fkp-mock-360grid-180bin.exe power-fkp-mock-360grid-180bin.f

for i in $(seq -f "%02g" 1 $n); do
    for letter in "a" "b" "c" "d"; do
        mock_file=$datadir$mock_prefix$i$letter$mock_suffix
        echo $mock_file
        if [[ ! -a $mock_file ]]; then 
            idl -e "build_mock_wcp_assign,"$i",'"$letter"'"
        fi 
        echo $i
        echo $letter
        jobcnt=(`jobs -p`)
        if [ ${#jobcnt[@]} -gt 3 ]; then
            wait
        fi 
    done;
    wait 
done; 
wait

FFTrandname=$FFTdir$FFT$rand_file".delta.grid"$grid".P0"$P0".box"$box
if [[ ! -a $FFTrandname ]]; then
    ./FFT-fkp-delta-corr-lasdamas-mock-360grid.exe $Rbox 1 $P0  $randdir$rand_file $FFTrandname
    echo $FFTrandname
fi

for i in $(seq -f "%02g" 1 $n); do 
    for letter in "a" "b" "c" "d"; do 
        echo $i
        echo $letter
        mock_file=$datadir$mock_prefix$i$letter$mock_suffix
        FFTname=$FFTdir$FFT$mock_prefix$i$letter$mock_suffix".delta.grid"$grid".P0"$P0".box"$box
        powername=$powerdir$power$mock_prefix$i$letter$mock_suffix".delta.grid"$grid".P0"$P0".box"$box
         
        ./FFT-fkp-delta-corr-lasdamas-mock-360grid.exe $Rbox 0 $P0 $mock_file $FFTname
        echo $FFTname
        ./power-fkp-mock-360grid-180bin.exe $FFTname $FFTrandname $powername $sscale
        echo $powername
    done 
done
