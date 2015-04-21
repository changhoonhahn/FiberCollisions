#!/bin/bash
mock_prefix="sdssmock_gamma_lrgFull_zm_oriana"
mock_suffix="_no.rdcz.dat"
rand_file="sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat"

P0=20000; sscale=3600.0; Rbox=1800.0; box="3600"; FFT="FFT_"; power="power_"; grid="960"

datadir="/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/"
randdir="/mount/chichipio2/rs123/MOCKS/randoms/"
FFTdir="/mount/riachuelo1/hahn/FFT/las_damas/fiber_collision/"
powerdir="/mount/riachuelo1/hahn/power/las_damas/fiber_collision/"

echo "n="; read n 
ifort -O3 -o FFT-fkp-lasdamas-mock-960grid.exe FFT-fkp-lasdamas-mock-960grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o power-fkp-mock-960grid-480bin.exe power-fkp-mock-960grid-480bin.f

FFTrandname=$FFTdir$FFT$rand_file".grid"$grid".P0"$P0".box"$box
if [[ ! -a $FFTrandname ]]; then
    ./FFT-fkp-lasdamas-mock-960grid.exe $Rbox 1 $P0  $randdir$rand_file $FFTrandname
    echo $FFTrandname
fi

for i in $(seq 1 $n); do 
    for letter in "a" "b" "c" "d"; do 
        i=$( printf '%02d' $i ) 
        echo $i
        echo $letter
        mock_file=$datadir$mock_prefix$i$letter$mock_suffix
        FFTname=$FFTdir$FFT$mock_prefix$i$letter$mock_suffix".grid"$grid".P0"$P0".box"$box
        powername=$powerdir$power$mock_prefix$i$letter$mock_suffix".grid"$grid".P0"$P0".box"$box
         
        ./FFT-fkp-lasdamas-mock-960grid.exe $Rbox 0 $P0 $mock_file $FFTname
        echo $FFTname
        ./power-fkp-mock-960grid-480bin.exe $FFTname $FFTrandname $powername $sscale
        echo $powername
    done 
done
