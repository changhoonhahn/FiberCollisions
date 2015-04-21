#!/bin/bash
name0="cmass_dr11_north_ir4"; nameend=".v11.0.wghtv.txt"
nbarfname="nbar-cmass-dr11may22-N-Anderson.dat"
sscale=3600.0; Rbox=1800.0; box="3600"
FFT="FFT_"; power="power_"
grid="360"; P0=20000

echo "n="; read n 

ifort -O3 -o FFT-fkp-pthalo-v11p0-upweight-360grid.exe FFT-fkp-pthalo-v11p0-upweight-360grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o power-fkp-pthalo-dr11-v11p0.exe power-fkp-pthalo-dr11-v11p0.f
datadir="/mount/riachuelo1/hahn/data/manera_mock/v11p0/"
FFTdir="/mount/riachuelo1/hahn/FFT/manera_mock/v11p0/"
powerdir="/mount/riachuelo1/hahn/power/manera_mock/v11p0/"

randname="cmass_dr11_north_600_randoms_ir4_combined.v11.0.txt"
FFTrandname=$FFTdir$FFT"cmass_dr11_north_600_randoms_ir4_combined.v11.0.upweight.grid"$grid".P0"$P0".box"$box

if [[ ! -a $datadir$randname ]]; then 
    echo $datadir$randname 
    ipython -noconfirm_exit rand_combine-dr11-v11p0-north.py $n
fi
if [[ ! -a $FFTrandname ]]; then 
    echo $FFTrandname 
    ./FFT-fkp-pthalo-v11p0-upweight-360grid.exe $Rbox 1 $P0 $datadir$nbarfname $datadir$randname $FFTrandname
fi

for i in $(seq -f "%03g" 1 $n); do 
    echo $i
    echo $fname
    fname=$datadir$name0$i$nameend
    FFTname=$FFTdir$FFT$name0$i".v11.0.upweight.grid"$grid".P0"$P0".box"$box
    powername=$powerdir$power$name0$i".v11.0.600randoms.upweight.grid"$grid".P0"$P0".box"$box
    echo $FFTname 
    if [[ ! -a $FFTname ]]; then 
        ./FFT-fkp-pthalo-v11p0-upweight-360grid.exe $Rbox 0 $P0 $datadir$nbarfname $fname $FFTname
    fi
    echo $powername
    ./power-fkp-pthalo-dr11-v11p0.exe $FFTname $FFTrandname $powername $sscale
done 
