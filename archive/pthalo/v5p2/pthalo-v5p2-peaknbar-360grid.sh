#!/bin/bash
name0="cmass_dr10_north_ir4"; nameend=".v5.2.wghtv.txt"
nameend_peak=".v5.2.peakcorr.txt"
nbarfname="nbar-dr10v5-N-Anderson.dat"
sscale=3600.0; Rbox=1800.0; box="3600"
FFT="FFT_"; power="power_"
grid="360"; P0=20000
sigma="5.44"; fpeak="1.0"

echo "n="; read n 
ifort -O3 -o mock-cp-dlos-peak-pthalo-v11p0.exe mock-cp-dlos-peak-pthalo-v11p0.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o FFT-fkp-pthalo-v11p0-peaknbar-corr-360grid.exe FFT-fkp-pthalo-v11p0-peaknbar-corr-360grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o power-fkp-pthalo-v11p0.exe power-fkp-pthalo-v11p0.f
datadir="/mount/chichipio2/hahn/data/manera_mock/v5p2/"
fibcolldir="/mount/riachuelo1/hahn/data/manera_mock/v11p0/fibcoll/"
FFTdir="/mount/riachuelo1/hahn/FFT/manera_mock/v11p0/fibcoll/"
powerdir="/mount/riachuelo1/hahn/power/manera_mock/v11p0/fibcoll/"

randname="cmass_dr11_north_600_randoms_ir4_combined.v11.0.txt"
FFTrandname="/mount/riachuelo1/hahn/FFT/manera_mock/v11p0/"$FFT"cmass_dr11_north_600_randoms_ir4_combined.v11.0.upweight.grid"$grid".P0"$P0".box"$box

if [[ ! -a $datadir$randname ]]; then 
    echo $datadir$randname 
    ipython -noconfirm_exit rand_combine-dr11-v11p0-north.py 600 
fi
if [[ ! -a $FFTrandname ]]; then 
    echo $FFTrandname 
    ./FFT-fkp-pthalo-v11p0-peaknbar-corr-360grid.exe $Rbox 1 $P0 $datadir$nbarfname $datadir$randname $FFTrandname
fi

for i in $(seq -f "%03g" 1 $n); do 
    echo $i 
    fname=$datadir$name0$i$nameend
    corrfile=$fibcolldir$name0$i$nameend_peak
    ./mock-cp-dlos-peak-pthalo-v11p0.exe $datadir$nbarfname $fname $corrfile &
    jobcnt=(`jobs -p`)
    if [ ${#jobcnt[@]} -gt 4 ]; then 
        wait
    fi
done 
wait

for i in $(seq -f "%03g" 1 $n); do 
    echo $i
    corrfile=$fibcolldir$name0$i$nameend_peak
    FFTname=$FFTdir$FFT$name0$i".v11.0.peaknbar.sigma"$sigma"fpeak"$fpeak".grid"$grid".P0"$P0".box"$box
    powername=$powerdir$power$name0$i".v11.0.600randoms.peaknbar.sigma"$sigma"fpeak"$fpeak".grid"$grid".P0"$P0".box"$box
    echo $fname
    ./FFT-fkp-pthalo-v11p0-peaknbar-corr-360grid.exe $Rbox 0 $P0 $datadir$nbarfname $corrfile $FFTname
    echo $FFTname 
    ./power-fkp-pthalo-v11p0.exe $FFTname $FFTrandname $powername $sscale
    echo $powername
done 
