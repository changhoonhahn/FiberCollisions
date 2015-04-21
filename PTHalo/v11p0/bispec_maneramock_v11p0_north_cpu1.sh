#!/bin/bash
prefix="FFT_cmass_dr11_north_ir4"
bispecpre="BISP_cmass_dr11_north_ir4"
suffix=".v11.0.upweight"
sscale=3600.0
P0=20000
box="3600"
grid="360"
nmax="40"
nstep="3"
ncut="3"
iflag=2
FFTdir="/mount/riachuelo1/hahn/FFT/manera_mock/v11p0/"
BISdir="/mount/riachuelo1/hahn/bispec/manera_mock/v11p0/"

ifort -O4 -o bisp_maneramock_v5p2.exe bisp_maneramock_v5p2.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw

randFFTname=$FFTdir"FFT_cmass_dr11_north_600_randoms_ir4_combined.v11.0.grid"$grid".P0"$P0".box"$box
for i in $(seq -f "%03g" 1 100); do 
    echo $i 
    FFTname=$FFTdir$prefix$i$suffix".grid"$grid".P0"$P0".box"$box
    bispec=$BISdir$bispecpre$i$suffix".grid"$grid".nmax"$nmax".nstep"$nstep".P0"$P0".box"$box
    echo $FFTname
    echo $randFFTname 
    echo $bispec

    ./bisp_maneramock_v5p2.exe $iflag $randFFTname $FFTname $bispec
done
