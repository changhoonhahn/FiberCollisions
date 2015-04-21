#!/bin/bash
# P(k) for the tiling mock catalog with NO weighting scheme
mock_file="cmass-boss5003sector-icoll012.dat"
rand_file="randoms-boss5003-icoll012-vetoed.dat"
nbar_file="nbar-cmass-boss5003sector-icoll012.dat"

P0=20000; sscale=4000.0; Rbox=2000.0; box="4000"; grid="360"

FFT="FFT_"; power="power_"

datadir="/mount/riachuelo1/hahn/data/tiling_mocks/"
randdir="/mount/riachuelo1/hahn/data/tiling_mocks/"
FFTdir="/mount/riachuelo1/hahn/FFT/tiling_mocks/"
powerdir="/mount/riachuelo1/hahn/power/tiling_mocks/"

ifort -O3 -o FFT-fkp-tiling-mock-360grid.exe FFT-fkp-tiling-mock-360grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o power-fkp-mock-360grid-180bin.exe power-fkp-mock-360grid-180bin.f

FFTrandname=$FFTdir$FFT$rand_file".grid"$grid".P0"$P0".box"$box
./FFT-fkp-tiling-mock-360grid.exe $Rbox 1 $P0 $datadir$nbar_file $randdir$rand_file $FFTrandname
echo $FFTrandname

FFTname=$FFTdir$FFT$mock_file".grid"$grid".P0"$P0".box"$box
powername=$powerdir$power$mock_file".grid"$grid".P0"$P0".box"$box

./FFT-fkp-tiling-mock-360grid.exe $Rbox 0 $P0 $datadir$nbar_file $datadir$mock_file $FFTname
echo $FFTname
./power-fkp-mock-360grid-180bin.exe $FFTname $FFTrandname $powername $sscale
echo $powername
