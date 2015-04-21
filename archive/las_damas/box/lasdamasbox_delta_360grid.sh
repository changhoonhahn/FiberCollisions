#!/bin/bash
data_prefix="mock_gamma_lrg21p2_zmo_Oriana_1"
data_suffix="_z0p342_fof_b0p2.zdist.ff"
mock_suffix="_z0p342_fof_b0p2.zdist.fibcoll.dat"
wght_suffix="_z0p342_fof_b0p2.zdist.fibcoll.delta.dat"
parm_suffix="_z0p342_fof_b0p2.zdist.param.dat"
rand_file="sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat"

sscale=3600.0; Rbox=1800.0; box="3600"; FFT="FFT_"; power="power_" 
grid="360"; nbin="180"; Lx=360; Ly=360; Lz=360

datadir="/mount/chichipio2/rs123/MOCKS/LRG21p2_zm_box/gaussian/zspace/"
mockdir="/mount/riachuelo1/hahn/data/las_damas/box/"
FFTdir="/mount/riachuelo1/hahn/FFT/las_damas/box/"
powerdir="/mount/riachuelo1/hahn/power/las_damas/box/"

echo "n="; read n 
#ifort -O4 -o export_lasdamas_ff.exe export_lasdamas_ff.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw
ifort -O4 -o FFT_lasdamasbox_delta_aniso.exe FFT_lasdamasbox_delta_aniso.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw
ifort -O4 -o power_lasdamasbox_aniso.exe power_lasdamasbox_aniso.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw

for i in $(seq -f "%03g" 1 $n); do 
    echo $i
    data_file=$datadir$data_prefix$i$data_suffix
    mock_file=$mockdir$data_prefix$i$mock_suffix
    parm_file=$mockdir$data_prefix$i$parm_suffix
    FFT_file=$FFTdir$FFT$data_prefix$i$wght_suffix".grid"$grid".P0"$P0".box"$box
    power_file=$powerdir$power$data_prefix$i$wght_suffix".grid"$grid".P0"$P0".box"$box
    #./export_lasdamas_ff.exe $data_file $mock_file $parm_file
    ./FFT_lasdamasbox_delta_aniso.exe $Lx $Ly $Lz $parm_file $mock_file $FFT_file
    ./power_lasdamasbox_aniso.exe $FFT_file $power_file $nbin
done
