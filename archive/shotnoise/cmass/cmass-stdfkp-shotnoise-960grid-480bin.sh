#!/bin/bash
# P(k) with standard FKP shotnoise and normalization  
read -p "Data Release # and Version # (eg. dr12v1)" -e drv  # Read data release 
read -p "North or South (N or S)" -e ns

file_prefix="cmass-"$drv"-"
file_ext=".dat"; file_random_ext=".ran.dat"

if [[ $drv = "dr12v4" ]]; then 
    file_flag="-Reid-weights-zlim"
else 
    file_flag="-Anderson-weights-zlim"
fi 
output_flag="-stdfkp-shotnoise"

sscale=3600.0; Rbox=1800.0; box="3600"; P0=20000

fft_ext="-ngalsys-"$box"lbox-960grid.dat"
fft_ranext="-ngalsys-"$box"lbox-960grid.ran.dat"
power_ext="-ngalsys-"$box"lbox-960grid-480bin.dat"
FFT="FFT-"
power="power-"

ifort -O3 -o exe/FFT-stdfkp-shotnoise-960grid.exe FFT-stdfkp-shotnoise-960grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/power-stdfkp-shotnoise-960grid-480bin.exe power-stdfkp-shotnoise-960grid-480bin.f

datadir="/mount/riachuelo1/hahn/data/"
FFTdir="/mount/riachuelo1/hahn/FFT/"
powerdir="/mount/riachuelo1/hahn/power/"

cmass_file=$datadir$file_prefix$ns$file_flag$file_ext
random_file=$datadir$file_prefix$ns$file_flag$file_random_ext
ls -ltr $cmass_file
ls -ltr $random_file
#if [[ ! -a $cmass_file ]]; then 
#    # First convert .fits data and random data in local directory to ASCII file 
#    if [ $ns = "N" ]; then 
#        idl -e "build_cmass_ns_data,'"$drv"',/north"
#    fi 
#    if [ $ns = "S" ]; then 
#        idl -e "build_cmass_ns_data,'"$drv"',/south"
#    fi
#fi

# Compute FFT for CMASS data and random data
FFT_file=$FFTdir$FFT$file_prefix$ns$file_flag$output_flag$fft_ext
exe/FFT-stdfkp-shotnoise-960grid.exe $Rbox 0 $P0 $cmass_file $FFT_file
echo $FFT_file

FFT_rand_file=$FFTdir$FFT$file_prefix$ns$file_flag$fft_ranext
if [[ ! -a $FFT_rand_file ]]; then 
    exe/FFT-stdfkp-shotnoise-960grid.exe $Rbox 1 $P0 $random_file $FFT_rand_file
    echo $FFT_rand_file
fi

# Compute power spectrum for CMASS
power_file=$powerdir$power$file_prefix$ns$file_flag$output_flag$power_ext
exe/power-stdfkp-shotnoise-960grid-480bin.exe $FFT_rand_file $FFT_file $power_file $sscale
echo $power_file
