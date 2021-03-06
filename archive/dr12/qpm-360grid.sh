#!/bin/bash
mock_prefix="sdssmock_gamma_lrgFull_zm_oriana"
mock_suffix="_no.rdcz.dat"; fibcoll_suffix="_no.rdcz.fibcoll.dat"
rand_file="sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat"

P0=20000; sscale=3600.0; Rbox=1800.0; box="3600"; FFT="FFT_"; power="power_"; grid="360"

datadir="/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/"     # data directory
fibcoll_datadir="/mount/riachuelo1/hahn/data/LasDamas/Geo/"    # fibercollided data directory  
randdir="/mount/chichipio2/rs123/MOCKS/randoms/"                    # random data directory
FFTdir="/mount/riachuelo1/hahn/FFT/LasDamas/Geo/"                   # FFT directory
powerdir="/mount/riachuelo1/hahn/power/LasDamas/Geo/"               # powerspectrum directory

read -p "Number of LasDamasGeo mocks: " -e n  
read -p "Correction method (true, delta, or peak): " -e correction  

ifort -O3 -o exe/peakcorr_ldg_fibcollmock.exe peakcorr_ldg_fibcollmock.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/FFT_ldg_fkp_now_360grid.exe FFT_ldg_fkp_now_360grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/FFT_ldg_fkp_w_360grid.exe FFT_ldg_fkp_w_360grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/power_ldg_fkp_360grid_180bin.exe power_ldg_fkp_360grid_180bin.f

FFTrandname=$FFTdir$FFT$rand_file".grid"$grid".P0"$P0".box"$box  # same FFT for all correction method   
if [[ ! -a $FFTrandname ]]; then
    ./exe/FFT_ldg_fkp_now_360grid.exe $Rbox 1 $P0  $randdir$rand_file $FFTrandname
    echo $FFTrandname
fi

if [ $correction = "true" ]; then 
    # Compute the true P(k) for LasDamasBox using no weight mocks 
    echo "True Correction Method" 
    for i in $(seq -f "%02g" 1 $n); do
        for letter in "a" "b" "c" "d"; do
            echo $i$letter
            mock_file=$datadir$mock_prefix$i$letter$mock_suffix
            FFTname=$FFTdir$FFT$mock_prefix$i$letter$mock_suffix".grid"$grid".P0"$P0".box"$box
            ./exe/FFT_ldg_fkp_now_360grid.exe $Rbox 0 $P0 $mock_file $FFTname & 
            echo $FFTname
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 3 ]; then
                wait
            fi 
        done
    done
    for i in $(seq -f "%02g" 1 $n); do
        for letter in "a" "b" "c" "d"; do
            echo $i$letter
            FFTname=$FFTdir$FFT$mock_prefix$i$letter$mock_suffix".grid"$grid".P0"$P0".box"$box
            powername=$powerdir$power$mock_prefix$i$letter$mock_suffix".grid"$grid".P0"$P0".box"$box
            ./exe/power_ldg_fkp_360grid_180bin.exe $FFTname $FFTrandname $powername $sscale
            echo $powername
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 3 ]; then
                wait
            fi 
        done
    done
elif [ $correction = "delta" ] || [ $correction = "peak" ]; then 
    echo $correction" Correction Method" 
    for i in $(seq -f "%02g" 1 $n); do
        for letter in "a" "b" "c" "d"; do
            fibcoll_file=$fibcoll_datadir$mock_prefix$i$letter$fibcoll_suffix
            echo $i$letter
            if [[ ! -a $fibcoll_file ]]; then       # Generate fibercollided mocks, if not already there
                idl -e "ldg_fibcollmock_wcp_assign,"$i",'"$letter"'" &
            fi 
            echo $fibcoll_file
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 3 ]; then
                wait
            fi 
        done
    done
    if [ $correction = "delta" ]; then 
        echo "Delta Correction Method" 
        for i in $(seq -f "%02g" 1 $n); do 
            for letter in "a" "b" "c" "d"; do 
                echo $i$letter
                fibcoll_file=$fibcoll_datadir$mock_prefix$i$letter$fibcoll_suffix
                FFTname=$FFTdir$FFT$mock_prefix$i$letter$fibcoll_suffix".delta.grid"$grid".P0"$P0".box"$box
                ./exe/FFT_ldg_fkp_w_360grid.exe $Rbox 0 $P0 $fibcoll_file $FFTname
                echo $FFTname
                jobcnt=(`jobs -p`)
                if [ ${#jobcnt[@]} -gt 3 ]; then
                    wait
                fi 
            done 
        done
        for i in $(seq -f "%02g" 1 $n); do 
            for letter in "a" "b" "c" "d"; do 
                echo $i$letter
                FFTname=$FFTdir$FFT$mock_prefix$i$letter$fibcoll_suffix".delta.grid"$grid".P0"$P0".box"$box
                powername=$powerdir$power$mock_prefix$i$letter$fibcoll_suffix".delta.grid"$grid".P0"$P0".box"$box
                ./exe/power_ldg_fkp_360grid_180bin.exe $FFTname $FFTrandname $powername $sscale
                echo $powername
                jobcnt=(`jobs -p`)
                if [ ${#jobcnt[@]} -gt 3 ]; then
                    wait
                fi 
            done 
        done
    else 
        echo "Peak Correction Method" 
        sigma=5.3       # hardcoded peak sigma value
        sigma_suffix="_no.rdcz.fibcoll.dat.peak.sigma"$sigma".fpeak1.0"
        nbar_file="/mount/riachuelo1/hahn/data/manera_mock/dr11/nbar-dr10v5-N-Anderson.dat"
        for i in $(seq -f "%02g" 1 $n); do
            for letter in "a" "b" "c" "d"; do
                fibcoll_file=$fibcoll_datadir$mock_prefix$i$letter$fibcoll_suffix
                peakcorr_file=$fibcoll_datadir$mock_prefix$i$letter$sigma_suffix
                echo $i$letter
                if [[ ! -a $peakcorr_file ]]; then
                    ./exe/peakcorr_ldg_fibcollmock.exe $sigma $nbar_file $fibcoll_file $peakcorr_file &
                fi
                jobcnt=(`jobs -p`)
                if [ ${#jobcnt[@]} -gt 3 ]; then
                    wait
                fi
            done
        done
        for i in $(seq -f "%02g" 1 $n); do
            for letter in "a" "b" "c" "d"; do
                echo $i$letter
                peakcorr_file=$fibcoll_datadir$mock_prefix$i$letter$sigma_suffix
                FFTname=$FFTdir$FFT$mock_prefix$i$letter$sigma_suffix".grid"$grid".P0"$P0".box"$box
                ./exe/FFT_ldg_fkp_w_360grid.exe $Rbox 0 $P0 $peakcorr_file $FFTname &
                echo $FFTname
                jobcnt=(`jobs -p`)
                if [ ${#jobcnt[@]} -gt 3 ]; then
                    wait
                fi
            done
        done
        for i in $(seq -f "%02g" 1 $n); do
            for letter in "a" "b" "c" "d"; do
                echo $i$letter
                FFTname=$FFTdir$FFT$mock_prefix$i$letter$sigma_suffix".grid"$grid".P0"$P0".box"$box
                powername=$powerdir$power$mock_prefix$i$letter$sigma_suffix".grid"$grid".P0"$P0".box"$box
                ./exe/power_ldg_fkp_360grid_180bin.exe $FFTname $FFTrandname $powername $sscale &
                echo $powername
                jobcnt=(`jobs -p`)
                if [ ${#jobcnt[@]} -gt 3 ]; then
                    wait
                fi
            done
        done
    fi
else 
    echo "Enter correct correction method: true, delta, peak. Entry is case sensititve!"
fi 

