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
read -p "Correction method (true, upweight, or peak): " -e correction  

ifort -O3 -o exe/peakcorr_ldg_fibcollmock.exe peakcorr_ldg_fibcollmock.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/FFT_ldg_fkp_now_360grid.exe FFT_ldg_fkp_now_360grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/FFT_ldg_fkp_w_360grid.exe FFT_ldg_fkp_w_360grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/power_ldg_fkp_360grid_180bin.exe power_ldg_fkp_360grid_180bin.f

if [ $correction = "true" ]; then 
    # Compute true random FFT (only one random file)  
    FFTrandname=$FFTdir$FFT$rand_file".grid"$grid".P0"$P0".box"$box  # same FFT for all correction method   
    if [[ ! -a $FFTrandname ]]; then
        ./exe/FFT_ldg_fkp_now_360grid.exe $Rbox 1 $P0  $randdir$rand_file $FFTrandname
        echo $FFTrandname
    fi

    # Compute the true P(k) for LasDamasBox using no weight mocks 
    echo "True Correction Method" 
    # Compute FFT for true mock file 
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
    # Compute P(k) for true mock file  
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
# Fiber collision corrections 
elif [ $correction = "upweight" ] || [ $correction = "peak" ]; then 
    echo $correction" Correction Method" 
    
    # Apply fiber collisions on the true mock file 
    # detects fiber collisions on 62'' angular scales 
    # then imposes upweighting scheme 
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

    # for upweight correction method 
    if [ $correction = "upweight" ]; then 
        echo "Delta Correction Method" 
        # Compute upweight random FFT: 
        rand_corr_flag=".1a1b1c1d2a2b2c2d.upweight.fibcoll"             # upweight random catalog currently hard-coded
        rand_corr_file=$randdir$rand_file$rand_corr_flag
        FFTrandname=$FFTdir$FFT$rand_corr_file".grid"$grid".P0"$P0".box"$box 
        if [[ ! -a $FFTrandname ]]; then
            ./exe/FFT_ldg_fkp_now_360grid.exe $Rbox 1 $P0  $rand_corr_file $FFTrandname
            echo $FFTrandname" WRITTEN"
        fi
        
        # Compute upweight mock FFT
        for i in $(seq -f "%02g" 1 $n); do 
            for letter in "a" "b" "c" "d"; do 
                echo $i$letter
                fibcoll_file=$fibcoll_datadir$mock_prefix$i$letter$fibcoll_suffix
                FFTname=$FFTdir$FFT$mock_prefix$i$letter$fibcoll_suffix".upweight.grid"$grid".P0"$P0".box"$box
                ./exe/FFT_ldg_fkp_w_360grid.exe $Rbox 0 $P0 $fibcoll_file $FFTname
                echo $FFTname
                jobcnt=(`jobs -p`)
                if [ ${#jobcnt[@]} -gt 3 ]; then
                    wait
                fi 
            done 
        done

        # Compute upweight corrected P(k)
        for i in $(seq -f "%02g" 1 $n); do 
            for letter in "a" "b" "c" "d"; do 
                echo $i$letter
                FFTname=$FFTdir$FFT$mock_prefix$i$letter$fibcoll_suffix".upweight.grid"$grid".P0"$P0".box"$box
                powername=$powerdir$power$mock_prefix$i$letter$fibcoll_suffix".upweight.grid"$grid".P0"$P0".box"$box
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
        #read -p "i_peak (e.g. 0-10): " -e i_peak 
        for i_peak in $(seq 5 5); do 
            if [ $i_peak -ne 10 ]; then 
                peakfrac="0."$i_peak
            else 
                peakfrac="1.0"
            fi 
            echo $peakfrac
            sigma=5.3       # hardcoded best-fit peak sigma value
            corr_flag=".peak.sigma"$sigma".fpeak"$peakfrac
            sigma_suffix="_no.rdcz.fibcoll.dat"$corr_flag 

            nbar_file="/mount/riachuelo1/hahn/data/nbar-junk.dat"    # just an arbitrary nbar(z). only used for z values because lazy 
            # Correct upweight fibercollided mocks using the peak correction method 
            for i in $(seq -f "%02g" 1 $n); do
                for letter in "a" "b" "c" "d"; do
                    fibcoll_file=$fibcoll_datadir$mock_prefix$i$letter$fibcoll_suffix
                    peakcorr_file=$fibcoll_datadir$mock_prefix$i$letter$sigma_suffix
                    echo $i$letter
                    if [[ ! -a $peakcorr_file ]]; then
                        ./exe/peakcorr_ldg_fibcollmock.exe $sigma $peakfrac $nbar_file $fibcoll_file $peakcorr_file &
                        echo $peakcorr_file 
                    fi
                    jobcnt=(`jobs -p`)
                    if [ ${#jobcnt[@]} -gt 3 ]; then
                        wait
                    fi
                done
            done
            
            # assign corrected nbars to peak corrected mock catalogs  
            for i in $(seq -f "%02g" 1 $n); do
                for letter in "a" "b" "c" "d"; do
                    peakcorr_file=$fibcoll_datadir$mock_prefix$i$letter$sigma_suffix".corrnbar"
                    if [[ ! -a $peakcorr_file ]]; then
                        ipython /home/users/hahn/powercode/FiberCollision/build_append_corr_nbar.py "data" "lasdamasgeo" $i $letter "peaknbar" $sigma $peakfrac 
                        echo $peakcorr_file 
                    fi  

                    if [ ${#jobcnt[@]} -gt 3 ]; then
                        wait
                    fi
                done 
            done 

            rand_corr_flag=".allmocks"$corr_flag".corrnbar"             # peak corrected random catalog currently hard-coded
            rand_corr_file=$rand_file$rand_corr_flag
            # assign corrected nbars to random catalog 
            if [[ ! -a $fibcoll_datadir$rand_corr_file ]]; then 
                echo "WRITING    "$fibcoll_datadir$rand_corr_file
                ipython /home/users/hahn/powercode/FiberCollision/build_append_corr_nbar.py "random" "lasdamasgeo" $i $letter "peaknbar" $sigma $peakfrac
            else
                echo $fibcoll_datadir$rand_corr_file" ALREADY EXISTS!"
            fi 

            # Compute peak corrected random FFT: 
            FFTrandname=$FFTdir$FFT$rand_corr_file".grid"$grid".P0"$P0".box"$box 
            if [[ ! -a $FFTrandname ]]; then
                ./exe/FFT_ldg_fkp_w_360grid.exe $Rbox 1 $P0  $fibcoll_datadir$rand_corr_file $FFTrandname
                echo $FFTrandname" WRITTEN"
            else 
                echo $FFTrandname" ALREADY EXISTS!"
            fi
            
            # Compute peak corrected mock FFT
            for i in $(seq -f "%02g" 1 $n); do
                for letter in "a" "b" "c" "d"; do
                    echo $i$letter
                    peakcorr_file=$fibcoll_datadir$mock_prefix$i$letter$sigma_suffix".corrnbar"
                    FFTname=$FFTdir$FFT$mock_prefix$i$letter$sigma_suffix".corrnbar.grid"$grid".P0"$P0".box"$box
                    ./exe/FFT_ldg_fkp_w_360grid.exe $Rbox 0 $P0 $peakcorr_file $FFTname &
                    echo $FFTname
                    jobcnt=(`jobs -p`)
                    if [ ${#jobcnt[@]} -gt 3 ]; then
                        wait
                    fi
                done
            done
            # Compute peak corrected P(k) 
            for i in $(seq -f "%02g" 1 $n); do
                for letter in "a" "b" "c" "d"; do
                    echo $i$letter
                    FFTname=$FFTdir$FFT$mock_prefix$i$letter$sigma_suffix".corrnbar.grid"$grid".P0"$P0".box"$box
                    powername=$powerdir$power$mock_prefix$i$letter$sigma_suffix".corrnbar.grid"$grid".P0"$P0".box"$box
                    ./exe/power_ldg_fkp_360grid_180bin.exe $FFTname $FFTrandname $powername $sscale &
                    echo $powername
                    jobcnt=(`jobs -p`)
                    if [ ${#jobcnt[@]} -gt 3 ]; then
                        wait
                    fi
                done
            done 
        done
    fi
else 
    echo "Enter correct correction method: true, upweight, peak. Entry is case sensititve!"
fi 
