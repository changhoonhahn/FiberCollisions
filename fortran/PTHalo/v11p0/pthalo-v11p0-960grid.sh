#!/bin/bash
mock_prefix="cmass_dr11_north_ir4"; mock_suffix=".v11.0.wghtv.txt"
Nrandom=100         # Number of random catalog files to combine
rand_file="cmass_dr11_north_"$Nrandom"_randoms_ir4_combined.v11.0.txt"

P0=20000; sscale=3600.0; Rbox=1800.0; box="3600"; grid="960"
nmax="40"; nstep="3"; ncut="3"; iflag=2
FFT="FFT_"; power="power_"; bisp="bisp_"

datadir="/mount/riachuelo1/hahn/data/PTHalo/v11p0/"                 # data directory
fibcoll_datadir="/mount/riachuelo1/hahn/data/PTHalo/v11p0/"         # fibercollided data directory  
FFTdir="/mount/riachuelo1/hahn/FFT/PTHalo/v11p0/"                   # FFT directory
powerdir="/mount/riachuelo1/hahn/power/PTHalo/v11p0/"               # powerspectrum directory
bispdir="/mount/riachuelo1/hahn/bispec/PTHalo/v11p0/"               # bispectrum directory
    
nbar_file=$datadir"nbar-cmass-dr11may22-N-Anderson.dat"

read -p "PTHalo mocks: from " -e ni
read -p "PTHalo mocks: to " -e nf  
read -p "Correction method (noweight, wboss, upweight, or peak): " -e correction  
read -p "Powerspectrum or Bispectrum (pk or bk): " -e spectrum

ifort -O3 -o exe/fibcollmock_pth_peakcorr.exe fibcollmock_pth_peakcorr.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/FFT_pth_fkp_now_960grid.exe FFT_pth_fkp_now_960grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/FFT_pth_fkp_wboss_960grid.exe FFT_pth_fkp_wboss_960grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/FFT_pth_fkp_w_960grid.exe FFT_pth_fkp_w_960grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/power_pth_fkp_960grid_480bin.exe power_pth_fkp_960grid_480bin.f
ifort -O3 -o exe/bisp_pthalo_960grid.exe bisp_pthalo_960grid.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw

if [[ ! -a $datadir$rand_file  ]]; then                                     # If combined random file doesn't exist, build it 
   ipython rand_combine-pthalo-v11p0-north.py $Nrandom 
   echo $datadir$rand_file 
fi 
# FFT for random is the same for all fiber collision correction method 
FFTrandname=$FFTdir$FFT$rand_file".weight.grid"$grid".P0"$P0".box"$box 
if [[ ! -a $FFTrandname ]]; then
    echo "Random FFT"
    ./exe/FFT_pth_fkp_now_960grid.exe $Rbox 1 $P0 $nbar_file $datadir$rand_file $FFTrandname
    echo $FFTrandname
fi

if [ $correction = "noweight" ]; then 
    # Compute the true P(k) for PTHalo using only veto mask and no weights 
    echo "True Correction Method" 
    # Compute FFT 
    for i in $(seq -f "%03g" $ni $nf); do
        echo $i
        mock_file=$datadir$mock_prefix$i$mock_suffix
        FFTname=$FFTdir$FFT$mock_prefix$i$mock_suffix".noweights.grid"$grid".P0"$P0".box"$box
        ./exe/FFT_pth_fkp_now_960grid.exe $Rbox 0 $P0 $nbar_file $mock_file $FFTname & 
        echo $FFTname
        jobcnt=(`jobs -p`)
        if [ ${#jobcnt[@]} -gt 4 ]; then
            wait
        fi 
    done
    if [ $spectrum = "pk" ]; then 
        # compute P(k)
        for i in $(seq -f "%03g" $ni $nf); do
            echo $i
            FFTname=$FFTdir$FFT$mock_prefix$i$mock_suffix".noweights.grid"$grid".P0"$P0".box"$box
            powername=$powerdir$power$mock_prefix$i$mock_suffix".noweights."$Nrandom"randoms.grid"$grid".P0"$P0".box"$box
            ./exe/power_pth_fkp_960grid_480bin.exe $FFTname $FFTrandname $powername $sscale & 
            echo $powername
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 4 ]; then
                wait
            fi 
        done
    elif [ $spectrum = "bk" ]; then 
        # compute B(k)
        for i in $(seq -f "%03g" $ni $nf); do
            echo $i
            FFTname=$FFTdir$FFT$mock_prefix$i$mock_suffix".noweights.grid"$grid".P0"$P0".box"$box
            bispname=$bispdir$bisp$mock_prefix$i$mock_suffix".noweights."$Nrandom"randoms.grid"$grid".nmax"$nmax".nstep"$nstep".P0"$P0".box"$box
            ./exe/bisp_pthalo_960grid.exe $iflag $FFTrandname $FFTname $bispname & 
            echo $bispname
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 4 ]; then
                wait
            fi 
        done
    fi 
elif [ $correction = "wboss" ]; then 
    echo "wboss ONLY" 
    for i in $(seq -f "%03g" $ni $nf); do 
        echo $i
        mock_file=$datadir$mock_prefix$i$mock_suffix
        FFTname=$FFTdir$FFT$mock_prefix$i$mock_suffix".wboss.grid"$grid".P0"$P0".box"$box
        ./exe/FFT_pth_fkp_wboss_960grid.exe $Rbox 0 $P0 $nbar_file $mock_file $FFTname & 
        echo $FFTname
        jobcnt=(`jobs -p`)
        if [ ${#jobcnt[@]} -gt 4 ]; then
            wait
        fi 
    done 
    if [ $spectrum = "pk" ]; then 
        # compute P(k) 
        for i in $(seq -f "%03g" $ni $nf); do 
            echo $i
            FFTname=$FFTdir$FFT$mock_prefix$i$mock_suffix".wboss.grid"$grid".P0"$P0".box"$box
            powername=$powerdir$power$mock_prefix$i$mock_suffix".wboss."$Nrandom"randoms.grid"$grid".P0"$P0".box"$box
            ./exe/power_pth_fkp_960grid_480bin.exe $FFTname $FFTrandname $powername $sscale & 
            echo $powername
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 4 ]; then
                wait 
            fi 
        done
    elif [ $spectrum = "bk" ]; then 
        # compute B(k)
        for i in $(seq -f "%03g" $ni $nf); do
            echo $i
            FFTname=$FFTdir$FFT$mock_prefix$i$mock_suffix".wboss.grid"$grid".P0"$P0".box"$box
            bispname=$bispdir$bisp$mock_prefix$i$mock_suffix".wboss."$Nrandom"randoms.grid"$grid".nmax"$nmax".nstep"$nstep".P0"$P0".box"$box
            ./exe/bisp_pthalo_960grid.exe $iflag $FFTrandname $FFTname $bispname & 
            echo $bispname
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 4 ]; then
                wait
            fi 
        done
    fi
elif [ $correction = "upweight" ]; then 
    echo "Upweight Correction Method" 
    for i in $(seq -f "%03g" $ni $nf); do 
        echo $i
        mock_file=$datadir$mock_prefix$i$mock_suffix
        FFTname=$FFTdir$FFT$mock_prefix$i$mock_suffix".upweight.grid"$grid".P0"$P0".box"$box
        ./exe/FFT_pth_fkp_w_960grid.exe $Rbox 0 $P0 $nbar_file $mock_file $FFTname & 
        echo $FFTname
        jobcnt=(`jobs -p`)
        if [ ${#jobcnt[@]} -gt 4 ]; then
            wait
        fi 
    done 
    if [ $spectrum = "pk" ]; then 
        # compute P(k) 
        for i in $(seq -f "%03g" $ni $nf); do 
            echo $i
            FFTname=$FFTdir$FFT$mock_prefix$i$mock_suffix".upweight.grid"$grid".P0"$P0".box"$box
            powername=$powerdir$power$mock_prefix$i$mock_suffix".upweight."$Nrandom"randoms.grid"$grid".P0"$P0".box"$box
            ./exe/power_pth_fkp_960grid_480bin.exe $FFTname $FFTrandname $powername $sscale & 
            echo $powername
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 4 ]; then
                wait 
            fi 
        done
    elif [ $spectrum = "bk" ]; then 
        # compute B(k)
        for i in $(seq -f "%03g" $ni $nf); do
            echo $i
            FFTname=$FFTdir$FFT$mock_prefix$i$mock_suffix".upweight.grid"$grid".P0"$P0".box"$box
            bispname=$bispdir$bisp$mock_prefix$i$mock_suffix".upweight."$Nrandom"randoms.grid"$grid".nmax"$nmax".nstep"$nstep".P0"$P0".box"$box
            ./exe/bisp_pthalo_960grid.exe $iflag $FFTrandname $FFTname $bispname & 
            echo $bispname
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 4 ]; then
                wait
            fi 
        done
    fi 
elif [ $correction = 'peak' ]; then  
    echo "Peak Correction Method" 
    sigma=5.5       # hardcoded peak sigma value based on first 10 PTHalo mocks
    sigma_suffix=$mock_suffix".peak.sigma"$sigma".fpeak1.0"
    for i in $(seq -f "%03g" $ni $nf); do
        mock_file=$datadir$mock_prefix$i$mock_suffix
        peakcorr_file=$datadir$mock_prefix$i$sigma_suffix
        echo $i
        if [[ ! -a $peakcorr_file ]]; then
            ./exe/fibcollmock_pth_peakcorr.exe $sigma $nbar_file $mock_file $peakcorr_file &
        fi
        jobcnt=(`jobs -p`)
        if [ ${#jobcnt[@]} -gt 4 ]; then
            wait
        fi
    done
    for i in $(seq -f "%03g" $ni $nf); do
        echo $i
        peakcorr_file=$datadir$mock_prefix$i$sigma_suffix
        FFTname=$FFTdir$FFT$mock_prefix$i$sigma_suffix".grid"$grid".P0"$P0".box"$box
        ./exe/FFT_pth_fkp_w_960grid.exe $Rbox 0 $P0 $nbar_file $peakcorr_file $FFTname &
        echo $FFTname
        jobcnt=(`jobs -p`)
        if [ ${#jobcnt[@]} -gt 4 ]; then
            wait
        fi
    done
    if [ $spectrum = "pk" ]; then 
        for i in $(seq -f "%03g" $ni $nf); do
            echo $i
            FFTname=$FFTdir$FFT$mock_prefix$i$sigma_suffix".grid"$grid".P0"$P0".box"$box
            powername=$powerdir$power$mock_prefix$i$sigma_suffix"."$Nrandom"randoms.grid"$grid".P0"$P0".box"$box
            ./exe/power_pth_fkp_960grid_480bin.exe $FFTname $FFTrandname $powername $sscale &
            echo $powername
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 4 ]; then
                wait
            fi
        done
    elif [ $spectrum = "bk" ]; then 
        # compute B(k)
        for i in $(seq -f "%03g" $ni $nf); do
            echo $i
            FFTname=$FFTdir$FFT$mock_prefix$i$sigma_suffix".grid"$grid".P0"$P0".box"$box
            bispname=$bispdir$bisp$mock_prefix$i$sigma_suffix"."$Nrandom"randoms.grid"$grid".nmax"$nmax".nstep"$nstep".P0"$P0".box"$box
            ./exe/bisp_pthalo_960grid.exe $iflag $FFTrandname $FFTname $bispname & 
            echo $bispname
            jobcnt=(`jobs -p`)
            if [ ${#jobcnt[@]} -gt 4 ]; then
                wait
            fi 
        done
    fi 
else 
    echo "Enter correct correction method: noweights, upweight, peak. Entry is case sensititve!"
fi 
