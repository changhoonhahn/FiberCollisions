#!/bin/bash
mock_file="cmass-boss5003sector-icoll012.dat"
rand_file="randoms-boss5003-icoll012-vetoed.dat"
nbar_file="nbar-cmass-boss5003sector-icoll012.dat"

P0=20000; sscale=4000.0; Rbox=2000.0; box="4000"; FFT="FFT_"; power="power_"; grid="360"; BISP="BISP_"
# bispectrum inputs 
nstep="3"; ncut="3"; iflag=2
# directories
datadir="/mount/riachuelo1/hahn/data/tiling_mocks/"     # data directory 
randdir="/mount/riachuelo1/hahn/data/tiling_mocks/"     # random data directory 
FFTdir="/mount/riachuelo1/hahn/FFT/tiling_mocks/"       # FFT directory
powerdir="/mount/riachuelo1/hahn/power/tiling_mocks/"   # P(k) directory 
bispdir="/mount/riachuelo1/hahn/bispec/tiling_mock/"

# bash inputs
read -p "Correction method (true, delta, or peak): " -e correction  
read -p "Powerspectrum or Bispectrum (pk or bk): " -e spectrum  

# compile FORTRAN files
ifort -O3 -o exe/peakcorr_tilingmock_fibcollmock.exe peakcorr_tilingmock_fibcollmock.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/FFT-fkp-tilingmock-360grid.exe FFT-fkp-tilingmock-360grid.f -L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw
ifort -O3 -o exe/power-fkp-tilingmock-360grid-180bin.exe power-fkp-tilingmock-360grid-180bin.f
ifort -O4 -o exe/bisp_tilingmock_360grid.exe bisp_tilingmock_360grid.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw

# random FFT
FFTrandname=$FFTdir$FFT$rand_file".grid"$grid".P0"$P0".box"$box  # same FFT for all correction method   
if [[ ! -a $FFTrandname ]]; then
    ./exe/FFT-fkp-tilingmock-360grid.exe $Rbox 1 $P0 $datadir$nbar_file $randdir$rand_file $FFTrandname
    echo $FFTrandname
fi

if [ $correction = "true" ]; then 
    # Compute the true P(k) for Tiling Mock using no weight mocks 
    echo "True Correction Method" 
    # Compute FFT 
    FFTname=$FFTdir$FFT$mock_file".grid"$grid".P0"$P0".box"$box
    ./exe/FFT-fkp-tilingmock-360grid.exe $Rbox 0 $P0 $datadir$nbar_file $datadir$mock_file $FFTname 
    echo $FFTname

    if [ $spectrum = "pk" ]; then 
        # Compute P(k)  
        powername=$powerdir$power$mock_file".grid"$grid".P0"$P0".box"$box
        ./exe/power-fkp-tilingmock-360grid-180bin.exe $FFTname $FFTrandname $powername $sscale
        echo $powername
    elif [ $spectrum = "bk" ]; then 
        # compute B(k) 
        bispname=$bispdir$BISP$mock_file".grid"$grid".nmax"$nmax".nstep"$nstep".P0"$P0".box"$box
        ./exe/bisp_tilingmock_360grid.exe $iflag $FFTrandname $FFTname $bispname
    fi
elif [ $correction = "delta" ] || [ $correction = "peak" ]; then 
    # build wcp assigned mocks 
    fibcoll_file="cmass-boss5003sector-icoll012.fibcoll.dat"
    if [[ ! -a $datadir$fibcoll_file ]]; then       
        echo "Building w_cp assigned mocks" 
        idl -e "build_mock_wcp_assign"
    fi 
    echo $datadir$fibcoll_file

    if [ $correction = "delta" ]; then 
        echo "Delta Correction Method" 
        # Compute FFT
        FFTname=$FFTdir$FFT$fibcoll_file".delta.grid"$grid".P0"$P0".box"$box
        ./exe/FFT-fkp-tilingmock-360grid.exe $Rbox 0 $P0 $datadir$nbar_file $datadir$fibcoll_file $FFTname
        echo $FFTname
        if [ $spectrum = "pk" ]; then 
            # Compute P(k) 
            powername=$powerdir$power$fibcoll_file".delta.grid"$grid".P0"$P0".box"$box
            ./exe/power-fkp-tilingmock-360grid-180bin.exe $FFTname $FFTrandname $powername $sscale
            echo $powername
        elif [ $spectrum = "bk" ]; then 
            # compute B(k)
            bispname=$bispdir$BISP$fibcoll_file".delta.grid"$grid".nmax"$nmax".nstep"$nstep".P0"$P0".box"$box
            ./exe/bisp_tilingmock_360grid.exe $iflag $FFTrandname $FFTname $bispname
        fi
    # Peak Correction Method
    else 
        echo "Peak Correction Method" 
        sigma=5.98       # hardcoded peak sigma value
        sigma_suffix=".peak.sigma"$sigma".fpeak1.0"

        # Build peak corrected mock
        peakcorr_file=$datadir$fibcoll_file$sigma_suffix
        if [[ ! -a $peakcorr_file ]]; then
            ./exe/peakcorr_tilingmock_fibcollmock.exe $sigma $datadir$nbar_file $datadir$fibcoll_file $peakcorr_file 
            echo $peakcorr_file 
        fi
        # compute FFT
        FFTname=$FFTdir$FFT$fibcoll_file$sigma_suffix".grid"$grid".P0"$P0".box"$box
        ./exe/FFT-fkp-tilingmock-360grid.exe $Rbox 0 $P0 $datadir$nbar_file $peakcorr_file $FFTname 
        echo $FFTname
        if [ $spectrum = "pk" ]; then 
           # compute P(k)  
            powername=$powerdir$power$fibcoll_file$sigma_suffix".grid"$grid".P0"$P0".box"$box
            ./exe/power-fkp-tilingmock-360grid-180bin.exe $FFTname $FFTrandname $powername $sscale
            echo $powername
        elif [ $spectrum = "bk" ]; then 
            # compute B(k) 
            bispname=$bispdir$BISP$fibcoll_file$sigma_suffix".grid"$grid".nmax"$nmax".nstep"$nstep".P0"$P0".box"$box
            ./exe/bisp_tilingmock_360grid.exe $iflag $FFTrandname $FFTname $bispname
        fi
    fi
else 
    echo "Enter correct correction method: true, delta, peak. Entry is case sensititve!"
fi 

