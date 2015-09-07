#!/bin/bash
echo "Taring CMASS LOWZ combined catalog for Ariel"

datadir="/mount/riachuelo1/hahn/data/CMASS/dr12v5/"
cd $datadir
# tar: (c)reate g(z)ip (v)erbose (f)ile [filename.tar.gz] [contents]...
tar -zcvf "DR12v5_CMASSLOWZ_LOW.tar.gz" "galaxy_DR12v5_CMASSLOWZ_LOW_North.fidcosmo.dat" "galaxy_DR12v5_CMASSLOWZE2_LOW_North.fidcosmo.dat" "galaxy_DR12v5_CMASSLOWZE3_LOW_North.fidcosmo.dat" "galaxy_DR12v5_CMASSLOWZ_LOW_North.fidcosmo.gauss.dlospeak.sigma6.9.fpeak0.72.dat" "galaxy_DR12v5_CMASSLOWZE2_LOW_North.fidcosmo.gauss.dlospeak.sigma6.9.fpeak0.72.dat" "galaxy_DR12v5_CMASSLOWZE3_LOW_North.fidcosmo.gauss.dlospeak.sigma6.9.fpeak0.72.dat" "random0_DR12v5_CMASSLOWZ_LOW_North.ran.dat" "random0_DR12v5_CMASSLOWZE2_LOW_North.ran.dat" "random0_DR12v5_CMASSLOWZE3_LOW_North.ran.dat"

tar -zcvf "DR12v5_CMASSLOWZ_HIGH.tar.gz" "galaxy_DR12v5_CMASSLOWZ_HIGH_North.fidcosmo.dat" "galaxy_DR12v5_CMASSLOWZE2_HIGH_North.fidcosmo.dat" "galaxy_DR12v5_CMASSLOWZE3_HIGH_North.fidcosmo.dat" "galaxy_DR12v5_CMASSLOWZ_HIGH_North.fidcosmo.gauss.dlospeak.sigma6.3.fpeak0.7.dat" "galaxy_DR12v5_CMASSLOWZE2_HIGH_North.fidcosmo.gauss.dlospeak.sigma6.3.fpeak0.7.dat" "galaxy_DR12v5_CMASSLOWZE3_HIGH_North.fidcosmo.gauss.dlospeak.sigma6.3.fpeak0.7.dat" "random0_DR12v5_CMASSLOWZ_HIGH_North.ran.dat" "random0_DR12v5_CMASSLOWZE2_HIGH_North.ran.dat" "random0_DR12v5_CMASSLOWZE3_HIGH_North.ran.dat"
