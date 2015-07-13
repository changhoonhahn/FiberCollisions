pro build_cmass_nbar_cosmo, fiducial=fiducial 
; Build corrected nbar(z) file for CMASS DR12v4 depending on the choice of cosmology 

    ; Read in cmass nbar(z) fil 
    data_dir = '/mount/riachuelo1/hahn/data/'
    data_file = data_dir + 'nbar-cmass-dr12v4-N-Reid.dat'
    readcol, data_file, z_mid, z_low, z_high, nbar, wfkp, shell_voll, ngal_tot
    return 
end 
