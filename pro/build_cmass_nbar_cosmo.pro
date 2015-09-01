pro build_cmass_nbar_cosmo, fiducial=fiducial 
;NAME:
;   build_cmass_nbar_cosmo
;PURPOSE:
;   Build corrected nbar(z) file for CMASS DR12v4 depending on the choice of cosmology 
;CALLING SEQUENCE:
;
;OUTPUTS:
;
;COMMENTS:
;
;REVISION HISTORY:
;13-07-2015  Written by ChangHoon Hahn, NYU.
    ; Read in cmass nbar(z) file
    data_dir = '/mount/riachuelo1/hahn/data/'
    data_file = data_dir + 'nbar-cmass-dr12v4-N-Reid.dat'
    readcol, data_file, z_mid, z_low, z_high, nbar, wfkp, shell_vol, ngal_tot

    survey_sqrdeg = 6.851419E+03  ; square degrees of dr12v4
    total_sqrdeg = 41253.0  ; square degrees of entire sky 
    survey_frac = survey_sqrdeg / total_sqrdeg
    print, 'Survey Fraction = ', survey_frac

    if keyword_set(fiducial) then begin 
        omega_m = 0.31
        omega_l = 1.0 - omega_m 
        file_str = 'fidcosmo'
    endif else begin 
        print, 'Not Yet Implemented'
        return 
    endelse 
    
    P0 = 20000.0
    new_nbar = fltarr(n_elements(z_mid))
    new_wfkp = fltarr(n_elements(z_mid))
    new_shell = fltarr(n_elements(z_mid))
    for iz=0L,n_elements(z_mid)-1L do begin 
        ; comoving volume of redshift shell 
        V_shell = (4.0/3.0) * !PI * survey_frac * (lf_comvol(z_high[iz], omega0=omega_m, omegal0=omega_l) - $
            lf_comvol(z_low[iz], omega0=omega_m, omegal0=omega_l)) 

        new_nbar[iz] = ngal_tot[iz]/V_shell 
        new_wfkp[iz] = 1.0/(1.0 + new_nbar[iz]*P0)
        new_shell[iz] = V_shell
    endfor 
    
    ; write nbar(z) file 
    nbar_file = data_dir + 'nbar-cmass-dr12v4-N-Reid-'+file_str+'.dat'
    openw, lun, nbar_file, /get_lun 
    for iz=0L,n_elements(z_mid)-1L do begin 
        printf, lun, z_mid[iz], z_low[iz], z_high[iz], $
            new_nbar[iz], new_wfkp[iz], new_shell[iz], ngal_tot[iz], $
            format='(f,f,f,e11.4,e11.4,e11.4,e11.4)'
    endfor 
    free_lun, lun 
return 
end 
