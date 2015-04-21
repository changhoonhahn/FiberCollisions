pro build_nbar, tilingmock=tilingmock 
; calcultes nbar(z) for tiling mock galaxy mock catalog using comoving volume
; calculations of CMASS nbar(z) in order to consistently compare tiling mock 
; P(k) with CMASS P(k)
    datadir = '/mount/riachuelo1/hahn/data/manera_mock/dr11/'
    cmass_nbar = 'nbar-cmass-dr11may22-N-Anderson.dat'
    readcol, datadir+cmass_nbar, z_mid,z_low,z_high,nbarz,wfkp,shell_vol,gal_tot
    cmass_area = 6.3077345956682420E+03     ; hardcoded square degree area of the survey

    if keyword_set(tilingmock) then begin 
    ; import tiling mock geometry
        tm_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'
        tm_poly = 'geometry-boss5003-icoll012.ply'
        read_mangle_polygons, tm_dir+tm_poly, tm_polygon
        tm_area = total(tm_polygon.str, /double)*(180.0/!PI)^2 ; in square degrees
        f_area = tm_area/cmass_area     ; factor to scale the shell volumes 
        tm_shell_vol = f_area*shell_vol

    ; import tiling mock data 
        tm_file = 'cmass-boss5003sector-icoll012.dat'
        readcol, tm_dir+tm_file, ra, dec, z, wcomp 

    ; calculate tiling mock nbar(z) 
        tm_nbar_file = 'nbar-'+tm_file
        openw, lun, tm_dir+tm_nbar_file, /get_lun
        for i=0L,n_elements(z_mid)-1L do begin 
            zbin = where((z GE z_low[i]) AND (z LT z_high[i]),n_zbin)   ; number of galaxies in each bin 
            if (n_zbin NE 0) then ngal_zbin = total(wcomp[zbin]) $
                else ngal_zbin = 0.0 
            tm_nbarz = float(ngal_zbin)/tm_shell_vol[i]
            printf, lun, z_mid[i],z_low[i],z_high[i],tm_nbarz,float(wfkp[i]),float(tm_shell_vol[i]),ngal_zbin, $
                format='(4F,3E)'
        endfor
        free_lun, lun 
    endif 
return 
end 
