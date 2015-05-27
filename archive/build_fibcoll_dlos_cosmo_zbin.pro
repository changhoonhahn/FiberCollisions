pro build_fibcoll_dlos_cosmo_zbin
; build line-of-sight separation distances for fibercollided pairs in Mock Catalogs
; spherematch is used to compare upweighted galaxies with unweighted galaxies
; USING PROPER COSMOLOGY! AND DIVIDED INTO REDSHIFT BINS
; to look at the redshift dependence of d_LOS
; only coded for QPM 
    fib_angscale = 0.01722
   
    data_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'
    for i_file=1,10L do begin 
        mock_file = data_dir+'a0.6452_'+strmid(strtrim(string(i_file+10000),2),1)+'.dr12d_cmass_ngc.vetoed.fibcoll.dat' 
        print, mock_file     
        ;if (strlowcase(catalog) NE 'qpm') then STOP  
            
        readcol, mock_file, mock_ra, mock_dec, mock_redshift, w_fkp, w_cp
        
        omega_m = 0.31
        for i_z=0L,2L do begin 
            if (i_z EQ 0) then begin 
                zmin = 0.43
                zmax = 0.5

                zbin_str = 'low' 
            endif else if (i_z EQ 1L) then begin 
                zmin = 0.5
                zmax = 0.6

                zbin_str = 'mid' 
            endif else begin 
                zmin = 0.6
                zmax = 0.7

                zbin_str = 'high' 
            endelse
            dlos_file = data_dir+'dLOS_a0.6452_'+strmid(strtrim(string(i_file+10000),2),1)+'.dr12d_cmass_ngc.vetoed.fibcoll.'+zbin_str+'.dat'
        
            zbin_index = where((mock_redshift GE zmin) AND (mock_redshift LT zmax), n_zbin) 

            zbin_ra         = mock_ra[zbin_index]
            zbin_dec        = mock_dec[zbin_index]
            zbin_redshift   = mock_redshift[zbin_index]
            zbin_wfkp       = w_fkp[zbin_index]
            zbin_wcp        = w_cp[zbin_index]

            ;Galaxies with up-weighted w_cp
            upcp_indx = where(zbin_wcp gt 1)
            ra_upcp  = zbin_ra[upcp_indx]
            dec_upcp = zbin_dec[upcp_indx]
            red_upcp = zbin_redshift[upcp_indx]

            nocp_indx = where(zbin_wcp eq 0)
            ra_nocp  = zbin_ra[nocp_indx]
            dec_nocp = zbin_dec[nocp_indx]
            red_nocp = zbin_redshift[nocp_indx]

            spherematch, ra_nocp, dec_nocp, ra_upcp, dec_upcp, fib_angscale, m_nocp, m_upcp, d12, maxmatch=0

            gal_upcp = m_upcp[uniq(m_upcp, sort(m_upcp))]
            print, dlos_file 

            openw, lun, dlos_file, /get_lun

            for i=0L,n_elements(gal_upcp)-1L do begin
                ngal=gal_upcp[i]
                collision_indx = where( m_upcp eq ngal )

                targ_indx = upcp_indx[ngal]
                targ_phi = ra_upcp[ngal]
                targ_theta = 90.0-dec_upcp[ngal]
                targ_red = red_upcp[ngal]
                targ_r = 3000.0*comdis(targ_red, omega_m, 1.0-omega_m)

                neigh_indx = nocp_indx[m_nocp[collision_indx]]
                neigh_phi = ra_nocp[m_nocp[collision_indx]]
                neigh_theta = 90.0-dec_nocp[m_nocp[collision_indx]]
                neigh_red = red_nocp[m_nocp[collision_indx]]
                neigh_r = 3000.0*comdis(neigh_red, omega_m, 1.0-omega_m)

                angles_to_xyz, targ_r, targ_phi, targ_theta, targ_x, targ_y, targ_z
                angles_to_xyz, neigh_r, neigh_phi, neigh_theta, neigh_x, neigh_y, neigh_z

                targ_mag= targ_x^2+targ_y^2+targ_z^2
                targ_dot_neigh= targ_x*neigh_x+targ_y*neigh_y+targ_z*neigh_z
                proj = targ_dot_neigh/targ_mag

                LOS_x = proj*targ_x
                LOS_y = proj*targ_y
                LOS_z = proj*targ_z
                LOS_mag = LOS_x^2+LOS_y^2+LOS_z^2

                d_los = dblarr(n_elements(LOS_mag))
                for j=0L,n_elements(LOS_mag)-1L do begin
                    if LOS_mag[j] ge targ_mag then begin
                        d_los[j] = sqrt((LOS_x[j]-targ_x)^2+(LOS_y[j]-targ_y)^2+(LOS_z[j]-targ_z)^2)
                    endif else begin
                        d_los[j] = -sqrt((LOS_x[j]-targ_x)^2+(LOS_y[j]-targ_y)^2+(LOS_z[j]-targ_z)^2)
                    endelse
                    printf, lun, d_los[j], format='(f)'
                endfor
            endfor  
            free_lun, lun 
        endfor 
    endfor 
return 
end 
;    if (strlowcase(catalog) EQ 'pthalo') then begin 
;        omega_m = 0.274
;
;        readcol,mock_file, mock_ra, mock_dec, mock_redshift, ipoly, $
;            w_boss, w_cp, w_red, redtrue, flag, m1, m2, veto
;
;        ;Galaxies with up-weighted w_cp
;        upcp_indx = where(w_cp gt 1)
;        ra_upcp  = mock_ra[upcp_indx]
;        dec_upcp = mock_dec[upcp_indx]
;        red_upcp = mock_redshift[upcp_indx]
;
;        nocp_indx = where(w_cp eq 0)
;        ra_nocp  = mock_ra[nocp_indx]
;        dec_nocp = mock_dec[nocp_indx]
;        red_nocp = mock_redshift[nocp_indx]
;
;        spherematch, ra_nocp, dec_nocp, ra_upcp, dec_upcp, fib_angscale, m_nocp, m_upcp, d12, maxmatch=0
;    endif else if (strlowcase(catalog) EQ 'tilingmock') then begin 
;        omega_m = 0.274
;        readcol, mock_file, mock_ra, mock_dec, mock_redshift, mock_w
;        
;        upcp_indx = range(0,n_elements(mock_ra)-1L) 
;        nocp_indx = range(0,n_elements(mock_ra)-1L) 
;        ra_upcp  = mock_ra
;        dec_upcp = mock_dec
;        red_upcp = mock_redshift
;
;        ra_nocp  = mock_ra
;        dec_nocp = mock_dec
;        red_nocp = mock_redshift
;
;        spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, match1, match2, d12, maxmatch=0     
;
;        overlap_indx = where(d12 GT 0.0, overlapcount)
;        m_nocp = match1[overlap_indx]
;        m_upcp = match2[overlap_indx]
;    endif else if (strlowcase(catalog) EQ 'lasdamas') then begin 
;        omega_m = 0.25
;        readcol, mock_file, mock_ra, mock_dec, mock_redshift
;        mock_redshift = temporary(mock_redshift)/299800.0
;
;        upcp_indx = range(0,n_elements(mock_ra)-1L)
;        nocp_indx = range(0,n_elements(mock_ra)-1L)
;        ra_upcp  = mock_ra
;        dec_upcp = mock_dec
;        red_upcp = mock_redshift
;
;        ra_nocp  = mock_ra
;        dec_nocp = mock_dec
;        red_nocp = mock_redshift
;
;        spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, match1, match2, d12, maxmatch=0
;
;        overlap_indx = where(d12 GT 0.0, overlapcount)
;        m_nocp = match1[overlap_indx]
;        m_upcp = match2[overlap_indx]



;    endif else if (strlowcase(catalog) EQ 'cmass') then begin 
;        omega_m = 0.274
;        readcol, mock_file, mock_ra, mock_dec, mock_redshift, mock_wsys, mock_wnoz, mock_wfc, mock_nbar, mock_comp
;
;        upcp_indx = where(mock_wfc GT 1.) 
;        ra_upcp  = mock_ra[upcp_indx]
;        dec_upcp = mock_dec[upcp_indx]
;        red_upcp = mock_redshift[upcp_indx]
;
;        nocp_indx = range(0,n_elements(mock_ra)-1L)
;        ra_nocp  = mock_ra 
;        dec_nocp = mock_dec
;        red_nocp = mock_redshift
;
;        spherematch, ra_nocp, dec_nocp, ra_upcp, dec_upcp, fib_angscale, match1, match2, d12, maxmatch=0
;        
;        overlap_indx = where(d12 GT 0.0, overlapcount)
;        m_nocp = match1[overlap_indx]
;        m_upcp = match2[overlap_indx]
;    endif 
