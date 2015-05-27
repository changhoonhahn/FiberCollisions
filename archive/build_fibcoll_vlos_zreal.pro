function z_to_v, z 
    c_speed = 299792.0
    return (c_speed*z^2+2*c_speed*z)/(z^2+2*z+2)
end 

pro build_fibcoll_vlos_zreal, catalog, mock_file, vlos_file
; build line-of-sight separation distances for fibercollided pairs in Mock Catalogs
; spherematch is used to compare upweighted galaxies with unweighted galaxies
    fib_angscale = 0.01722
    c_speed = 299792.0 
    print, mock_file     
    if (strlowcase(catalog) EQ 'qpm') then begin 
        omega_m = 0.31 
        readcol,mock_file, mock_ra, mock_dec, mock_redshift, w_fkp, w_cp

        ; read in mock .info file
        readcol, mock_file+'.info', junk1, junk2, mock_redshift_real, junk3, junk4, junk5 ;, skipline=3

        if (n_elements(junk1) NE n_elements(mock_ra)) then begin
            print, 'not matching'
            stop
        endif 

        ;Galaxies with up-weighted w_cp
        upcp_indx = where(w_cp gt 1)
        ra_upcp  = mock_ra[upcp_indx]
        dec_upcp = mock_dec[upcp_indx]
        red_upcp = mock_redshift[upcp_indx]

        nocp_indx = where(w_cp eq 0)
        ra_nocp  = mock_ra[nocp_indx]
        dec_nocp = mock_dec[nocp_indx]
        red_nocp = mock_redshift[nocp_indx]

        spherematch, ra_nocp, dec_nocp, ra_upcp, dec_upcp, fib_angscale, m_nocp, m_upcp, d12, maxmatch=0
    endif else begin 
        print, 'blaksdjflakjsdfkjasdf'
        stop
    endelse

    gal_upcp = m_upcp[uniq(m_upcp, sort(m_upcp))]
    
    print, vlos_file 
    openw, lun, vlos_file, /get_lun

    for i=0L,n_elements(gal_upcp)-1L do begin
        ngal=gal_upcp[i]
        collision_indx = where( m_upcp eq ngal )

        targ_indx = upcp_indx[ngal]
        targ_phi = ra_upcp[ngal]
        targ_theta = 90.0-dec_upcp[ngal]
        targ_red = red_upcp[ngal]

        neigh_indx = nocp_indx[m_nocp[collision_indx]]
        neigh_phi = ra_nocp[m_nocp[collision_indx]]
        neigh_theta = 90.0-dec_nocp[m_nocp[collision_indx]]
        neigh_red = red_nocp[m_nocp[collision_indx]]

        v_los = dblarr(n_elements(neigh_red))
        for j=0L,n_elements(neigh_red)-1L do begin
            v_los[j] = c_speed*(neigh_red[j]-targ_red)
            printf, lun, v_los[j], format='(f)'
        endfor
    endfor  
    free_lun, lun 
return 
end 
