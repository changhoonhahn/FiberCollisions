function z_to_v, z 
    c_speed = 299792.0
    return, (c_speed*z^2+2*c_speed*z)/(z^2+2*z+2)
end 

function v_pec, z1, z2
    c_speed = 299792.0
    return, c_speed*(z2-z1)/(1.0+z1) 
end

pro build_fibcoll_vlos, catalog, mock_file, vlos_file
; build line-of-sight separation distances for fibercollided pairs in Mock Catalogs
; spherematch is used to compare upweighted galaxies with unweighted galaxies
    fib_angscale = 0.01722
    c_speed = 299792.0 
    print, mock_file     
    if (strlowcase(catalog) EQ 'pthalo') then begin 
        readcol,mock_file, mock_ra, mock_dec, mock_redshift, ipoly, $
            w_boss, w_cp, w_red, redtrue, flag, m1, m2, veto

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
    endif else if (strlowcase(catalog) EQ 'tilingmock') then begin 
        readcol, mock_file, mock_ra, mock_dec, mock_redshift, mock_w
        
        upcp_indx = range(0,n_elements(mock_ra)-1L) 
        nocp_indx = range(0,n_elements(mock_ra)-1L) 
        ra_upcp  = mock_ra
        dec_upcp = mock_dec
        red_upcp = mock_redshift

        ra_nocp  = mock_ra
        dec_nocp = mock_dec
        red_nocp = mock_redshift

        spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, match1, match2, d12, maxmatch=0     

        overlap_indx = where(d12 GT 0.0, overlapcount)
        m_nocp = match1[overlap_indx]
        m_upcp = match2[overlap_indx]
    endif else if (strlowcase(catalog) EQ 'lasdamasgeo') then begin 
        readcol, mock_file, mock_ra, mock_dec, mock_redshift

        upcp_indx = range(0,n_elements(mock_ra)-1L)
        nocp_indx = range(0,n_elements(mock_ra)-1L)
        ra_upcp  = mock_ra
        dec_upcp = mock_dec
        red_upcp = mock_redshift/c_speed

        ra_nocp  = mock_ra
        dec_nocp = mock_dec
        red_nocp = mock_redshift/c_speed

        spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, match1, match2, d12, maxmatch=0

        overlap_indx = where(d12 GT 0.0, overlapcount)
        m_nocp = match1[overlap_indx]
        m_upcp = match2[overlap_indx]
        print, n_elements(m_upcp) 
    endif else if (strlowcase(catalog) EQ 'qpm') then begin 
        readcol,mock_file, mock_ra, mock_dec, mock_redshift, w_fkp, w_cp

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
    endif else if (strlowcase(catalog) EQ 'cmass') then begin 
        readcol, mock_file, mock_ra, mock_dec, mock_redshift, mock_wsys, mock_wnoz, mock_wfc, mock_nbar, mock_comp

        upcp_indx = where(mock_wfc GT 1.) 
        ra_upcp  = mock_ra[upcp_indx]
        dec_upcp = mock_dec[upcp_indx]
        red_upcp = mock_redshift[upcp_indx]

        nocp_indx = range(0,n_elements(mock_ra)-1L)
        ra_nocp  = mock_ra 
        dec_nocp = mock_dec
        red_nocp = mock_redshift

        spherematch, ra_nocp, dec_nocp, ra_upcp, dec_upcp, fib_angscale, match1, match2, d12, maxmatch=0
        
        overlap_indx = where(d12 GT 0.0, overlapcount)
        m_nocp = match1[overlap_indx]
        m_upcp = match2[overlap_indx]
    endif 

    gal_upcp = m_upcp[uniq(m_upcp, sort(m_upcp))]
    
    print, vlos_file 
    openw, lun, vlos_file, /get_lun

    for i=0L,n_elements(gal_upcp)-1L do begin
        ngal=gal_upcp[i]
        collision_indx = where( m_upcp eq ngal )

        targ_red = red_upcp[ngal]

        neigh_red = red_nocp[m_nocp[collision_indx]]

        v_los = dblarr(n_elements(neigh_red))
        for j=0L,n_elements(neigh_red)-1L do begin
            printf, lun, c_speed*(neigh_red[j]-targ_red), z_to_v(neigh_red[j])-z_to_v(targ_red), v_pec(targ_red, neigh_red[j]), $
                format='(f, f, f)'
        endfor
    endfor  
    free_lun, lun 
return 
end 
