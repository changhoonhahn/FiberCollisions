pro build_fibcoll_catalog, catalog, mock_file, mock_dlos_file
; build catalog data that includes line-of-sight separation distances, upweighted neighbor, etc for Mock Catalogs
; spherematch is used to compare upweighted galaxies with unweighted galaxies
; USING PROPER COSMOLOGY!
    fib_angscale = 0.01722
    print, mock_file     

    ; currently only coded for QPM (which already has wcp values) 
    if (strlowcase(catalog) EQ 'qpm') then begin 
        omega_m = 0.31
        readcol, mock_file, mock_ra, mock_dec, mock_redshift, w_fkp, w_cp, w_comp

        n_gal = n_elements(w_comp) 

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
    endif 
    gal_upcp = m_upcp[uniq(m_upcp, sort(m_upcp))]

    d_fc = fltarr(n_gal)    ; LOS displacement between fiber-collided pairs
    ; Fiber-collided pair redshift, ra, and dec (for galaxies that have been downweighted)  
    fc_pair_z   = fltarr(n_gal)  
    fc_pair_ra  = fltarr(n_gal)  
    fc_pair_dec = fltarr(n_gal)  
    fc_pair_index = lonarr(n_gal)

    for i=0L,n_elements(gal_upcp)-1L do begin
        ngal=gal_upcp[i]
        collision_indx = where( m_upcp eq ngal )

        targ_ra  = ra_upcp[ngal] 
        targ_dec = dec_upcp[ngal] 
        targ_red = red_upcp[ngal]           ; target redshift
        targ_Dc  = 3000.0*comdis(targ_red, omega_m, 1.0-omega_m)         ; Mpc/h 

        neigh_ra  = ra_nocp[m_nocp[collision_indx]]
        neigh_dec = dec_nocp[m_nocp[collision_indx]]
        neigh_red = red_nocp[m_nocp[collision_indx]]    ; neighbor redshift(s)
        neigh_Dc  = 3000.0*comdis(neigh_red, omega_m, 1.0-omega_m)       ; Mpc/h

        n_neigh = n_elements(neigh_red)     ; number of neighbors 
        for j=0L,n_neigh-1L do begin
            d_fc[nocp_indx[m_nocp[collision_indx[j]]]] = neigh_Dc[j] - targ_Dc
            fc_pair_z[nocp_indx[m_nocp[collision_indx[j]]]]   = targ_red
            fc_pair_ra[nocp_indx[m_nocp[collision_indx[j]]]]  = targ_ra
            fc_pair_dec[nocp_indx[m_nocp[collision_indx[j]]]] = targ_dec
            fc_pair_index[nocp_indx[m_nocp[collision_indx[j]]]] = upcp_indx[ngal]
        endfor 
    endfor  

    ;if ((strlowcase(catalog) EQ 'lasdamasgeo') or (strlowcase(catalog) EQ 'tilingmock')) then begin 
    ;    ; to account for double counting the fc pairs 
    ;    print, '#dlos with repeats', n_elements(all_dlos)
    ;    abs_dlos = abs(all_dlos)
    ;    sorted_indx = sort(abs_dlos)
    ;    sorted_abs_dlos = abs_dlos[sorted_indx]
    ;    uniq_sorted_indx = uniq(sorted_abs_dlos)

    ;    all_dlos = all_dlos[sorted_indx[uniq_sorted_indx]]          ; sorts the dlos but who cares
    ;    print, '#dlos without repeats', n_elements(all_dlos)
    ;endif 
    openw, lun, mock_dlos_file, /get_lun
    for i_gal=0L,n_gal-1L do begin 
        printf, lun, mock_ra[i_gal], mock_dec[i_gal], mock_redshift[i_gal], $
            w_fkp[i_gal], w_cp[i_gal], w_comp[i_gal],$
            d_fc[i_gal], fc_pair_ra[i_gal], fc_pair_dec[i_gal], fc_pair_z[i_gal], fc_pair_index[i_gal],$
            format='(f,f,f,f,f,f,f,f,f,f,i)'
    endfor 
    free_lun, lun 
return 
end 
;else if (strlowcase(catalog) EQ 'tilingmock') then begin 
;    omega_m = 0.274
;    readcol, mock_file, mock_ra, mock_dec, mock_redshift, mock_w
;    
;    upcp_indx = range(0,n_elements(mock_ra)-1L) 
;    nocp_indx = range(0,n_elements(mock_ra)-1L) 
;    ra_upcp  = mock_ra
;    dec_upcp = mock_dec
;    red_upcp = mock_redshift
;                                                                                                               
;    ra_nocp  = mock_ra
;    dec_nocp = mock_dec
;    red_nocp = mock_redshift
;                                                                                                               
;    spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, match1, match2, d12, maxmatch=0     
;                                                                                                               
;    overlap_indx = where(d12 GT 0.0, overlapcount)
;    m_nocp = match1[overlap_indx]
;    m_upcp = match2[overlap_indx]
;endif else if (strlowcase(catalog) EQ 'lasdamasgeo') then begin 
;    omega_m = 0.25
;    readcol, mock_file, mock_ra, mock_dec, mock_redshift
;    mock_redshift = temporary(mock_redshift)/299800.0
;                                                                                                               
;    upcp_indx = range(0,n_elements(mock_ra)-1L)
;    nocp_indx = range(0,n_elements(mock_ra)-1L)
;    ra_upcp  = mock_ra
;    dec_upcp = mock_dec
;    red_upcp = mock_redshift
;                                                                                                               
;    ra_nocp  = mock_ra
;    dec_nocp = mock_dec
;    red_nocp = mock_redshift
;                                                                                                               
;    spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, match1, match2, d12, maxmatch=0
;                                                                                                               
;    overlap_indx = where(d12 GT 0.0, overlapcount)
;    m_nocp = match1[overlap_indx]
;    m_upcp = match2[overlap_indx]
;endif else if (strlowcase(catalog) EQ 'cmass') then begin 
;    omega_m = 0.274
;    readcol, mock_file, mock_ra, mock_dec, mock_redshift, mock_wsys, mock_wnoz, mock_wfc, mock_nbar, mock_comp
;                                                                                                               
;    upcp_indx = where(mock_wfc GT 1.) 
;    ra_upcp  = mock_ra[upcp_indx]
;    dec_upcp = mock_dec[upcp_indx]
;    red_upcp = mock_redshift[upcp_indx]
;                                                                                                               
;    nocp_indx = range(0,n_elements(mock_ra)-1L)
;    ra_nocp  = mock_ra 
;    dec_nocp = mock_dec
;    red_nocp = mock_redshift
;                                                                                                               
;    spherematch, ra_nocp, dec_nocp, ra_upcp, dec_upcp, fib_angscale, match1, match2, d12, maxmatch=0
;    
;    overlap_indx = where(d12 GT 0.0, overlapcount)
;    m_nocp = match1[overlap_indx]
;    m_upcp = match2[overlap_indx]
;endif 
