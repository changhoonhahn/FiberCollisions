pro build_fibcoll_dlos_cosmo, catalog, mock_file, dlos_file
; build line-of-sight separation distances for fibercollided pairs in Mock Catalogs
; spherematch is used to compare upweighted galaxies with unweighted galaxies
; USING PROPER COSMOLOGY!

    fib_angscale = 0.01722
    print, mock_file     
    if (strlowcase(catalog) EQ 'pthalo') then begin 
        omega_m = 0.274

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
        omega_m = 0.274
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
        omega_m = 0.25
        readcol, mock_file, mock_ra, mock_dec, mock_redshift, w_cp 
        mock_redshift = mock_redshift

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

    endif else if (strlowcase(catalog) EQ 'ldgdownnz') then begin 
        
        omega_m = 0.25
        readcol, mock_file, mock_ra, mock_dec, mock_redshift, w_cp 
        mock_redshift = mock_redshift

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

    endif else if (strlowcase(catalog) EQ 'qpm') then begin 
        omega_m = 0.31
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
    
    endif else if (strlowcase(catalog) EQ 'nseries') then begin 
    
        omega_m = 0.31
        GOTO, nseries_jump 
        
        upcp_index = where(w_cp gt 1) 
        ra_upcp  = mock_ra[upcp_indx]
        dec_upcp = mock_dec[upcp_indx]
        red_upcp = mock_redshift[upcp_indx]

        nocp_indx = where(w_cp eq 0)
        ra_nocp  = mock_ra[nocp_indx]
        dec_nocp = mock_dec[nocp_indx]
        red_nocp = mock_redshift[nocp_indx]
        
        spherematch, ra_nocp, dec_nocp, ra_upcp, dec_upcp, fib_angscale, m_nocp, m_upcp, d12, maxmatch=0

    endif else if (strlowcase(catalog) EQ 'patchy') then begin 
        omega_m = 0.31
        readcol, mock_file, mock_ra, mock_dec, mock_redshift, nbar, w_cp 

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

    endif else if (strmatch(catalog, '*bigmd*', /fold_case) EQ 1) then begin 

        omega_m = 0.31
        
        readcol, mock_file, mock_ra, mock_dec, mock_redshift, nbar, w_cp 

        upcp_indx = where(w_cp gt 1)   ; upweighted galaxies  
        ra_upcp  = mock_ra[upcp_indx]
        dec_upcp = mock_dec[upcp_indx]
        red_upcp = mock_redshift[upcp_indx]

        nocp_indx = where(w_cp eq 0)    ; galaxies with no weight 
        ra_nocp  = mock_ra[nocp_indx]
        dec_nocp = mock_dec[nocp_indx]
        red_nocp = mock_redshift[nocp_indx]
        
        ; run sppherematch 
        spherematch, ra_nocp, dec_nocp, ra_upcp, dec_upcp, fib_angscale, m_nocp, m_upcp, d12, maxmatch=0

    endif else if (strmatch(catalog, '*cmass*', /fold_case) EQ 1) then begin 
        ; CMASS, CMASS+LOWZ Combined Sample 
        ; omega_m = 0.274
        omega_m = 0.31
        ; ra, dec, z, nbar, wsys, wnoz, wfc, wcomp
        readcol, mock_file, mock_ra, mock_dec, mock_redshift, mock_nz, mock_wsys, mock_wnoz, mock_wfc, mock_comp

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

    print, 'Writing ... ', dlos_file 
    openw, lun, dlos_file, /get_lun

    all_dlos = [] 
    all_targ_ra  = []
    all_targ_dec = []
    all_targ_red = []
    all_neigh_ra  = [] 
    all_neigh_dec = [] 
    all_neigh_red = [] 

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

        d_los = dblarr(n_neigh)
        for j=0L,n_neigh-1L do begin
            all_dlos = [all_dlos, neigh_Dc[j] - targ_Dc]    ; neighbor - target
            ; target ra, dec, z
            all_targ_ra  = [all_targ_ra, targ_ra]
            all_targ_dec = [all_targ_dec, targ_dec]
            all_targ_red = [all_targ_red, targ_red]
            ; neighbor ra, dec, z 
            all_neigh_ra  = [all_neigh_ra, neigh_ra[j]]
            all_neigh_dec = [all_neigh_dec, neigh_dec[j]]
            all_neigh_red = [all_neigh_red, neigh_red[j]]
        endfor 
    endfor  

    if (strlowcase(catalog) EQ 'tilingmock') then begin 
        ; to account for double counting the fc pairs 
        print, '#dlos with repeats', n_elements(all_dlos)
        abs_dlos = abs(all_dlos)
        sorted_indx = sort(abs_dlos)
        sorted_abs_dlos = abs_dlos[sorted_indx]
        uniq_sorted_indx = uniq(sorted_abs_dlos)

        all_dlos = all_dlos[sorted_indx[uniq_sorted_indx]]  ; sorts the dlos but who cares
        print, '#dlos without repeats', n_elements(all_dlos)
    endif 

    for j=0L,n_elements(all_dlos)-1L do begin
        printf, lun, all_dlos[j], all_targ_ra[j], all_targ_dec[j], all_targ_red[j], $
            all_neigh_ra[j], all_neigh_dec[j], all_neigh_red[j], format='(f, f, f, f, f, f, f)'   ; d_LOS, target_ra, target_dec, target redshift, neighbor_ra, neighbor_dec, neighbor redshift
    endfor
    
    free_lun, lun 
    return 

    nseries_jump: orig_file = strmid(mock_file, 0, strpos(mock_file, 'fibcoll.dat'))+'rdzwc'
    readcol, orig_file , mock_ra, mock_dec, mock_redshift, wfkp, w_fc, z_upw, upw_index 

    noupw_index = where(w_fc eq 0, n_noupw) 
     
    ;; check z_upw with upw_index 
    ;for i=0L,n_noupw-1L do begin 
    ;    if (z_upw[noupw_index[i]] NE mock_redshift[upw_index[noupw_index[i]]]) then begin 
    ;        print, z_upw[noupw_index[i]], mock_redshift[upw_index[noupw_index[i]]]
    ;    endif 
    ;endfor 

    dlos = 3000.0*(comdis(mock_redshift[noupw_index], omega_m, 1.-omega_m)-$
        comdis(z_upw[noupw_index], omega_m, 1.-omega_m))

    print, dlos_file 
    openw, lun, dlos_file, /get_lun

    for j=0L,n_noupw-1L do begin
        i_targ = upw_index[noupw_index[j]]
        ; d_LOS, target_ra, target_dec, target redshift, neighbor_ra, neighbor_dec, neighbor redshift
        printf, lun, dlos[j], mock_ra[i_targ], mock_dec[i_targ], mock_redshift[i_targ], $
            mock_ra[j], mock_dec[j], mock_redshift[j], format='(f, f, f, f, f, f, f)'   
    endfor
    free_lun, lun 
return 
end 
