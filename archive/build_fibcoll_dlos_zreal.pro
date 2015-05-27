pro build_fibcoll_dlos_zreal, catalog, mock_file, dlos_file
; build line-of-sight separation distances for fibercollided pairs in Mock Catalogs with REAL redshift WITHOUT RSD
; HARDCODED FOR QPM ONLY
; spherematch is used to compare upweighted galaxies with unweighted galaxies
    fib_angscale = 0.01722
    print, mock_file     
    if (strlowcase(catalog) EQ 'qpm') then begin 
        omega_m = 0.31
        omega_l = 0.69
        ; read in mock file 
        readcol,mock_file, mock_ra, mock_dec, mock_redshift_rsd, w_fkp, w_cp

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
        red_upcp = mock_redshift_real[upcp_indx]
        red_upcp_rsd = mock_redshift_rsd[upcp_indx]

        nocp_indx = where(w_cp eq 0)
        ra_nocp  = mock_ra[nocp_indx]
        dec_nocp = mock_dec[nocp_indx]
        red_nocp = mock_redshift_real[nocp_indx]
        red_nocp_rsd = mock_redshift_rsd[nocp_indx]

        spherematch, ra_nocp, dec_nocp, ra_upcp, dec_upcp, fib_angscale, m_nocp, m_upcp, d12, maxmatch=0
    endif else begin 
        ; HARDCODED FOR QPM ONLY
        print, 'only for QPM'
        stop 
    endelse

    gal_upcp = m_upcp[uniq(m_upcp, sort(m_upcp))]
    ;print,'gal_upcp=',n_elements(gal_upcp)
    ;print,'uniq(gal_upcp)=',n_elements(uniq(gal_upcp[sort(gal_upcp)]))
    ;print,'m_upcp=',n_elements(m_upcp)
    ;print,'uniq(m_upcp)=',n_elements(uniq(m_upcp[sort(m_upcp)]))
    ;print,'m_nocp=',n_elements(m_nocp)
    ;print,'uniq(m_nocp)=',n_elements(uniq(m_nocp[sort(m_nocp)]))
    ;print,'total(w_cp)=',total(w_cp)
    ;print,'total(w_upcp)=',total(w_cp[upcp_indx[m_upcp[gal_upcp]]])
    ;print,'total(w_nocp)=',total(w_cp[nocp_indx[uniq(m_nocp[sort(m_nocp)])]])

    print, dlos_file 
    openw, lun, dlos_file, /get_lun

    for i=0L,n_elements(gal_upcp)-1L do begin
        ngal=gal_upcp[i]
        collision_indx = where( m_upcp eq ngal )

        targ_red = red_upcp[ngal]
        targ_red_rsd = red_upcp_rsd[ngal]
        targ_Dc = 3000.0*comdis(targ_red, omega_m, omega_l)
        targ_Dc_rsd = 3000.0*comdis(targ_red_rsd, omega_m, omega_l)

        neigh_red = red_nocp[m_nocp[collision_indx]]
        neigh_red_rsd = red_nocp_rsd[m_nocp[collision_indx]]
        neigh_Dc = 3000.0*comdis(neigh_red, omega_m, omega_l)
        neigh_Dc_rsd = 3000.0*comdis(neigh_red_rsd, omega_m, omega_l)

        n_neigh = n_elements(neigh_red)

        d_los = dblarr(n_neigh)
        d_los_rsd = dblarr(n_neigh)
        for j=0L,n_neigh-1L do begin
            d_los[j] = neigh_Dc[j] - targ_Dc
            d_los_rsd[j] = neigh_Dc_rsd[j] - targ_Dc_rsd

            ;print, d_los[j], d_los_rsd[j]
            printf, lun, d_los[j], format='(f)'
        endfor
    endfor  
    free_lun, lun 
return 
end 
