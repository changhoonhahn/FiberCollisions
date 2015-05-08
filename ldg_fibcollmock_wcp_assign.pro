pro ldg_fibcollmock_wcp_assign, n, letter 
; build LasDamasGeo mock catalog with close pair fibercollision weights
; using spherematch with angular separation 
; and output dLOS for each fibercollided pairs
; read LasDamas Geo mock
    true_file = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+strmid(strtrim(string(n+100),1),1)+letter+'_no.rdcz.dat'
    readcol, true_file, mock_ra, mock_dec, mock_redshift
    Ngal = n_elements(mock_ra)
    print, Ngal, ' galaxies'

    ; Find fiber collided pairs 
    fib_angscale = 0.01722
    spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, m1, m2, d12, maxmatch=0
    overlap_index = d12 GT 0.0
    m1 = m1[where(overlap_index, overlap_index_count)] 
    m2 = m2[where(overlap_index)] 

    print, 'overlap_index_count= ', overlap_index_count
    gal1 = m1[uniq(m1, sort(m1))]
    gal2 = m2[uniq(m2, sort(m2))]
    wcp  = replicate(1.0, Ngal)
    print, 'unique galaxies in match', n_elements(gal1)

    for i=0L,n_elements(gal1)-1L do begin 
        if (wcp[gal1[i]] GT 0.0) then begin 
            collision_index = where(m1 EQ gal1[i]) 
            targ_index   = gal1[i]

            neigh_index  = m2[collision_index]
            neigh_wcp   = wcp[neigh_index]
            notzero = where(neigh_wcp GT 0.0, notzerocount)
            if (notzerocount GT 0) then begin 
                wcp[targ_index] = wcp[targ_index]+float(notzerocount)
                wcp[neigh_index[notzero]] = wcp[neigh_index[notzero]]-1.0
            endif 
        endif 
    endfor
    if (Ngal ne total(wcp)) then STOP 
    
    ; write mocks with fibercollision weights
    openw, lun, '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+strmid(strtrim(string(n+100),1),1)+$
        letter+'_no.rdcz.fibcoll.dat', /get_lun

    for i=0L,Ngal-1L do begin 
        printf, lun, mock_ra[i], mock_dec[i], mock_redshift[i], wcp[i],format='(4F)'
    endfor
    free_lun, lun 
end 
