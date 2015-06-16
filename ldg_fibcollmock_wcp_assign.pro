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
    ;overlap_index = d12 GT 0.0
    ;m1 = m1[where(overlap_index, overlap_index_count)] 
    ;m2 = m2[where(overlap_index)] 
    ;print, 'overlap_index_count= ', overlap_index_count
    gal1 = m1[uniq(m1, sort(m1))]

    wcp  = replicate(1.0, Ngal)         ; close pair weights (fiber collision weights) 
    z_upw = replicate(-1, Ngal)         ; redshift of upweighted galaxy 
    upw_index = replicate(-1, Ngal)     ; index of upweighted galaxy 
    print, 'unique galaxies in match', n_elements(gal1)

    for i=0L,n_elements(gal1)-1L do begin 

        targ_wcp = wcp[gal1[i]]

        if (targ_wcp GT 0.0) then begin 

            collision_index = where(m1 EQ gal1[i]) 
            targ_index   = gal1[i]

            neigh_index  = m2[collision_index]
            neigh_wcp   = wcp[neigh_index]

            notzero = where((neigh_wcp GT 0.0) AND (neigh_index NE targ_index), notzerocount)

            neigh_max_wcp = max(wcp[neigh_index[notzero]], neigh_max_index)   ; maximum wcp weight of neighbors

            if (notzerocount GT 0) then begin 

                ;print, 'B ', wcp[targ_index], '---', wcp[neigh_index[notzero]]
                if (targ_wcp GE neigh_max_wcp) then begin  
                    ; add up galaxy weights of galaxies that have weights
                    wcp[targ_index] = wcp[targ_index]+float(notzerocount)
                    
                    ; reduce weight of fiber collided galaxies 
                    wcp[neigh_index[notzero]] = wcp[neigh_index[notzero]]-1.0

                    ; redshift of upweighted galaxy
                    z_upw[neigh_index[notzero]] = mock_redshift[targ_index]
                    
                    ; index of upweighted galaxy 
                    upw_index[neigh_index[notzero]] = targ_index

                endif else begin 
                    ; Upweight the already upweighted galaxies (this is in the case there are 
                    ; triplets or 2 pairs close to each other 
                    
                    print, targ_wcp, neigh_max_wcp

                    wcp[targ_index] = wcp[targ_index] - 1.0 

                    wcp[neigh_index[notzero[neigh_max_index]]] = wcp[neigh_index[notzero[neigh_max_index]]] + 1.0

                    z_upw[targ_index] = mock_redshift[neigh_index[notzero[neigh_max_index]]]
                    upw_index[targ_index] = neigh_index[notzero[neigh_max_index]]
                endelse 
                ;print, 'A ', wcp[targ_index], '---', wcp[neigh_index[notzero]]
            endif 
        endif 
    endfor

    print, float(n_elements(wcp[where(wcp GT 1.0)]))/float(n_elements(wcp))

    if (Ngal ne total(wcp)) then STOP 
    
    ; write mocks with fibercollision weights
    openw, lun, '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+strmid(strtrim(string(n+100),1),1)+$
        letter+'_no.rdcz.fibcoll.dat', /get_lun

    for i=0L,Ngal-1L do begin 
        printf, lun, mock_ra[i], mock_dec[i], mock_redshift[i], wcp[i], z_upw[i], upw_index[i], format='(5F, I)'
    endfor
    free_lun, lun 
end 
