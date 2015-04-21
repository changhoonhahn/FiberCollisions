pro mock_lasdamasbox_shift, n
; Detects close pairs in x-y plane of LasDamas Box simulation and then 
; imposes fibercollision weights
    mock_dir = '/mount/riachuelo1/hahn/data/las_damas/box/'
    mock_file = 'mock_gamma_lrg21p2_zmo_Oriana_1'+$
        strmid(strtrim(string(n+1000),2),1)+'_z0p342_fof_b0p2.zdist.dat'
; import mock x,y,z data
    readcol, mock_dir+mock_file, mock_x, mock_y, mock_z
    ngal = n_elements(mock_x) 
    print, 'Number of galaxies =', ngal
    mock_w = replicate(1.0, ngal)
; use matchnd.pro to detect close pairs
    matchxy = dblarr(2,n_elements(mock_x))
    matchxy[0,*] = mock_x
    matchxy[1,*] = mock_y
; hardcoded match distance 
    matchdist = 0.2
    matchnd, matchxy, matchxy, matchdist, m1=m1, m2=m2, d12=d12, maxmatch=0 
    print, 'Number of total matches =', n_elements(m1)
    nooverlap = d12 GT 0.0d
    m1 = m1[where(nooverlap, nooverlap_count)]
    m2 = m2[where(nooverlap)]
    print, 'Number of matches that are not overlaped =', nooverlap_count
    gal1 = m1[uniq(m1, sort(m1))]
    gal2 = m2[uniq(m2, sort(m2))]
    dlos = [] 
    print, 'Number of matches with unique galaxies =', n_elements(gal1)
; open dLOS file
    data_dir = '/mount/riachuelo1/hahn/data/las_damas/box/'
    dlos_file = 'mock_gamma_lrg21p2_zmo_Oriana_1'+$
        strmid(strtrim(string(n+1000),2),1)+'_z0p342_fof_b0p2.zdist.dlos.dat'
    openw, lun, data_dir+dlos_file, /get_lun
    for i=0L, n_elements(gal1)-1L do begin 
        if (mock_w[gal1[i]] GT 0.0) then begin 
            fibcol_index = where(m1 EQ gal1[i])
            targ_index = gal1[i]
            neigh_index = m2[fibcol_index]
            
            neigh_wcp = mock_w[neigh_index]
            notzero = where(neigh_wcp GT 0.0, notzerocount)
            if (notzerocount GT 0) then begin 
                mock_w[targ_index] = mock_w[targ_index]+float(notzerocount)
                mock_w[neigh_index[notzero]] = mock_w[neigh_index[notzero]]-1.0
                dlos = [dlos, mock_z[neigh_index]-mock_z[targ_index]]
                for i_neigh=0L,notzerocount-1L do begin 
                ; write to dLOS file
                    printf, lun, mock_z[neigh_index[notzero[i_neigh]]]-mock_z[targ_index], format='(f)'
                endfor 
            endif
        endif 
    endfor
    print, 'dLOS n_elements =', n_elements(dlos)
    print, 'dLOS mean =', mean(dlos)
    print, 'dLOS median =', median(dlos)
    print, 'dLOS stddev =', stddev(dlos)
    print, '(# of dLOS)/Ngal = ', float(n_elements(dlos))/float(ngal) 
; close dLOS file
    print, 'number of close pairs=', n_elements(mock_w[where(mock_w gt 1.0)])
    free_lun, lun
    print, 'total weight = ', total(mock_w) ; should be same as Ngal 
; output data file with fibcol weight
    fibcol_file = 'mock_gamma_lrg21p2_zmo_Oriana_1'+$
        strmid(strtrim(string(n+1000),2),1)+'_z0p342_fof_b0p2.zdist.fibcol.dat'
    openw, lun, data_dir+fibcol_file, /get_lun
    for i=0L,n_elements(mock_x)-1L do begin 
        printf, lun, mock_x[i], mock_y[i], mock_z[i], mock_w[i], format='(f,f,f,f)'
    endfor 
    free_lun, lun
return
end 
