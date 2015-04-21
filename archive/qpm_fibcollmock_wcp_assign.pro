pro qpm_fibcollmock_wcp_assign, dr_version, n, ngc=ngc, sgc=sgc
; build QPM mocks with close pair fibercollision weights using spherematch with angular separation 
; and output dLOS for each fibercollided pairs
    if ((dr_version NE 'dr12') AND (dr_version NE 'dr11')) then begin
        print, "Data Release Version has to be 'dr11' or 'dr12'"
        stop 
    endif 
    if keyword_set(ngc) then ns_str = 'ngc'
    if keyword_set(sgc) then ns_str = 'sgc'
    ; read QPM mocks given Data Release version 
    mock_dir     = '/mount/riachuelo1/hahn/data/QPM/'+dr_version+'/'+ns_str+'/'
    mock_file    = 'a0.6452_'+strmid(strtrim(string(n+10000),1),1)+'.dr12c_'+ns_str+'.rdz'
    readcol, mock_dir+mock_file, mock_ra, mock_dec, mock_redshift           ; read QPM mocks
    Ngal = n_elements(mock_ra)
    print, Ngal, ' galaxies'
    ; use spherematch to find fiber collided pairs 
    fib_angscale = 0.01722          ; fiber collision angular scale 
    spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, m1, m2, d12, maxmatch=0
    overlap_index = d12 GT 0.0
    m1 = m1[where(overlap_index, overlap_index_count)] 
    m2 = m2[where(overlap_index)] 
    print, 'overlap_index_count= ', overlap_index_count
    gal1 = m1[uniq(m1, sort(m1))]
    gal2 = m2[uniq(m2, sort(m2))]
    wcp  = replicate(1.0, Ngal)
    print, 'unique galaxies in match', n_elements(gal1)
    ; output dLOS file
    openw, lun, mock_dir+'DLOS_'+mock_file+'.dat', /get_lun
    dlos = []
    for i=0L,n_elements(gal1)-1L do begin 
        if (wcp[gal1[i]] GT 0.0) then begin 
            collision_index = where(m1 EQ gal1[i]) 
            targ_index   = gal1[i]
            targ_phi    = mock_ra[gal1[i]]
            targ_theta  = 90.0-mock_dec[gal1[i]]
            targ_red    = mock_redshift[gal1[i]]
            targ_r      = 3000.0*comdis(targ_red,0.27,0.73)

            neigh_index  = m2[collision_index]
            neigh_wcp   = wcp[neigh_index]
            notzero = where(neigh_wcp GT 0.0, notzerocount)
            if (notzerocount GT 0) then begin 
                wcp[targ_index] = wcp[targ_index]+float(notzerocount)
                wcp[neigh_index[notzero]] = wcp[neigh_index[notzero]]-1.0
                neigh_phi   = mock_ra[neigh_index[notzero]]
                neigh_theta = 90.0-mock_dec[neigh_index[notzero]]
                neigh_red   = mock_redshift[neigh_index[notzero]]
                neigh_r     = 3000.0*comdis(neigh_red,0.27,0.73)

                angles_to_xyz, targ_r, targ_phi, targ_theta, targ_x, targ_y, targ_z
                angles_to_xyz, neigh_r, neigh_phi, neigh_theta, neigh_x, neigh_y, neigh_z
            ; calculate the target magnitude  
                targ_mag= targ_x^2+targ_y^2+targ_z^2
                targ_dot_neigh= targ_x*neigh_x+targ_y*neigh_y+targ_z*neigh_z
                proj = targ_dot_neigh/targ_mag
            ; calculate the projected x,y,z 
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
                dlos = [dlos, d_los]
            endif 
        endif 
    endfor
    print, 'dLOS n_elements =', n_elements(dlos)
    print, 'dLOS mean =', mean(dlos)
    print, 'dLOS median =', median(dlos)
    print, 'dLOS stddev =', stddev(dlos)
    print, '(# of dLOS)/(total # of galaxies)', float(n_elements(dlos))/float(Ngal)
    free_lun, lun 
    print, 'Total number of galaxies = ', Ngal 
    print, 'Total close pair weights = ', total(wcp)
    if (Ngal ne total(wcp)) then STOP 
    ; Write mocks with fibercollision weights
    fibcol_mock_file = mock_file+'.fibcoll.dat'
    openw, lun, mock_dir+fibcol_mock_file, /get_lun
    for i=0L,Ngal-1L do begin 
        printf, lun, mock_ra[i], mock_dec[i], mock_redshift[i], wcp[i], format='(4F)'
    endfor
    free_lun, lun 
return 
end 
