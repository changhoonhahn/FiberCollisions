pro build_mock_wcp_assign
; Builds a mock catalog using spherematch with fiber collision weights 
; Import mock catalog
    mock_dir     = '/mount/riachuelo1/hahn/data/tiling_mocks/'
    mock_file    = 'cmass-boss5003sector-icoll012.dat' 
    readcol, mock_dir+mock_file, mock_ra, mock_dec, mock_redshift, mock_w
    Ngal = n_elements(mock_ra)
    print, 'Number of mock galaxies = ', Ngal

; Fiber collision angular scales 
    fib_angscale = 0.01722
    spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, match1, match2, d12, maxmatch=0
    overlap_indx = where(d12 GT 0.0, overlap_count)
    match1 = match1[overlap_indx]
    match2 = match2[overlap_indx] 
    print, 'Number of non-overlapping matches = ', overlap_count

    gal1 = match1[uniq(match1, sort(match1))]
    gal2 = match2[uniq(match2, sort(match2))]

    wcp = replicate(1L,Ngal)
    for i=0L,n_elements(gal1)-1L do begin 
        if (wcp[gal1[i]] NE 0) then begin   ; If wcp[gal1[i]] = 0 then it's already been accounted for in a collision pair
            collision_indx = where(match1 EQ gal1[i]) 

; Index of target galaxy 
            targ_indx   = gal1[i]
; Index of neighboring galaxies
            neigh_indx  = match2[collision_indx]

            neigh_wcp = wcp[neigh_indx] 
            notzero = where( neigh_wcp NE 0, notzerocount)  ; So that we dont count galaxies that have already been counted
            wcp[targ_indx]  = wcp[targ_indx]+notzerocount
            wcp[neigh_indx] = 0L
        endif 
    endfor
    
    if (Ngal ne total(wcp)) then STOP 
    
    output_file = 'cmass-boss5003sector-icoll012.fibcoll.dat'
    openw, lun, mock_dir+output_file, /get_lun
    for i=0L,Ngal-1L do printf, lun, mock_ra[i], mock_dec[i], mock_redshift[i], wcp[i],format='(4F)'
    free_lun, lun 
end 
