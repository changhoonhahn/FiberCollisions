pro build_wcp_assign, catalog, catalog_param=catalog_parm 
; Assign fiber collision weights for mock catalogs taht do not have them 
    if strtrim(catalog,2) eq 'tilingmock' then begin 
        ; import tiling mock catalog
        mock_dir     = '/mount/riachuelo1/hahn/data/tiling_mocks/'
        mock_file    = 'cmass-boss5003sector-icoll012.zlim.dat' 
        readcol, mock_dir+mock_file, mock_ra, mock_dec, mock_redshift, mock_w
        Ngal = n_elements(mock_ra)
        print, Ngal, ' galaxies'

        ; Fiber collision angular scales 
        fib_angscale = 0.01722
        spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, match1, match2, d12, maxmatch=0
        overlap_indx = where(d12 GT 0.0, overlap_count)
        match1 = match1[overlap_indx]
        match2 = match2[overlap_indx] 

        gal1 = match1[uniq(match1, sort(match1))]
        gal2 = match2[uniq(match2, sort(match2))]

        wcp = replicate(1.0, Ngal)
        for i=0L,n_elements(gal1)-1L do begin 
            if (wcp[gal1[i]] NE 0) then begin   ; If wcp[gal1[i]] = 0 then it's already been accounted for in a collision pair
                collision_indx = where(match1 EQ gal1[i]) 

                ; index of target galaxy 
                targ_indx   = gal1[i]
                ; index of neighboring galaxies
                neigh_indx  = match2[collision_indx]

                neigh_wcp = wcp[neigh_indx] 
                notzero = where(neigh_wcp GT 0.0, notzerocount)  ; So that we dont count galaxies that have already been counted
                wcp[targ_indx]  = wcp[targ_indx]+notzerocount
                if (notzerocount GT 0) then $ 
                    wcp[neigh_indx[notzero]] = wcp[neigh_indx[notzero]]-replicate(1.0, notzerocount) 
            endif 
        endfor
        
        if (Ngal ne total(wcp)) then begin
            print, 'Ngal and total(wcp) does not match'
            print, 'Ngal = ', Ngal
            print, 'Total(wcp) = ', total(wcp) 
            STOP 
        endif 
        
        output_file = 'cmass-boss5003sector-icoll012.zlim.fibcoll.dat'
        openw, lun, mock_dir+output_file, /get_lun
        for i=0L,Ngal-1L do printf, lun, mock_ra[i], mock_dec[i], mock_redshift[i], wcp[i], format='(4F)'
        free_lun, lun 
        return 
    endif else if strtrim(catalog,2) eq 'lasdamasgeo' then begin 
        print, 'not yet coded'
        return 
    endif 
end 
