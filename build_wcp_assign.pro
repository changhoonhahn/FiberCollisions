pro build_wcp_assign, catalog, catalog_param=catalog_param 
; Assign fiber collision weights for mock catalogs taht do not have them 
    fib_angscale = 0.01722
    if strtrim(catalog,2) eq 'tilingmock' then begin 
        ; import tiling mock catalog
        mock_dir     = '/mount/riachuelo1/hahn/data/tiling_mocks/'
        mock_file    = 'cmass-boss5003sector-icoll012.zlim.dat' 
        readcol, mock_dir+mock_file, mock_ra, mock_dec, mock_redshift, mock_w
        Ngal = n_elements(mock_ra)
        print, Ngal, ' galaxies'

        ; Fiber collision angular scales 
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
    endif else if strtrim(catalog, 2) EQ 'nseries' then begin ; N series mocks 
        data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
        orig_file = 'CutskyN'+string(catalog_param)+'.rdzwc'
        info_file = 'CutskyN'+string(catalog_param)+'.mask_info'

        readcol, data_dir+orig_file, mock_ra, mock_dec, mock_z, dum, wfc, z_upw
        readcol, data_dir+info_file, mock_comp
        ngal = n_elements(mock_ra)

        no_fib = where(wfc eq 0, n_no_fib)      ; galaxies without fibers
        fib = where(wfc ge 1, n_fib)            ; galaxies with fibers

        spherematch, mock_ra[no_fib], mock_dec[no_fib], mock_ra[fib], mock_dec[fib], fib_angscale, $
            match_nofib, match_fib, d12, maxmatch=0

        gal_nofib = match_nofib[uniq(match_nofib, sort(match_nofib))]       ; unique no-fiber galaxies
    
        upw_ra = replicate(0.0, ngal)
        upw_dec = replicate(0.0, ngal)
        upw_z = replicate(0.0, ngal)

        for i=0L,n_elements(gal_nofib)-1L do begin 
            collision_index = where(match_nofib eq gal_nofib[i], n_coll)
            if n_coll gt 1 then begin
                print, n_coll, ' should be 1'
                print, wfc[fib[match_fib[collision_index]]]
                stop 
            endif 

            upw_ra[no_fib[match_nofib[collision_index]]] = mock_ra[fib[match_fib[collision_index]]]
            upw_dec[no_fib[match_nofib[collision_index]]] = mock_dec[fib[match_fib[collision_index]]]
            upw_z[no_fib[match_nofib[collision_index]]] = mock_z[fib[match_fib[collision_index]]]
        endfor 

        output_file  = 'CutskyN'+string(catalog_param)+'.rdzwc.dat'
        openw, lun, data_dir+output_file, /get_lun
        for i=0L,Ngal-1L do printf, lun, mock_ra[i], mock_dec[i], mock_z[i], wfc[i], z_upw[i], $
            upw_ra[i], upw_dec[i], upw_z[i], format='(8F)'
        free_lun, lun 
        return 
    endif else if strtrim(catalog,2) eq 'lasdamasgeo' then begin 
        print, 'not yet coded'
        return 
    endif 
end 
