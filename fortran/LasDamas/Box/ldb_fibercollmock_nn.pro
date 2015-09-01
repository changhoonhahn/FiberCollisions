function sample_cdf, yrand 
; obtains x value from input CDF and y value 
; using crude interpolation 
    common d_losvar, dlos_cdf, dlos_val 
    if (n_elements(cdf) NE n_elements(xval)) then STOP 
    return, interpol(dlos_cdf, dlos_val, yrand)  
end 

pro ldb_fibercollmock_nn, n
; LasDamasBox Fibercollided Mock Catalog Generator using Nearest Neighbor upweight
; Reads in the x,y,z coordinates for the LasDamas Box mocks
; and imposes fiber collisions. 
; upweights nearest neighbor after shifting galaxy by dLOS 
; value obtained from dLOS sampling.
    common d_losvar, dlos_cdf, dlos_val
    t00 = systime(1)
; import mock LasDamas box x,y,z file
    mock_dir  = '/mount/riachuelo1/hahn/data/las_damas/box/'
    mock_file = 'mock_gamma_lrg21p2_zmo_Oriana_1'+$
        strmid(strtrim(string(n+1000),2),1)+'_z0p342_fof_b0p2.zdist.dat'
    readcol, mock_dir+mock_file, mock_x, mock_y, mock_z 
    ngal = n_elements(mock_x) & print, ngal
    mock_w = replicate(1.0,ngal)  
; sample from dLOS distribution of LasDamas mocks to impose 
; fiber collision issues
    ; import dLOS distribution; currently hardcoded to LasDamasGeo 1a
    lasdamas_dir = '/mount/riachuelo1/hahn/data/las_damas/fibcoll/'
    letter = 'a'
    dlos_file = 'sdssmock_gamma_lrgFull_zm_oriana'+$
        strmid(strtrim(string(1+100),2),1)+letter+'_no.rdcz.disp_los.dat'
    readcol, lasdamas_dir+dlos_file, dlos 
    n_dlos = n_elements(dlos)  
    print, 'mean dlos', mean(dlos), n_dlos

    dlos_pdf = float(histo(dlos, dlos_val, binsize=0.1))/float(n_dlos)  ; pdf from histogram
    pdf_notzero = where(dlos_pdf GT 0.0, n_pdf_notzero) 
    dlos_pdf = dlos_pdf[pdf_notzero]
    dlos_val = dlos_val[pdf_notzero] 
    dlos_cdf = total(dlos_pdf, /cumulative) ;cdf from histogram 
 
    fibcol_fract = float(n_dlos)/float(113129)  ; ngal from LasDamasGeo 1a
    n_fibcol = round(float(ngal)*fibcol_fract) & print, n_fibcol
    sample = randomu(rseed, n_fibcol) 
    dlos_sample = fltarr(n_fibcol) 
    for i=0L,n_fibcol-1L do dlos_sample[i] = sample_cdf(sample[i])  
    print, mean(dlos_sample)    ; should be ~0.0 gets ~0.5
; fiber collided galaxies 
    fibcol_gal_index = cgrandomindices(ngal,n_fibcol)   ; randomly chosen galaxy indices
    fibcol_gal = fltarr(3,n_fibcol)
    fibcol_gal[0,*] = mock_x[fibcol_gal_index]
    fibcol_gal[1,*] = mock_y[fibcol_gal_index]
    fibcol_gal[2,*] = mock_z[fibcol_gal_index]+dlos_sample
; not fiber collided galaxies
    notfibcol_gal_index = cmset_op(lindgen(ngal),'and',/not2,fibcol_gal_index)
    all_gal = fltarr(3,n_elements(notfibcol_gal_index))
    all_gal[0,*] = mock_x[notfibcol_gal_index]
    all_gal[1,*] = mock_y[notfibcol_gal_index]
    all_gal[2,*] = mock_z[notfibcol_gal_index]
; find nearest non-fibercollided galaxies nearest to the fiber collided galaxies 
    matchnd, all_gal, fibcol_gal, 25.0, m1=m1, m2=m2, d12=d12, maxmatch=0
    isort   = sort(m2)
    iuniq   = uniq(m2[isort])
    if (n_elements(iuniq) LT n_fibcol) then begin ; ensures that every fibcol galaxy has a match 
        no_neighbor = cmset_op(lindgen(n_fibcol),'and',/not2,(m2[isort])[iuniq]) 
        print, n_elements(no_neighbor), ' galaxies do not have a neighbor within 25 Mpc/h'
    endif 
    istart = 0L
    tot_min_d = 0.0
    cp_dlos = []
    for i=0L,n_elements(iuniq)-1L do begin 
        iend = iuniq[i]
        icurr = isort[istart:iend] 
        min_d = min(d12[icurr],min_index)
        tot_min_d = tot_min_d+min_d

        upindex = m1[icurr[min_index]]
        fibindex = m2[icurr[0]] 

        targ_gal = fibcol_gal[*,fibindex]
        min_dlos = all_gal[2,upindex]-targ_gal[2]

        mock_w[fibcol_gal_index[fibindex]] = mock_w[fibcol_gal_index[fibindex]]-1.0
        mock_w[notfibcol_gal_index[upindex]] = mock_w[notfibcol_gal_index[upindex]]+1.0
        cp_dlos  = [cp_dlos, min_dlos]
        istart = iend+1L 
    endfor
    scales = [50.0, 100.0, 200.0, 300.0, 400.0, 500.0] 
    if (n_elements(no_neighbor) GT 0) then begin 
        for i=0L,n_elements(no_neighbor)-1L do begin 
            iscales = 0L
            n_near = 0L  
            while (n_near EQ 0L) do begin 
                near = where((all_gal[0,*] GT fibcol_gal[0,no_neighbor[i]]-scales[iscales]) AND $
                    (all_gal[0,*] LT fibcol_gal[0,no_neighbor[i]]+scales[iscales]) AND $
                    (all_gal[1,*] GT fibcol_gal[1,no_neighbor[i]]-scales[iscales]) AND $
                    (all_gal[1,*] LT fibcol_gal[1,no_neighbor[i]]+scales[iscales]) AND $
                    (all_gal[2,*] GT fibcol_gal[2,no_neighbor[i]]-scales[iscales]) AND $
                    (all_gal[2,*] LT fibcol_gal[2,no_neighbor[i]]+scales[iscales]), n_near) 
                iscales = iscales+1L 
            endwhile 
            neighbors = all_gal[*,near]
            targ_gal = fibcol_gal[*,no_neighbor[i]]
            d_neighbors = sqrt((neighbors[0,*]-targ_gal[0])^2+(neighbors[1,*]-targ_gal[1])^2+$
                (neighbors[2,*]-targ_gal[2])^2) 
            min_dist = min(d_neighbors, min_d_index) 
            min_dlos = neighbors[2, min_d_index]-targ_gal[2]

            mock_w[fibcol_gal_index[no_neighbor[i]]] = mock_w[fibcol_gal_index[no_neighbor[i]]]-1.0
            mock_w[notfibcol_gal_index[near[min_d_index]]] = mock_w[notfibcol_gal_index[near[min_d_index]]]+1.0
            cp_dlos  = [ cp_dlos, min_dlos]
        endfor 
    endif 
    print, 'average minimum distance', tot_min_d/float(n_elements(iuniq))
    print, total(mock_w[fibcol_gal_index])          ; should also be 0 
    print, total(abs(mock_w[fibcol_gal_index]))     ; weight sanity checks, should be 0
    print, total(mock_w), ngal, long(total(mock_w)) EQ ngal       ; make sure the weights add up right 

    dlos_output_file = 'mock_gamma_lrg21p2_zmo_Oriana_1'+$
        strmid(strtrim(string(n+1000),2),1)+'_z0p342_fof_b0p2.zdist.dlos.dat'
    openw, lun, mock_dir+dlos_output_file, /get_lun 
    for i=0L,n_elements(cp_dlos)-1L do begin 
        printf, lun, cp_dlos[i]
    endfor 
    free_lun, lun 

    output_file = 'mock_gamma_lrg21p2_zmo_Oriana_1'+$
        strmid(strtrim(string(n+1000),2),1)+'_z0p342_fof_b0p2.zdist.fibcoll.dat'
    openw, lun, mock_dir+output_file, /get_lun 
    for i=0L,ngal-1L do begin 
        printf, lun, mock_x[i], mock_y[i], mock_z[i], mock_w[i]
    endfor 
    free_lun, lun 
    print, 'Total Time = ', (systime(1)-t00)/60.0 
return
end
