pro build_pthalo_fibcoll_dlos 
; Simple hardcoded code to build a series of dLOS values for PTHalo mocks 
; perhaps later convert to other mock catalogs
    for i_mock=1L,10L do begin 
        build_fibcoll_dlos, mock='pthalo', version='v11p0', n_mock=i_mock 
    endfor 
return 
end 
