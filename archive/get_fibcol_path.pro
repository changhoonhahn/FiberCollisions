function get_fibcol_path, ext, data=data, fftee=fftee, power=power, $
    cmass=cmass, tiling=tiling, lasdamas=lasdamas, pthalo=pthalo, $
    ria=ria, chi=chi 
; function to manage the paths for fiber collision codes
    if keyword_set(ria) then path = '/mount/riachuelo1/hahn/'
    if keyword_set(chi) then path = 'chichichihchi'     ; edit later
    if ((keyword_set(ria) EQ 0) AND (keyword_set(chi) EQ 0)) then $
        path = '/mount/riachuelo1/hahn/'
; data, FFT, or power? 
    if keyword_set(data) then path = path+'data/'$
        else if keyword_set(fftee) then path = path+'FFT/'$
        else if keyword_set(power) then path = path+'power/'$
        else STOP  
; using which data/simulation? 
    if keyword_set(cmass) then path = path 
    if keyword_set(tiling) then path = path+'tiling_mocks/'
    if keyword_set(lasdamas) then path = path+'las_damas/'
    if keyword_set(pthalo) then path = path+'manera_mock/'
; for extra extensions 
    if keyword_set(ext) then path = path+ext
    return, path 
end 
