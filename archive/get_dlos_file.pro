function get_dlos_file, mock=mock, version=version, n_mock=n_mock
; Gets files once Mock Catalog, Version, and # of mock is specified
    if (strlowcase(mock) EQ 'pthalo') then begin 
        pthalo_mock_dir = '/mount/riachuelo1/hahn/data/PTHalo/'
        if ((version NE 'v11p0') AND (version NE 'v11p1') AND (version NE 'v11p2')) then begin 
            print, 'CANNOT FIND VERSION '+version
            stop
        endif else begin 
            pthalo_mock_dir = pthalo_mock_dir+version+'/'
            version_corr = strjoin(strsplit(version, 'p', /extract), '.') 
        endelse 
        pthalo_file_name = 'DLOS_cmass_dr11_north_ir4'+strmid(strtrim(string(n_mock+1000),1),1)+'.'+version_corr+'.wghtv.txt'

        mock_file_name = pthalo_mock_dir + pthalo_file_name 
    endif
    if (strlowcase(mock) EQ 'lasdamas') then begin
        if (version EQ 'geo') then begin 
            lasdamas_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
            lasdamas_n_mock = strmid(n_mock, 0, stregex(n_mock, '(a|b|c|d)'))
            lasdamas_letter = strmid(n_mock, stregex(n_mock, '(a|b|c|d)'))
            lasdamas_file = 'DLOS_sdssmock_gamma_lrgFull_zm_oriana'+strmid(strtrim(string(lasdamas_n_mock+100),1),1)+$
                lasdamas_letter+'_no.rdcz.dat'
            mock_file_name = lasdamas_dir+lasdamas_file
        endif  
    endif 
    if (strlowcase(mock) EQ 'tilingmock') then begin 
        tilingmock_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/' 
        tilingmock_file = 'DLOS_cmass-boss5003sector-icoll012.dat'
        mock_file_name = tilingmock_dir+tilingmock_file
    endif 
    return, mock_file_name
end 
