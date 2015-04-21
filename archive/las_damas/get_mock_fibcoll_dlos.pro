pro get_mock_fibcoll_dlos,n,letter
    mock_dir     = '/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/'
    mock_file    = 'sdssmock_gamma_lrgFull_zm_oriana'+strmid(strtrim(string(n+100),1),1)+$
        letter+'_no.rdcz.dat'

    readcol, mock_dir+mock_file, mock_ra, mock_dec, mock_redshift
    Ngal = n_elements(mock_ra)
    print, Ngal, ' galaxies'
    mock_redshift = mock_redshift/299800.0

    fib_angscale = 0.01722
    spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, m1, m2, d12, maxmatch=0
    overlap_indx = d12 GT 0.0
    m1 = m1[where(overlap_indx, overlap_indx_count)]
    m2 = m2[where(overlap_indx)]
    print, 'Number of fiber collided matches =', overlap_indx_count

    gal1 = m1[uniq(m1,sort(m1))]
    gal2 = m2[uniq(m2,sort(m2))]

    disp_los = []
    for i=0L,n_elements(gal1)-1L do begin 
        collision_indx = where(m1 EQ gal1[i])

        targ_indx   = gal1[i]
        targ_phi    = mock_ra[gal1[i]] 
        targ_theta  = 90.0-mock_dec[gal1[i]]
        targ_red    = mock_redshift[gal1[i]]
        targ_r      = 3000.0*comdis(targ_red,0.27,0.73)

        neigh_indx  = m2[collision_indx]
        neigh_phi   = mock_ra[m2[collision_indx]]
        neigh_theta = 90.0-mock_dec[m2[collision_indx]]
        neigh_red   = mock_redshift[m2[collision_indx]]
        neigh_r     = 3000.0*comdis(neigh_red,0.27,0.73)

        angles_to_xyz, targ_r, targ_phi, targ_theta, targ_x, targ_y, targ_z
        angles_to_xyz, neigh_r, neigh_phi, neigh_theta, neigh_x, neigh_y, neigh_z

        targ_mag= targ_x^2+targ_y^2+targ_z^2
        targ_dot_neigh= targ_x*neigh_x+targ_y*neigh_y+targ_z*neigh_z
        proj = targ_dot_neigh/targ_mag

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
        endfor
        disp_los = [disp_los, d_los]
    endfor
    
    datadir = '/mount/riachuelo1/hahn/data/las_damas/fibcoll/'
    
    openw, lun, datadir+'sdssmock_gamma_lrgFull_zm_oriana'+strmid(strtrim(string(n+100),1),1)+letter+$
        '_no.rdcz.disp_los.dat', /get_lun
    for i=0L,n_elements(disp_los)-1L do printf, lun, disp_los[i], format='(f)'
    free_lun, lun 
end 
