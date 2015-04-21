pro get_mock_fibcoll_dlos
; Gets the line of sight displacements of fiber collided galaxy pairs 

; Read the mock catalog
    mock_dir     = '/mount/riachuelo1/hahn/data/tiling_mocks/'
    mock_file    = 'cmass-boss5003sector-icoll012.dat'
    readcol, mock_dir+mock_file, mock_ra, mock_dec, mock_redshift, mock_w
    print, mock_dir+mock_file
    Ngal = n_elements(mock_ra)
    print, Ngal, ' galaxies'

; Fiber collision angular sacle 
    fib_angscale = 0.01722
; Spherematch for fibercollision angular scale
    spherematch, mock_ra, mock_dec, mock_ra, mock_dec, fib_angscale, match1, match2, d12, maxmatch=0
    overlap_indx = where(d12 GT 0.0, overlapcount)
    match1 = match1[overlap_indx]
    match2 = match2[overlap_indx]
    print, 'Number of fiber collision matches = ', overlapcount
    
    iuniq1 = match1[uniq(match1,sort(match1))]
    iuniq2 = match2[uniq(match2,sort(match2))] 

    disp_los = []
    for i=0L,n_elements(iuniq1)-1L do begin
        collision_indx = where( match1 EQ iuniq1[i] )

        targ_indx = iuniq1[i]
        targ_phi = mock_ra[targ_indx]
        targ_theta = 90.0-mock_dec[targ_indx]
        targ_redshift = mock_redshift[targ_indx]
        targ_r = 3000.0*comdis(targ_redshift,0.27,0.73)

        neigh_indx = match2[collision_indx]
        neigh_phi = mock_ra[neigh_indx]
        neigh_theta = 90.0-mock_dec[neigh_indx]
        neigh_redshift = mock_redshift[neigh_indx]
        neigh_r = 3000.0*comdis(neigh_redshift,0.27,0.73)

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
    
    output_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'
    output_file = 'cmass-boss5003sector-icoll012.disp_los.dat' 
    openw, lun, output_dir+output_file, /get_lun
    for i=0L,n_elements(disp_los)-1L do printf, lun, disp_los[i], format='(f)'
    free_lun, lun 
return 
end 
