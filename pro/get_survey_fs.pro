pro get_survey_fs, mock_file
    ; Calculate f_s of a survey by 
    
    fib_angscale = 0.01722
    
    readcol, mock_file, ra, dec, z, w_cp, comp 
    

    spherematch, ra, dec, ra, dec, fib_angscale, m1, m2, d12, maxmatch=0

    no_overlap = where(d12 GT 0.0)
    d12 = d12[no_overlap]
    print, n_elements(w_cp[where(w_cp eq 0)])
    print, n_elements(d12)
   
    nocoll = where(w_cp eq 1)
    spherematch, ra[nocoll], dec[nocoll], ra[nocoll], dec[nocoll], fib_angscale, m1, m2, d12, maxmatch=0

    no_overlap = where(d12 GT 0.0)
    d12 = d12[no_overlap]
    print, n_elements(d12)/2

return 
end 
