pro plot_power_cmass_shotnoise_comparison, drversion, NS
; CMASS P(k) comparison for different shot noise methods
    ;shotnoises = ['stdfkp', 'our', 'florian', 'cole']
    shotnoises = ['our', 'florian', 'cole', 'real florian']
    n_sn = n_elements(shotnoises)
    ;sn_colors = ['black', 'blue', 'green', 'red']
    sn_colors = ['blue', 'green', 'red', 'black']
    sn_labels = strarr(n_sn)
    for i=0L,n_sn-1L do sn_labels[i] = 'CMASS P(k) '+shotnoises[i]+' shotnoise'

    ; configure plot
    fig_dir = '/home/users/hahn/research/figures/boss/powerspectrum/'
    fig_file = fig_dir+'power-cmass-'+NS+'-'+strjoin(shotnoises,'-')+$
        '-shotnoise-3600lbox-960grid-480bin.eps'
    im_plotconfig, 0, pos, psfile=fig_file, charsize=1.8

    xrange=[10.^(-3.), 10.^0.]
    yrange=[10.^3., 10.^5.5]
    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle=textoidl('k'), $
        ytitle=textoidl('P(k)'), /xlog, /ylog
    im_legend, sn_labels, linestyle=replicate(0,n_sn), /bottom, /left, box=0, $
        color=sn_colors, psym=replicate(16,n_sn), $
        symsize=replicate(1.0,n_sn)*1.7, symthick=8, spacing=2.7, charsize=1.7
    if (NS EQ 'N') then im_legend, 'North', /top, /right, box=0, charsize=1.7
    if (NS EQ 'S') then im_legend, 'South', /top, /right, box=0, charsize=1.7

    for iver=0L,n_sn-1L do begin 
        ; import CMASS P(k)
        cmass_dir = '/mount/riachuelo1/hahn/power/'
        if (drversion EQ 'dr12v4') then $
            person_str = 'Reid' $
            else $
            person_str = 'Anderson'
        if (shotnoises[iver] EQ 'our') then begin 
            cmass_file = 'power-cmass-'+drversion+'-'+NS+'-'+person_str+$
                '-weights-lim-ngalsys-3600lbox-960grid-480bin.dat' 
        endif else if (shotnoises[iver] EQ 'real florian') then begin            ; Florian's real data
            cmass_dir = '/mount/riachuelo1/hahn/literature/beutler_2013/'
            cmass_file = 'Beutler_et_al_2013_cmass_pk.dat'
        endif else begin 
            cmass_file = 'power-cmass-'+drversion+'-'+NS+'-'+person_str+$
                '-weights-zlim-'+shotnoises[iver]+'-shotnoise-ngalsys-3600lbox-960grid-480bin.dat' 
        endelse 
        
        if (shotnoises[iver] NE 'real florian') then begin     
            readcol, cmass_dir+cmass_file, k, Pk, v2, v3, v4, v5, v6, v7, v8, /silent
        endif else begin 
            readcol, cmass_dir+cmass_file, k, Pk, v2, v3, v4, v5, v6, v7, /silent
        endelse 
        
        psym = 16
        symsize = 1.0
        
        djs_oplot, k, Pk, psym=symcat(psym,thick=5), symsize=symsize, line=0, color=im_color(sn_colors[iver])
    endfor 
    im_plotconfig, /psclose, psfile=fig_file, pdf=pdf
    
    fig_file = fig_dir+'power-cmass-'+NS+'-'+strjoin(shotnoises,'-')+$
        '-shotnoise-3600lbox-960grid-480bin-residual.eps'
    im_plotconfig, 0, pos, psfile=fig_file, charsize=1.8

    xrange=[10.^(-1.), 10.^0.]
    yrange=[0.5,1.05]

    ratio_colors = ['green', 'red', 'black']
    ratio_labels = strarr(n_sn-1)
    ;for i=0L,n_sn-2L do ratio_labels[i] = 'CMASS (P(k) '+shotnoises[i+1L]+' shotnoise)/(Standard FKP P(k))'
    for i=0L,n_sn-2L do ratio_labels[i] = 'CMASS (P(k) '+shotnoises[i+1L]+' shotnoise)/(Our FKP P(k))'
    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle=textoidl('k'), $
        ytitle=textoidl('P_{SN}(k)/P_{our FKP}(k)'), /xlog
        ;ytickinterval=0.025,$
    im_legend, ratio_labels, linestyle=replicate(0,n_sn-1), /bottom, /left, box=0, $
        color=ratio_colors, psym=replicate(16,n_sn-1), $
        symsize=replicate(1.0,n_sn)*1.7, symthick=8, spacing=2.7, charsize=1.7
    if (NS EQ 'N') then im_legend, 'North', /top, /right, box=0, charsize=1.7
    if (NS EQ 'S') then im_legend, 'South', /top, /right, box=0, charsize=1.7

    for iver=0L,n_sn-1L do begin 
        ; import CMASS P(k)
        cmass_dir = '/mount/riachuelo1/hahn/power/'
        if (drversion EQ 'dr12v4') then $
            person_str = 'Reid' $
            else $
            person_str = 'Anderson'
        if (shotnoises[iver] EQ 'our') then begin 
            cmass_file = 'power-cmass-'+drversion+'-'+NS+'-'+person_str+$
                '-weights-zlim-ngalsys-3600lbox-960grid-480bin.dat' 
        endif else if (shotnoises[iver] EQ 'real florian') then begin            ; Florian's real data
            cmass_dir = '/mount/riachuelo1/hahn/literature/beutler_2013/'
            cmass_file = 'Beutler_et_al_2013_cmass_pk.dat'
        endif else begin 
            cmass_file = 'power-cmass-'+drversion+'-'+NS+'-'+person_str+$
                '-weights-zlim-'+shotnoises[iver]+'-shotnoise-ngalsys-3600lbox-960grid-480bin.dat' 
        endelse 
        
        if (shotnoises[iver] NE 'real florian') then begin     
            readcol, cmass_dir+cmass_file, k, Pk, v2, v3, v4, v5, v6, v7, v8, /silent
        endif else begin 
            readcol, cmass_dir+cmass_file, k, Pk, v2, v3, v4, v5, v6, v7, /silent
        endelse 
       
        psym = 16
        symsize = 1.0
        
        if (shotnoises[iver] EQ 'our') then begin 
            our_Pk = Pk 
        endif else begin 
            if (shotnoises[iver] NE 'real florian') then begin 
                djs_oplot, k, Pk/our_Pk, $
                    psym=symcat(psym,thick=5), symsize=symsize, line=0, color=im_color(sn_colors[iver])
            endif 
        endelse 
    endfor 
    im_plotconfig, /psclose, psfile=fig_file, pdf=pdf
return 
end 
