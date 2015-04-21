pro plot_power_lasdamasbox_aniso, n 
; plot P(k) for LasDamas Box simulation using different 
; weighting schemes
    fig_dir = '/home/users/hahn/research/figures/boss/powerspectrum/las_damas/box/'
    lasdamas_dir = get_fibcol_path('box/', /ria, /power, /lasdamas)
; import LasDamas Box mocks P(k)
    lasdamas_file = 'power_mock_gamma_lrg21p2_zmo_Oriana_1'+$
        strmid(strtrim(string(1000+n),2),1,4)+'_z0p342_fof_b0p2.zdist.dat.grid360.P0.box3600'
    pknb_file = 'power_mock_gamma_lrg21p2_zmo_Oriana_1'+$   ; peak nbar correction file
        strmid(strtrim(string(1000+n),2),1,4)+'_z0p342_fof_b0p2.zdist.fibcoll.peak.dat.grid360.P0.box3600'
    del_file = 'power_mock_gamma_lrg21p2_zmo_Oriana_1'+$    ; delta function correction file
        strmid(strtrim(string(1000+n),2),1,4)+'_z0p342_fof_b0p2.zdist.fibcoll.delta.dat.grid360.P0.box3600'
    readcol, lasdamas_dir+lasdamas_file, ld_k, ld_Pk, v2, v3, v4, /silent
    readcol, lasdamas_dir+pknb_file, pknb_k, pknb_Pk, v2, v3, v4, /silent
    readcol, lasdamas_dir+del_file, del_k, del_Pk, v2, v3, v4, /silent
; import cmass DR11v2 P(k) 
    readcol, '/mount/riachuelo1/hahn/power/power-cmass-dr11v2-N-Anderson-weights-zlim-ngalsys-3600lbox-360grid-180bin.dat',$
        cmass_k, cmass_Pk, v2, v3, v4, v5, v6, v7, v8, /silent
; configure plot
    psfile = fig_dir+'power-lasdamasbox-correct-comp-3600lbox-360grid-180bin.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8

    ld_color='grey' & ld_psym=16 & ld_symsize=1.0 & ld_label='Las Damas Box'
    pknb_color='red' & pknb_psym=16 & pknb_symsize=1.0 & pknb_label='Peak Correction'
    del_color='blue' & del_psym=16 & del_symsize=1.0 & del_label=textoidl('\Delta')+' Correction'
    cmass_color='black' & cmass_psym=16 & cmass_symsize=1.0 & cmass_label='CMASS DR11 v11.2'

    xrange=[10.^(-3.),10.^0.]
    yrange=[10.^2.,10.^5.5]
    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle=textoidl('k'), $
        ytitle=textoidl('P(k)'), /xlog, /ylog
    im_legend, [ld_label, del_label, pknb_label, cmass_label], linestyle=[0,0,0,0], /bottom, /left, box=0, $
        color=[ld_color, del_color, pknb_color, cmass_color], psym=[ld_psym, del_psym, pknb_psym, cmass_psym], $
        symsize=[ld_symsize, del_symsize, pknb_symsize, cmass_symsize]*1.7, symthick=8, spacing=2.7, charsize=1.7
    im_legend, 'Las Damas Box', /top, /right, box=0, charsize=1.7

    djs_oplot, ld_k, factor*ld_Pk, $
        psym=symcat(ld_psym,thick=5), symsize=ld_symsize, line=0, color=im_color(ld_color)
    djs_oplot, pknb_k, factor*pknb_Pk, $
        line=2, thick=6, color=im_color(pknb_color)
    djs_oplot, del_k, factor*del_Pk, $
        line=2, thick=6, color=im_color(del_color)
    djs_oplot, cmass_k, cmass_Pk, $
        line=2, thick=6, color=im_color(cmass_color)
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf

    psfile = fig_dir+'power-lasdamasbox-correct-comp-ratio-3600lbox-360grid-180bin.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8
    yrange = [0.8, 1.3]
    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle=textoidl('k'), $
        ytitle=textoidl('P(k)/P(k)_{True}'), /xlog
    im_legend, [pknb_label, del_label], linestyle=[0,0], $
        /left, /bottom, box=0, color=[pknb_color, del_color], $
        psym=[pknb_psym, del_psym], symsize=[pknb_symsize, del_symsize]*1.7,$
        symthick=8, spacing=2.7, charsize=2.1
    im_legend, 'Las Damas Box', /top, /right, box=0, charsize=1.7
    djs_oplot, pknb_k, pknb_Pk/ld_Pk, $
        line=2, thick=6, color=im_color(pknb_color)
    djs_oplot, del_k, del_Pk/ld_Pk, $
        line=2, thick=6, color=im_color(del_color)
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
return 
end 
