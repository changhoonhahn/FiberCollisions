pro plot_power_tiling_mocks
    fig_dir = '/home/users/hahn/research/figures/boss/powerspectrum/tiling_mock/'
    tm_dir = '/mount/riachuelo1/hahn/power/tiling_mocks/'
    cmass_dir   = '/mount/riachuelo1/hahn/power/'

; import tiling mock P(k)
    tm_file = tm_dir+$
        'power_cmass-boss5003sector-icoll012.dat.grid360.P020000.box4000'
    pknb_file = tm_dir+$    ; peak nbar correction file
        'power_cmass-boss5003sector-icoll012.fibcoll.dat.peaknbarcorr.grid360.P020000.box4000'
    del_file = tm_dir+$    ; delta function correction file
        'power_cmass-boss5003sector-icoll012.fibcoll.dat.deltacorr.grid360.P020000.box4000'
    readcol, tm_file, tm_k, tm_Pk, v2, v3, v4, v5, v6, v7, v8, /silent
    readcol, pknb_file, pknb_k, pknb_Pk, v2, v3, v4, v5, v6, v7, v8, /silent
    readcol, del_file, del_k, del_Pk, v2, v3, v4, v5, v6, v7, v8, /silent
; import cmass DR11v2 P(k) 
    readcol, '/mount/riachuelo1/hahn/power/power-cmass-dr11v2-N-Anderson-weights-zlim-ngalsys-3600lbox-360grid-180bin.dat',$
        cmass_k, cmass_Pk, v2, v3, v4, v5, v6, v7, v8, /silent
; configure plot
    psfile = fig_dir+'power-tilingmocks-correct-comp-4000lbox-360grid-180bin.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8

    tm_color='grey' & tm_psym=16 & tm_symsize=0.75 & tm_label='Tiling Mocks NGC'
    pknb_color='red' & pknb_psym=16 & pknb_symsize=1.5 & pknb_label='Peak '+$
        textoidl("n(z)")+' Correction NGC'
    del_color='blue' & del_psym=16 & del_symsize=1.5 & del_label='Delta Function Correction NGC'
    cmass_color='black' & cmass_psym=0 & cmass_symsize=0 & cmass_label='CMASS DR11 v11.2'
    tm_del_residual_color = 'black' & tm_del_residual_label=textoidl('P(k)_{\Delta}-P(k)_{TM}')

    xrange=[10.^(-3.),10.^0.]
    yrange=[10.^2.,10.^5.5]
    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle=textoidl('k'), $
        ytitle=textoidl('P(k)'), /xlog, /ylog
    im_legend, [tm_label, pknb_label, del_label, cmass_label, tm_del_residual_label], linestyle=[0,0,0,2,0], $
        position=[0.001, 10.^4], box=0, color=[tm_color, pknb_color, del_color, cmass_color,tm_del_residual_color], $
        psym=[tm_psym, pknb_psym, del_psym, cmass_psym, 16], symsize=[tm_symsize, pknb_symsize, $
        del_symsize, cmass_symsize, 1.5]*1.7, symthick=8, spacing=2.7, charsize=1.7

    djs_oplot, tm_k, tm_Pk, $
        psym=symcat(tm_psym,thick=5), symsize=tm_symsize, line=0, color=im_color(tm_color)
    djs_oplot, pknb_k, pknb_Pk, $
        line=2, thick=6, color=im_color(pknb_color)
;        psym=symcat(pknb_psym,thick=5), symsize=pknb_symsize, line=0, color=im_color(pknb_color)
    djs_oplot, del_k, del_Pk, $
        line=2, thick=6, color=im_color(del_color)
;        psym=symcat(del_psym,thick=5), symsize=del_symsize, line=0, color=im_color(del_color)
    djs_oplot, del_k, del_Pk-tm_Pk, $
        line=0, thick=6, color=im_color('black')
    djs_oplot, cmass_k, cmass_Pk, $
        line=2, thick=4, color=im_color(cmass_color)
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf

    psfile = fig_dir+'power-tilingmocks-correct-delta-comp-4000lbox-360grid-180bin.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8

    yrange = [0.8, 1.3]
    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle=textoidl('k'), $
        ytitle=textoidl('P(k)/P(k)_{tm}'), /xlog
    im_legend, [pknb_label, del_label], linestyle=[0,0], $
        /left, /bottom, box=0, color=[pknb_color, del_color], $
        psym=[pknb_psym, del_psym], symsize=[pknb_symsize, del_symsize]*1.7,$
        symthick=8, spacing=2.7, charsize=2.1
    djs_oplot, pknb_k, pknb_Pk/tm_Pk, $
        line=2, thick=6, color=im_color(pknb_color)
;        psym=symcat(pknb_psym,thick=5), symsize=pknb_symsize, line=0, color=im_color(pknb_color)
    djs_oplot, del_k, del_Pk/tm_Pk, $
        line=2, thick=6, color=im_color(del_color)
;        psym=symcat(del_psym,thick=5), symsize=del_symsize, line=0, color=im_color(del_color)
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
return 
end 
