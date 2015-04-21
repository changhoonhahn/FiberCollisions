pro plot_nbar, tilingmock=tilingmock
; plots nbar(z) for CMASS and BOSS5003 Tiling mock
    fig_dir = '/home/users/hahn/research/figures/boss/powerspectrum/tiling_mock/'
    psfile = fig_dir+'fig-nbar-cmass-tilingmock-comp.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8

    tm_color='dodger blue' & tm_psym=16 & tm_symsize=0.75 & tm_label='Tiling Mocks NGC'
    cmass_color='black' & cmass_psym=0 & cmass_symsize=0 & cmass_label='CMASS DR11 May22'

    xrange = [0.0, 1.0]
    yrange = [0.0, 0.0005]
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle=textoidl('k'), $
        ytitle=textoidl('\bar{n}(k)')
    im_legend, [tm_label,cmass_label], linestyle=[0,2], $
        /right, /top, box=0, color=[tm_color, cmass_color], $
        psym=[tm_psym, cmass_psym], symsize=[tm_symsize, cmass_symsize]*1.7, $
        symthick=8, spacing=2.7, charsize=2.1

; import CMASS nbar(z) file
    cmass_dir = '/mount/riachuelo1/hahn/data/manera_mock/dr11/'
    cmass_file = 'nbar-cmass-dr11may22-N-Anderson.dat'
    readcol, cmass_dir+cmass_file, zmid,zlow,zhigh,nbarz,wfkp,shell_vol,galtot
; import Tiling mock nbar(z) file
    tm_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'
    tm_file = 'nbar-cmass-boss5003sector-icoll012.dat'
    readcol, tm_dir+tm_file, tm_zmid,tm_zlow,tm_zhigh,tm_nbarz,tm_wfkp,tm_shell_vol,tm_galtot

    djs_oplot, zmid, nbarz, $
        line=2, thick=6, color=im_color(cmass_color)
    djs_oplot, tm_zmid, tm_nbarz,$
      psym=symcat(tm_psym,thick=5),symsize=tm_symsize, color=im_color(tm_color)
return
end 
