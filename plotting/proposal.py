'''

Figures made for proposal.  


'''
import pickle
import numpy as np 
import scipy as sp
import os.path
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from scipy.interpolate import interp1d
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages

# --- Local --- 
from dlos.dlos import Dlos
from fourier_corr import fourier_corr as four_corr 
from corr_spec.corr_corrdata import CorrCorrData
from corr_spec.corr_average import CorrAvgSpec as AvgSpec

# --- plotting ---
from ChangTools.plotting import prettyplot 
from ChangTools.plotting import prettycolors 

# colors for the mocks 
qpm_label = 'QPM' 
qpm_color = prettycolors()[1]
qpm_hatch = 'X'
bmd_label = 'BigMultiDark' 
bmd_color = prettycolors()[7]#'maroon' #prettycolors()[9] #prettycolors()[0]
nsr_label = 'Nseries'
nsr_color = prettycolors()[3]#7]
nsr_hatch = '..'
cms_label = 'BOSS CMASS'
cms_color = 'black' # prettycolors()[9]



def DelP_fc_demo():
    ''' Compare the P(k) residual from FC correction methods to upweighted 
    '''
    # tophat convolution correction parameters
    fs = 0.6 
    rc = 0.43
    fold = 10
    rebin = 20

    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(12,3.8))
    cat_color = nsr_color
    bkgd = fig.add_subplot(111, frameon=False)

    # monopole
    sub_mono = fig.add_subplot(131)

    dlospeak_dict = {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 3.9, 'fpeak': 0.68}
    tru_catdict = {'name': 'nseries', 'n_mock': 1}
    tru_corrdict = {'name': 'true'}
    tru_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': 960, 'ell': 0}
    tru_cat_corr = {'catalog': tru_catdict, 
            'correction': tru_corrdict, 'spec': tru_specdict}
    upw_cat_corr = {'catalog': tru_catdict, 
            'correction': {'name': 'upweight'}, 'spec': tru_specdict}
    dlospeak_cat_corr = {'catalog': tru_catdict, 
            'correction': dlospeak_dict, 'spec': tru_specdict}

    tru_avg_spec = AvgSpec(84, 'pk', tru_cat_corr)
    tru_avg_spec.read(rebin=6) # rebin to 6x
            
    upw_avg_spec = AvgSpec(84, 'pk', upw_cat_corr)
    upw_avg_spec.read(rebin=6)
    
    dlospeak_avg_spec = AvgSpec(84, 'pk', dlospeak_cat_corr)
    dlospeak_avg_spec.read(rebin=6) # rebin to 6x
            
    k_arr = tru_avg_spec.k
    tru_avg_pk = getattr(tru_avg_spec, 'p0k')
    true_mono_interp = interp1d(k_arr, tru_avg_pk, kind='cubic')
    upw_avg_pk = getattr(upw_avg_spec, 'p0k')
    dlospeak_avg_pk = getattr(dlospeak_avg_spec, 'p0k')

    pk_err = tru_avg_spec.stddev(rebin=6)

    samp_var = sub_mono.fill_between(
            k_arr, 
            np.zeros(len(k_arr)) - pk_err, 
            np.zeros(len(k_arr)),# + pk_err, 
            lw=0, 
            color='grey',
            alpha=0.6
            )
    sub_mono.plot((10**-3, 10**0), (0.0,0.0), 'k--', lw=2.5)

    sub_mono.scatter(k_arr, tru_avg_pk- dlospeak_avg_pk, lw=0, color=cat_color)
    
    mock = sub_mono.scatter(k_arr, upw_avg_pk - tru_avg_pk, marker='x', s=20, lw=0.75, color='k') 
    # Axes 
    sub_mono.set_xscale('log')
    sub_mono.set_xlim([5*10**-2,10**0])
    sub_mono.text(0.065, -450., r"$\ell=0$", fontsize=24)
    sub_mono.text(0.6, -450., r"$\mathtt{a)}$", fontsize=24)
    sub_mono.set_ylim([-500.,100.])
    sub_mono.set_yticks([-400, -200, 0])
    sub_mono.minorticks_on()
    sub_mono.set_ylabel(r"$\mathtt{P_l(k)}$ Residuals", fontsize=24)
    
    #sub_mono.legend([mock], 
    #        [r"$\mathtt{\overline{P_l^{FC}} - \overline{P_l^{true}}}$"], 
    #        loc='lower right', scatterpoints=1, 
    #        prop={'size':20}, borderpad=0.5, handletextpad=-0.25, 
    #        markerscale=3., scatteryoffsets=[0.75]) 
        
    # quadrupole 
    sub_quad = fig.add_subplot(132)

    # tophatconv correction parameters hardcoded for 
    # each of the mock catalogs
    cat_label = nsr_label
    cat_color = nsr_color
    # p^upw(k)
    upw_catdict = {'name': 'nseries', 'n_mock': 1}
    upw_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': 960, 'ell': 2}
    upw_cat_corr = {
            'catalog': upw_catdict, 
            'correction': {'name': 'upweight'},  
            'spec': upw_specdict
            }
    upw_avg_spec = AvgSpec(84, 'pk', upw_cat_corr)
    upw_avg_spec.read(rebin=6)
    # p^true(k)
    true_cat_corr = {
            'catalog': upw_catdict, 
            'correction': {'name': 'true'},  
            'spec': upw_specdict
            }
    true_avg_spec = AvgSpec(84, 'pk', true_cat_corr)
    true_avg_spec.read(rebin=6)

    k_arr = upw_avg_spec.k
    upw_avg_pk = getattr(upw_avg_spec, 'p2k')
    true_avg_pk = getattr(true_avg_spec, 'p2k')
    
    pk_err = true_avg_spec.stddev(rebin=6)
    
    samp_var = sub_quad.fill_between(
            k_arr, 
            np.zeros(len(k_arr)) - pk_err, 
            np.zeros(len(k_arr)),# + pk_err, 
            lw=0, color='grey', alpha=0.6)
    sub_quad.plot((10**-3, 10**0), (0.0,0.0), 'k--', lw=2.5)
    # Del P Uncorrelated
    k_uncorr, delP_uncorr = four_corr.DelPuncorr(2, fs=fs, rc=rc, k_arr=k_arr)

    # Del P Correlated
    k_corr, delP_corr = four_corr.DelPcorr_pkmu(2, fs=fs, rc=rc, 
            fold=fold, rebin=rebin, dqdmu=True, dmudq=False)
    corr_interp = interp1d(k_corr, delP_corr, kind='cubic')

    k_range = np.where((k_arr > k_corr[0]) & (k_arr < k_corr[-1]))

    delP_tophat = delP_uncorr[k_range] + corr_interp(k_arr[k_range])
    tophat_plot = sub_quad.scatter(k_arr[k_range], 
            delP_tophat - (upw_avg_pk[k_range] - true_avg_pk[k_range]), 
            lw=0, c=cat_color)
    
    mock = sub_quad.scatter(k_arr, true_avg_pk - upw_avg_pk , 
            marker='x', s=20, lw=1, color='k') #0.75, color='k') # , edgecolor='k') #lw=0.3, facecolors='none', 

    # Axes 
    sub_quad.set_xscale('log')
    sub_quad.set_xlim([5.*10**-2,10**0])
    sub_quad.text(0.065, -450., r"$\ell=2$", fontsize=24)
    sub_quad.text(0.6, -450., r"$\mathtt{b)}$", fontsize=24)
    sub_quad.set_ylim([-500.,100.])
    sub_quad.set_yticklabels([])
    sub_quad.set_yticks([-400, -200, 0])
    sub_quad.minorticks_on()
    #sub_mono.legend(
    #        [samp_var], #[r'$\sigma_{\mathtt{P_l}}$'], 
    #        ['Sample\n Variance'], 
    #        loc='upper right', scatterpoints=1, 
    #        prop={'size':25}, handletextpad=0.2, borderpad=0.0, bbox_to_anchor=(1.1, 1.1)) 
    
    #sub_quad.legend(
    #        [samp_var, tophat_plot, mock], #[r'$\sigma_{\mathtt{P_l}}$'], 
    #        ['\,\, Sample\n\,\, Variance', r"Hahn+(2016)", 'Fib.Collided'], 
    #        loc='upper right', scatterpoints=1, frameon=True,
    #        prop={'size':25}, handletextpad=-0.4, borderpad=0.1, 
    #        ncol=2, columnspacing=-0.1, bbox_to_anchor=(1.0, 1.2)) 
    #samp_leg = sub_quad.legend(
    #        [tophat_plot], 
    #        [r"Hahn+(2016)"],
    #        loc='upper right', scatterpoints=1, 
    #        prop={'size':25}, handletextpad=-0.5, borderpad=0.3, 
    #        markerscale=3., scatteryoffsets=[0.4]) 
    #plt.gca().add_artist(samp_leg)
    #
    #sub_quad.legend([mock], 
    #        ['Fib.Collided'], #r"$\mathtt{\overline{P_l^{FC}} - \overline{P_l^{true}}}$"], 
    #        loc='lower right', scatterpoints=1, 
    #        prop={'size':25}, borderpad=0.2, handletextpad=-0.5, 
    #        markerscale=3., scatteryoffsets=[0.4]) 

    # read in P(k) neutrino mass pickle
    sub_neut = fig.add_subplot(133)
    sub_neut.plot((10**-3, 10**0), (0.0,0.0), 'k--', lw=2.5)
    
    pk_neut_dict = pickle.load(open('/mount/riachuelo1/hahn/power/pk_neut.p', 'rb'))
    pk_no_neut = pk_neut_dict[str(0.0)]
    pk_no_neut_interp = interp1d(pk_no_neut['k'], pk_no_neut['pk'][0], kind='cubic') # interpolate
    bias_no_neut = true_mono_interp(0.05)/pk_no_neut_interp(0.05)
    ks = k_arr[np.where(k_arr > 0.05)]
    pk_no_neut_k_arr = pk_no_neut_interp(ks)
    print len(k_arr)
    
    k_text = [0.065, 0.08, 0.08, 0.25]
    pk_text = [-120., -280., -440., -400.]
    for i_mnu, mnu in enumerate([0.1, 0.2, 0.3, 0.4]): 
        pk_neut = pk_neut_dict[str(mnu)]
        #print pk_neut['pk'][0]
        pk_neut_interp = interp1d(pk_neut['k'], pk_neut['pk'][0], kind='cubic') # interpolate
        bias_neut = true_mono_interp(0.05)/pk_neut_interp(0.05)
        print bias_no_neut, bias_neut
        sub_neut.plot(ks,  bias_neut * pk_neut_interp(ks) - bias_no_neut * pk_no_neut_k_arr, 
                c=pretty_colors[2*i_mnu+1], lw=2)
        if i_mnu == 0: 
            sub_neut.text(k_text[i_mnu], pk_text[i_mnu], r'$\Sigma m_\nu = '+str(mnu)+'\mathtt{eV}$', 
                    fontsize=19, color=pretty_colors[2*i_mnu+1])
        else: 
            sub_neut.text(k_text[i_mnu], pk_text[i_mnu], r'$'+str(mnu)+'\mathtt{eV}$', 
                    fontsize=20, color=pretty_colors[2*i_mnu+1])
    sub_neut.text(0.07, 10, r'Impact\,of\,$\Sigma m_\nu$\,on\,$\mathtt{P_m}$', fontsize=20)

    sub_neut.set_xscale('log')
    sub_neut.set_xlim([5.*10**-2,10**0])
    sub_neut.text(0.6, -450., r"$\mathtt{c)}$", fontsize=24)
    sub_neut.set_ylim([-500.,100.])
    sub_neut.set_yticklabels([])
    sub_neut.set_yticks([-400, -200, 0])
    sub_neut.minorticks_on()

    fig.subplots_adjust(wspace=0.05, hspace=0.)
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_xlabel('k (h/Mpc)', fontsize=24)

    fig_file = ''.join(['figure/fc_paper/proposal_pk_fc_demo.png'])
    fig.savefig(fig_file, bbox_inches="tight", dpi=500)
    plt.close()
    return None


if __name__=='__main__': 
    DelP_fc_demo()
