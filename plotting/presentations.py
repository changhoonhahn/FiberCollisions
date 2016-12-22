'''

Figures for talks and posters

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


def Pk_residual_NN(catalogs='nseries', n_mocks=84, Ngrid=960, rebin=None): 
    ''' Figure that compares the fiber collision residuals Del P_l(k) 
    of the NN upweighting method for l = 0, 2, 4 for a specified 
    mock catalog
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14, 5), dpi=1000)
    all_fig = fig.add_subplot(111, frameon=False)

    for i_ell, ell in enumerate([0, 2, 4]):  # monopole or quadrupole
        sub = fig.add_subplot(1, 3, i_ell+1)
        if catalogs == 'nseries': 
            cat_label = nsr_label 
            cat_color = nsr_color 
        else: 
            raise NotImplementedError 

        tru_catdict = {'name': catalogs, 'n_mock': 1}
        tru_corrdict = {'name': 'true'}
        upw_corrdict = {'name': 'upweight'}
        tru_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell}
        tru_cat_corr = {
                'catalog': {'name': tru_catdict['name'], 'n_mock': 1}, 
                'correction': tru_corrdict, 
                'spec': tru_specdict
                }
        upw_cat_corr = {
                'catalog': {'name': tru_catdict['name'], 'n_mock': 1}, 
                'correction': upw_corrdict, 
                'spec': tru_specdict
                }
        now_cat_corr = {
            'catalog': {'name': tru_catdict['name'], 'n_mock': 1}, 
            'correction': {'name': 'noweight'}, 
            'spec': tru_specdict
            }
        tru_avg_spec = AvgSpec(n_mocks, 'pk', tru_cat_corr)
        tru_avg_spec.read(rebin=rebin)
        
        upw_avg_spec = AvgSpec(n_mocks, 'pk', upw_cat_corr)
        upw_avg_spec.read(rebin=rebin)

        now_avg_spec = AvgSpec(n_mocks, 'pk', now_cat_corr)
        now_avg_spec.read(rebin=rebin)
        now_avg_pk = getattr(now_avg_spec, ''.join(['p', str(ell), 'k']))
    
        k_arr = tru_avg_spec.k
        tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
        upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
        
        # plot the sample variance
        if n_mocks > 1: 
            pk_err = tru_avg_spec.stddev(rebin=rebin)

            samp_var = sub.fill_between(
                    k_arr, 
                    np.zeros(len(k_arr)) - pk_err, 
                    np.zeros(len(k_arr)) + pk_err, 
                    lw=0, 
                    color='grey',
                    alpha=0.6
                    )
        sub.plot((10**-3, 10**0), (0.0,0.0), 'k--')

        if ell in [0,2]: 
            plt_label = None
        elif ell == 4: 
            plt_label = cat_label

        mock = sub.scatter(k_arr, upw_avg_pk - tru_avg_pk, lw=0, color=cat_color)
        # P(k) with no weights 
        #noweight_plot, = sub.plot(k_arr, now_avg_pk - tru_avg_pk, 
        #        lw=1.2, ls='--', color='k', dashes=[10,5])
        noweight_plot = sub.scatter(k_arr, now_avg_pk - tru_avg_pk, 
                marker='x', s=16, lw=0.75, color='k') # , edgecolor='k') #lw=0.3, facecolors='none', 
        noweight_label = r"$\Delta \mathtt{P_l^{NoW}(k)}$"

        # Axes 
        sub.set_xlim([10**-2,10**0])
        sub.set_xscale('log')

        sub.text(0.2, 0.85, r"$\ell\mathtt{= "+str(ell)+"}$", fontsize=24, 
                ha='center', va='center', transform=sub.transAxes)
        if ell == 0: 
            sub.set_ylim([-750.,750.])
            plt.xticks([10**-2, 10**-1])
            plt.sca(sub)
            plt.yticks(list(np.arange(-750., 1000., 250.)) )
        elif ell == 2: 
            sub.set_ylim([-1000.,1000.])
            plt.yticks(list(np.arange(-1000., 1500., 500.)) )
            plt.xticks([10**-2, 10**-1])
        elif ell == 4: 
            sub.set_ylim([-500.,500.])
            plt.sca(sub)
            #plt.yticks(list(np.arange(-500., 1000., 500.)) )
            plt.xticks([10**-2, 10**-1])
            sub.text(0.7, -400., 'Nseries', fontsize=20, horizontalalignment='right') 
        plt.sca(sub)
        plt.xticks([10**-2, 10**-1, 10**0])
                
    all_fig.set_xticklabels([])
    all_fig.set_yticklabels([])
    all_fig.set_ylabel(r"$\Delta \mathtt{P_l(k)}$", fontsize=24, labelpad=50)
    all_fig.set_xlabel('k (h/Mpc)', fontsize=24, labelpad=20)
    #r"$\mathtt{\overline{P_l^{NN}} - \overline{P_l^{true}}(k)}$", 
    fig.subplots_adjust(wspace=0.35, hspace=0.0)
    
    if rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(rebin)+'x'
    fig_file = ''.join(['figure/fc_paper/',
        'Pk_residual_NN', rebin_str, '.catalog_', catalogs, '.png'])
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()
    return None

# LOS Reconstruction Residual plot
def Pk_residual_LOSrecon(catalogs='nseries', n_mocks=84, Ngrid=960, florian_offset=None, rebin=None): 
    ''' Figure that compares the P_l(k) residual of the line-of-sight reconstructed 
    power spectrum monopole and quadrupole of specified mock catalogs. 
    
    P_l(k)^dlospeak - P_l(k)^true l = 0, 2 for 

    Parameters
    ----------
    catalogs : list
        string that specify the mock catalogs to include. 
        (e.g. 'nseries')

    n_mocks : list 
        int of mocks to average the P_l(k) over. Same length as 
        the list of catalog names. 

    Ngrid : int 
        int that specify the Ngrids of the power spectrum 
        monopole and quadrupole. Ngrid list has two elements. 
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(10,5))
    all_fig = fig.add_subplot(111, frameon=False)

    dlospeak_dict = {
            'qpm': {
                'name': 'dlospeak', 
                'fit': 'gauss', 
                'sigma': 4.4, 
                'fpeak': 0.62
                }, 
            'bigmd': {
                'name': 'dlospeak', 
                'fit': 'gauss', 
                'sigma': 5.5, 
                'fpeak': 0.6
                },
            'nseries': {
                'name': 'dlospeak', 
                'fit': 'gauss', 
                'sigma': 3.9, 
                'fpeak': 0.68
                }
            }

    for i_ell, ell in enumerate([0, 2]):  # monopole or quadrupole
        sub = fig.add_subplot(1, 2, i_ell+1)
        
        # dLOS peak correction parameters hardcoded for 
        # each of the mock catalogs
        if catalogs == 'qpm': 
            cat_label = qpm_label 
            cat_color = qpm_color
        elif catalogs == 'bigmd': 
            cat_label = bmd_label 
            cat_color = bmd_color
        elif catalogs == 'nseries': 
            cat_label = nsr_label
            cat_color = nsr_color
        else: 
            raise NotImplementedError 
        dlospeak_corrdict = dlospeak_dict[catalogs]  

        tru_catdict = {'name': catalogs, 'n_mock': 1}
        tru_corrdict = {'name': 'true'}
        tru_specdict = {
                'P0': 20000,
                'Lbox': 3600, 
                'Ngrid': Ngrid, 
                'ell': ell 
                }
        tru_cat_corr = {
                'catalog': tru_catdict, 
                'correction': tru_corrdict, 
                'spec': tru_specdict
                }
        dlospeak_cat_corr = {
                'catalog': tru_catdict, 
                'correction': dlospeak_corrdict, 
                'spec': tru_specdict
                }
        florian_cat_corr = {
                'catalog': tru_catdict, 
                'correction': {'name': 'floriansn'}, 
                'spec': tru_specdict
                }
        hector_cat_corr = {
                'catalog': tru_catdict, 
                'correction': {'name': 'hectorsn'}, 
                'spec': tru_specdict
                }

        tru_avg_spec = AvgSpec(n_mocks, 'pk', tru_cat_corr)
        tru_avg_spec.read(rebin=rebin)
        
        dlospeak_avg_spec = AvgSpec(n_mocks, 'pk', dlospeak_cat_corr)
        dlospeak_avg_spec.read(rebin=rebin)
        
        k_arr = tru_avg_spec.k
        tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
        dlospeak_avg_pk = getattr(dlospeak_avg_spec, ''.join(['p', str(ell), 'k']))
    
        # plot sample variance 
        if n_mocks > 1: 
            pk_err = tru_avg_spec.stddev(rebin=rebin)

            samp_var = sub.fill_between(
                    k_arr, 
                    np.zeros(len(k_arr)) - pk_err, 
                    np.zeros(len(k_arr)) + pk_err, 
                    lw=0, 
                    color='grey',
                    alpha=0.6
                    )
        sub.plot((10**-3, 10**0), (0.0,0.0), 'k--')

        if ell == 0: 
            plt_label = None
        elif ell == 2: 
            plt_label = cat_label
        
        # residual 
        mock = sub.scatter(k_arr, dlospeak_avg_pk - tru_avg_pk, lw=0, color=cat_color)

        if ell == 0: 
            florian_avg_spec = AvgSpec(n_mocks, 'pk', florian_cat_corr)
            florian_avg_spec.read(rebin=rebin)
            florian_avg_pk = getattr(florian_avg_spec, ''.join(['p', str(ell), 'k']))
            #florian, = sub.plot(k_arr, florian_avg_pk - tru_avg_pk, 
            #        lw=1., ls='--', color='k', dashes=[10,7])
            if florian_offset is None: 
                florian_offset = 0 

            #florian = sub.scatter(k_arr, florian_avg_pk - tru_avg_pk - florian_offset, 
            #        marker='+', s=20, lw=0.7, color='k') # , edgecolor='k') #lw=0.3, facecolors='none', 

            hector_avg_spec = AvgSpec(n_mocks, 'pk', hector_cat_corr)
            hector_avg_spec.read(rebin=rebin)
            hector_avg_pk = getattr(hector_avg_spec, ''.join(['p', str(ell), 'k']))
            #hector = sub.scatter(k_arr, hector_avg_pk - tru_avg_pk, 
            #        s=8, lw=0.5, facecolors='none', edgecolor='k')#, dashes=[10,5])
            hector, = sub.plot(k_arr, hector_avg_pk - tru_avg_pk, 
                    lw=1., ls='--', color='k', dashes=[10,7])
        # Axes 
        sub.set_xlim([10**-2,10**0])
        sub.set_xscale('log')

        sub.text(0.2, 0.85, r"$\ell\mathtt{= "+str(ell)+"}$", fontsize=24, 
                ha='center', va='center', transform=sub.transAxes)

        if ell == 0: 
            sub.set_ylim([-750.,750.])
            plt.sca(sub)
            plt.yticks(list(np.arange(-750., 1000., 250.)) )
            plt.xticks([10**-2, 10**-1, 10**0])
        elif ell == 2: 
            sub.set_ylim([-1000.,1000.])
            plt.yticks(list(np.arange(-1000., 1500., 500.)) )
            plt.xticks([10**-2, 10**-1, 10**0])
            sub.text(0.8, -850., cat_label, fontsize=20, horizontalalignment='right') 

    all_fig.set_xticklabels([])
    all_fig.set_yticklabels([])
    all_fig.set_ylabel(r"$\Delta \mathtt{P_l(k)}$", fontsize=24, labelpad=50)
    all_fig.set_xlabel(r"k (h/Mpc)", fontsize=24, labelpad=20)

    fig.subplots_adjust(wspace=0.4, hspace=0.0)
    
    if florian_offset is None: 
        florian_str = ''
    else: 
        florian_str = '.floran_offset'+str(round(florian_offset,1))
    if rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(rebin)+'x'
    fig_file = ''.join(['figure/fc_paper/', 
        'Plk_residual_LOSrecon', florian_str, rebin_str, 
        '.catalog_', catalogs, '.png'])
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()

def DelP_untrusted_EffWin(ktrust, n_mocks=[84, 84], Ngrid=[960, 480], order=10, k_rebin=None, figdata=False): 
    ''' Figure that compares Del P^corr integrated from q = ktrust to inf
    using the polynomial expansion of P(k, mu) to the empirical Del P 
    minus Del P^uncorr minus Del P^corr integrated from q = 0 to ktrust
    using P(k, mu). 
    
    Parameters
    ----------
    n_mocks : list 
        Length 2 list of ints that specify the number of mocks to average 
        the P_l(k) over. 

    Ngrid : list 
        Length 2 list of ints that specify the Ngrids of the power spectrum 
        monopole and quadrupole, respectively
    '''
    # tophat convolution correction parameters
    fs = 0.6 
    rc = 0.43
    fold = 10
    rebin = 20

    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14,6))

    fig_data = [] 
    for i_ell, ell in enumerate([0, 2]):  # monopole or quadrupole
        sub = fig.add_subplot(1, 2, i_ell+1)

        # TophatConv correction parameters hardcoded for 
        # each of the mock catalogs
        cat_label = nsr_label
        cat_color = nsr_color
        tophat_conv_corrdict = {
                'name': 'tophat_conv', 
                }
        # P^upw(k)
        upw_catdict = {'name': 'nseries', 'n_mock': 1}
        upw_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid[i_ell], 'ell': ell}
        upw_cat_corr = {
                'catalog': upw_catdict, 
                'correction': {'name': 'upweight'},  
                'spec': upw_specdict
                }
        upw_avg_spec = AvgSpec(n_mocks[i_ell], 'pk', upw_cat_corr)
        upw_avg_spec.read(rebin=k_rebin)
        # P^true(k)
        true_cat_corr = {
                'catalog': upw_catdict, 
                'correction': {'name': 'true'},  
                'spec': upw_specdict
                }
        true_avg_spec = AvgSpec(n_mocks[i_ell], 'pk', true_cat_corr)
        true_avg_spec.read(rebin=k_rebin)

        k_arr = upw_avg_spec.k
        upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
        true_avg_pk = getattr(true_avg_spec, ''.join(['p', str(ell), 'k']))
        
        delP_emp = upw_avg_pk - true_avg_pk 
        for imock in range(1, n_mocks[i_ell]+1): 
            upwP_i = upw_avg_spec.spec_i(imock, rebin=k_rebin)[1]
            truP_i = true_avg_spec.spec_i(imock, rebin=k_rebin)[1]
            delP_i = upwP_i - truP_i 
           
            if imock == 1:
                delP_var = (delP_emp - delP_i)**2
            else: 
                delP_var += (delP_emp - delP_i)**2
            
        delP_stddev = np.sqrt(delP_var)
        
        # Del P Uncorrelated
        k_uncorr, delP_uncorr = four_corr.DelPuncorr(ell, fs=fs, rc=rc, k_arr=k_arr)

        # Del P Corr Trusted Range
        k_corr_kt, delP_corr_kt = four_corr.DelPcorr_pkmu(
                ell, fs=fs, rc=rc, fold=fold, rebin=rebin, 
                dqdmu=True, dmudq=False, ktrust=ktrust)
        corr_kt_interp = interp1d(k_corr_kt, delP_corr_kt, kind='cubic')
        k_range = np.where((k_arr > k_corr_kt[0]) & (k_arr < k_corr_kt[-1]))

        # Empirical Del P
        delP_emp_untrust = delP_emp[k_range] - delP_uncorr[k_range] - corr_kt_interp(k_arr[k_range])

        #- delP_uncorr[k_range] - corr_kt_interp(k_arr[k_range])
        samp_var = sub.fill_between(
                k_arr[k_range], 
                delP_emp_untrust - delP_stddev[k_range], 
                delP_emp_untrust + delP_stddev[k_range], 
                lw=0, 
                color='grey',
                alpha=0.6)
        emp_plot, = sub.plot(k_arr[k_range], delP_emp_untrust, lw=2, c='k', ls='--')

        # Del P untrusted from integrated coefficients and polynomial
        lp_lim = np.where(k_arr < 0.35)
        delpcorr_untrust = four_corr.DelPcorr_untrust_poly(
                k_arr, ell=ell, ktrust=ktrust, order='all', 
                fs=fs, rc=rc, fold=10, rebin=50, max_order=order)         
        #tophat_plot = sub.scatter(k_arr[lp_lim], delpcorr_untrust[lp_lim], 
        #        lw=0, s=8, c=cat_color)
    
        delpcorr_untrust_orderlp = four_corr.DelPcorr_untrust_poly(
                k_arr, ell=ell, ktrust=ktrust, order='upto2', 
                fs=fs, rc=rc, fold=10, rebin=50, max_order=order)

        poly2_plot, = sub.plot(k_arr[lp_lim], delpcorr_untrust_orderlp[lp_lim], 
                lw=2, ls='-', c=pretty_colors[1])
            
        # x-axis 
        sub.set_xscale('log')
        sub.set_xlim([10**-2,0.5])
        sub.set_xlabel('k (h/Mpc)', fontsize=24)
        # y-axis 
        if ell == 0: 
            sub.text(1.5*10**-2, 125., r"$\ell{=} 0$", fontsize=24)
            #sub.set_ylim([-600., 0.])
            sub.set_ylim([-1100., 300.])
            sub.set_ylabel(
                    r"$\Delta \mathtt{P^{corr}_l(k)} |_\mathtt{q=k_{trust}}^{\mathtt{q=}\infty}$", 
                    fontsize=24)
        elif ell == 2: 
            sub.text(1.5*10**-2, 900., r"$\ell{=} 2$", fontsize=24)
            #sub.set_ylim([-500., 500.])
            sub.set_ylim([-1200., 1200.])
            #sub.set_ylabel(r"$\mathtt{P_2(k)}$ Residuals", fontsize=24)#, labelpad=20)
            sub.yaxis.tick_right()
            #sub.yaxis.set_label_position('right')
    if k_rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(k_rebin)+'x'

    fig.subplots_adjust(wspace=0.15, hspace=0.15)

    fig_file = 'figure/fc_paper/DelPlk_untrusted_EffWin_'+str(order)+rebin_str+'.png' 
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()
    return None


if __name__=='__main__': 
    #Pk_residual_NN(catalogs='nseries', n_mocks=84, Ngrid=960, rebin=6)
    #Pk_residual_LOSrecon(catalogs='nseries', n_mocks=84, Ngrid=960, florian_offset=250., rebin=6)
    DelP_untrusted_EffWin(0.3, n_mocks=[84, 84], Ngrid=[960, 960], order=18, k_rebin=6)
