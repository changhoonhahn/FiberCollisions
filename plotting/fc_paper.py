'''

Figures for the Fiber Collisions paper

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

def fibcol_fraction(catalogs=['nseries'], n_mocks=[84]): 
    ''' Calculate the fiber collision fraction of the different 
    catalogs
    '''
    for i_cat, cat in enumerate(catalogs):
        n_tot = 0.
        n_fc = 0. 

        for i_mock in range(1, n_mocks[i_cat]+1): 
            cat_corr = {
                    'catalog': {'name': cat, 'n_mock': i_mock},
                    'correction': {'name': 'upweight'}
                    }
            gal_data = CorrCorrData(cat_corr)
            gal_data.read()
    
            n_tot += np.float(len(gal_data.z))
            n_fc += np.sum(gal_data.wfc == 0)
    
        print cat
        print n_fc / n_tot 
    
    return None 

def fig_mock_catalogs(catalogs=['nseries'], n_mocks=[84]): 
    ''' Figure that compares the redshift distribution of mock catalogs 
    with the BOSS DR12 CMASS sample galaxies. 

    Parameters
    ----------
    catalogs : list
        List of strings that specify the mock catalogs to include. 
        (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    n_mocks : list 
        List of ints that specify the number of mocks for each of the 
        mock catalogs 
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(8,7))
    sub = fig.add_subplot(111)
    
    # CMASS DR12v4
    z_mid, z_low, z_high, cmass_ngal_eff = np.loadtxt(
            '/mount/riachuelo1/hahn/data/nbar-cmass-dr12v4-N-Reid.dat', 
            unpack=True, 
            usecols=[0,1,2,6]
            )
    fig_data = [] 
    for i_cat, cat in enumerate(catalogs):
        for i_mock in range(1, n_mocks[i_cat]+1): 
            if cat == 'qpm': 
                cat_label = qpm_label 
                cat_color = qpm_color
            elif cat == 'bigmd': 
                cat_label = bmd_label
                cat_color = bmd_color
            elif cat == 'nseries': 
                cat_label = nsr_label
                cat_color = nsr_color
            else: 
                raise NotImplementedError 

            cat_corr = {
                    'catalog': {'name': cat, 'n_mock': i_mock},
                    'correction': {'name': 'true'}
                    }
            gal_data = CorrCorrData(cat_corr)
            gal_data.read()
    
            if i_mock == 1: 
                ngal_eff = np.zeros(len(z_mid))

            for i_z in xrange(len(z_mid)):
                zbin = np.where( (gal_data.z > z_low[i_z]) & (gal_data.z <= z_high[i_z]) )
                ngal_eff[i_z] += np.sum(gal_data.wfc[zbin])
        
        norm_ngal_eff = ngal_eff / np.sum(ngal_eff)
        sub.step(z_mid, norm_ngal_eff, lw=4, color=cat_color, label=cat_label, where='mid') 

        fig_data.append(
                {'catalog': cat, 'zmid': z_mid, 'norm_ngal_eff': norm_ngal_eff}
                )

    zlim = np.where((z_mid > 0.43) & (z_mid < 0.7))
    sub.step(
            z_mid[zlim], 
            cmass_ngal_eff[zlim] / np.sum(cmass_ngal_eff[zlim]), 
            lw=4, 
            color=cms_color, 
            label=cms_label, 
            where='mid'
            ) 

    sub.set_xlabel(r"Redshift $(\mathtt{z})$", fontsize=30) 
    sub.set_xlim([0.4, 0.75])
    sub.set_ylim([0.0, 0.045])
    sub.set_ylabel(r"$\mathtt{p(z)}$", fontsize=25)
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':25}) 
    
    fig_file = 'figure/fc_paper/mock_catalog_z_dist.png'
    fig_data_file = 'figure/fc_paper/mock_catalog_z_dist.figdata.p'
    pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()
    return None

def fig_p0k_p2k_mock_catalogs(catalogs=['nseries'], n_mocks=[84], Ngrid=[960], rebin=None): 
    ''' Figure that compares the powerspectrum  monopole and quadrupole of mock 
    catalogs 

    Parameters
    ----------
    catalogs : list 
        List of strings that specify the mock catalogs to include. (e.g. 
        ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    n_mocks : list
        List of ints that specify the number of mock catalogs to average over
    Ngrid : list 
        List of ints that specify Ngrid for ell terms. 
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14,8))
    mono = fig.add_subplot(121)
    quad = fig.add_subplot(122)
    
    if catalogs[-1] != 'cmass': 
        catalogs.append('cmass')
        n_mocks.append(1)

    fig_data = [] 
    for i_cat, cat in enumerate(catalogs):
        if cat == 'qpm': 
            cat_label = qpm_label 
            cat_color = qpm_color
            cat_hatch = qpm_hatch
        elif cat == 'bigmd': 
            cat_label = bmd_label 
            cat_color = bmd_color
        elif cat == 'nseries': 
            cat_label = nsr_label
            cat_color = nsr_color
            cat_hatch = nsr_hatch
        elif cat == 'cmass': 
            cat_label = cms_label
            cat_color = cms_color
        else: 
            raise NotImplementedError 
        
        for i_ell, ell in enumerate([0, 2]): 
            catdict = {'name': cat, 'n_mock': 1}
            corrdict = {'name': 'true'}
            specdict = {
                    'P0': 20000,
                    'Lbox': 3600, 
                    'Ngrid': Ngrid[i_ell], 
                    'ell': ell 
                    }
            cat_corr = {
                    'catalog': {'name': catdict['name'], 'n_mock': 1}, 
                    'correction': corrdict, 
                    'spec': specdict
                    }
            avg_spec = AvgSpec(n_mocks[i_cat], 'pk', cat_corr)
            avg_spec.read(rebin=rebin)
        
            k_arr = avg_spec.k
            avg_pk = getattr(avg_spec, ''.join(['p', str(ell), 'k']))
            if n_mocks[i_cat] > 1: 
                pk_err = avg_spec.stddev(rebin=rebin)
                pk_err_lower = pk_err 
                pk_err_lower[np.abs(pk_err) > np.abs(avg_pk)] = avg_pk[np.abs(pk_err) > np.abs(avg_pk)] * 0.99

                if ell == 0: 
                    mono.fill_between(k_arr, avg_pk - pk_err_lower, avg_pk + pk_err, 
                         color='none', linewidth=1, edgecolor=cat_color, hatch=cat_hatch)
                    mono.plot(k_arr, np.zeros(len(avg_pk)), lw=4, 
                            color=cat_color, label=None)
                elif ell == 2: 
                    quad.fill_between(k_arr, np.abs(avg_pk - pk_err_lower), np.abs(avg_pk + pk_err),
                         color='none', linewidth=1, edgecolor=cat_color, hatch=cat_hatch)
                    quad.plot(k_arr, np.zeros(len(avg_pk)), lw=4, color=cat_color, label=cat_label)
            else: 
                if ell == 0: 
                    #mono.scatter(k_arr, avg_pk, s=2, color=cat_color, label=None)
                    mono.plot(k_arr, avg_pk, lw=2, color=cat_color, label=None)
                    pk_err = None 
                elif ell == 2: 
                    #quad.scatter(k_arr, np.abs(avg_pk), s=2, color=cat_color, label=cat_label)
                    quad.plot(k_arr, np.abs(avg_pk), lw=2, color=cat_color, label=cat_label)
                    pk_err = None

            fig_data.append({
                        'catalog': cat, 
                        'cat_corr_dict': cat_corr, 
                        'k': k_arr, 
                        'p'+str(ell)+'k': avg_pk, 
                        'p'+str(ell)+'k_err': pk_err}
                        )
    
    mono.set_xlabel('k (h/Mpc)', fontsize=24)
    mono.set_xlim([10**-3,10**0])
    mono.set_xscale('log')
    mono.set_ylabel(r'$\mathtt{P_0(k)}$', fontsize=24)
    mono.set_ylim([10**2,10**5.5])
    mono.set_yscale('log')
    
    quad.set_xlabel('k (h/Mpc)', fontsize=24)
    quad.set_xlim([10**-3,10**0])
    quad.set_xscale('log')

    quad.yaxis.tick_right()
    quad.yaxis.set_label_position('right')
    quad.set_ylabel(r'$\mathtt{|P_2(k)|}$', fontsize=24)
    quad.set_ylim([10**2,10**5.5])
    quad.set_yscale('log')
    quad.legend(loc='upper right', scatterpoints=1, prop={'size':20}) 
    
    plt.sca(mono)
    plt.xticks([10**-3, 10**-2, 10**-1, 1])
    plt.sca(quad)
    plt.xticks([10**-2, 10**-1, 1])
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(rebin)+'x'
    fig_file = 'figure/fc_paper/mock_catalog_Plk'+rebin_str+'.png'
    fig_data_file = 'figure/fc_paper/mock_catalog_Plk'+rebin_str+'.figdata.p'
    pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()

def fig_NN_p4k_mocks_resid(catalogs=['nseries'], n_mocks=[84], Ngrid=960, rebin=None): 
    ''' Figure that compares the fiber collision residuals Del P_l(k) 
    of the NN upweighting method for l = 0, 2 for mock catalogs  

    Parameters
    ----------
    catalogs : 
        List of strings that specify the mock catalogs to include. 
        (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    n_mocks : 
        
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14,4), dpi=1000)
    all_fig = fig.add_subplot(111)

    if len(catalogs) < 3: 
        for i in xrange(3 - len(catalogs)): 
            catalogs.append(catalogs[0])
            n_mocks.append(n_mocks[0])
    if len(catalogs) != len(n_mocks): 
        raise ValueError
    
    fig_data = [] 
    for i_ell, ell in enumerate([4]):  # monopole or quadrupole
        mock_list, mock_label_list = [], [] 
        for i_cat, cat in enumerate(catalogs):
            sub = fig.add_subplot(1, 3, 3*i_ell + i_cat+1)

            if cat == 'qpm': 
                cat_label = qpm_label 
                cat_color = qpm_color
            elif cat == 'bigmd': 
                cat_label = bmd_label
                cat_color = bmd_color
            elif cat == 'nseries': 
                cat_label = nsr_label 
                cat_color = nsr_color 
            else: 
                raise NotImplementedError 

            tru_catdict = {'name': cat, 'n_mock': 1}
            tru_corrdict = {'name': 'true'}
            upw_corrdict = {'name': 'upweight'}
            tru_specdict = {
                    'P0': 20000,
                    'Lbox': 3600, 
                    'Ngrid': Ngrid[i_ell], 
                    'ell': ell 
                    }
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
            tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
            tru_avg_spec.read(rebin=rebin)
            
            upw_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', upw_cat_corr)
            upw_avg_spec.read(rebin=rebin)

            now_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', now_cat_corr)
            now_avg_spec.read(rebin=rebin)
            now_avg_pk = getattr(now_avg_spec, ''.join(['p', str(ell), 'k']))
        
            k_arr = tru_avg_spec.k
            tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
            upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
            
            if n_mocks[i_cat] > 1: 
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

            mock = sub.scatter(k_arr, upw_avg_pk - tru_avg_pk, lw=0, color=cat_color)
            # P(k) with no weights 
            #noweight_plot, = sub.plot(k_arr, now_avg_pk - tru_avg_pk, 
            #        lw=1.2, ls='--', color='k', dashes=[10,5])
            noweight_plot = sub.scatter(k_arr, now_avg_pk - tru_avg_pk, 
                    marker='x', s=16, lw=0.75, color='k') # , edgecolor='k') #lw=0.3, facecolors='none', 
            noweight_label = r"$\Delta \mathtt{P_l^{NoW}(k)}$"

            if ell == 2: 
                mock_list.append(mock)
                mock_label_list.append(plt_label)

            # Axes 
            sub.set_xlim([10**-2,10**0])
            sub.set_xscale('log')

            if ell == 0: 
                sub.set_ylim([-750.,750.])
                sub.set_xticklabels([])
                plt.sca(sub)
                #plt.yticks(list(np.arange(-1000., 1250., 500.)) )
                plt.yticks(list(np.arange(-750., 1000., 250.)) )
                if i_cat == len(catalogs)-1: 
                    cv_label = [r"$\sigma_\mathtt{l} \mathtt{(k)}$"]
                    sub.legend([samp_var], cv_label, loc='upper right', scatterpoints=1, 
                            prop={'size':25}, handletextpad=0.5, borderpad=0.5) 
                elif i_cat == 0: 
                    sub.text(1.5*10**-2, 500., r"$\mathtt{l = 0}$", fontsize=24)
                else: 
                    lgd = sub.legend(
                            [mock, noweight_plot], 
                            [r"$\mathtt{\overline{P_l^{NN}} - \overline{P_l^{true}}}$", noweight_label], 
                            loc='upper right', scatterpoints=1, 
                            prop={'size':20}, borderpad=0.5, handletextpad=-0.25, 
                            markerscale=3., scatteryoffsets=[0.75,0.75]) 
                    lgd.legendHandles[0].set_color('black')
                    # k_chi^2
                    #sub.annotate(r'$\mathtt{k}^\mathtt{NN}_{\chi^2}$', 
                    #        xy=(0.14, 0.0),
                    #        xytext=(0.14, 250), 
                    #        arrowprops=dict(facecolor='black', shrink=0.05)
                    #        )
            elif ell == 4: 
                #sub.set_ylim([-0.1,1.2])
                sub.set_ylim([-1000.,1000.])
                if i_cat == len(catalogs)-1: 
                    pass
                    #sub.legend([mock_list[i_cat]], [mock_label_list[i_cat]], 
                    #        loc='lower right', scatterpoints=1, 
                    #        prop={'size':20}, borderpad=0.1, handletextpad=0.5) 
                elif i_cat == 0: 
                    sub.text(1.5*10**-2, -750., r"$\mathtt{l = 2}$", fontsize=24)
                    plt.sca(sub)
                    plt.yticks(list(np.arange(-1000., 1000., 500.)) )
                    plt.xticks([10**-2, 10**-1])
                    #sub.legend([mock_list[i_cat]], [mock_label_list[i_cat]], 
                    #        loc='lower right', scatterpoints=1, 
                    #        prop={'size':20}, borderpad=0.1, handletextpad=0.5) 
                else: 
                    sub.set_xlabel('k (h/Mpc)', fontsize=24)
                    plt.sca(sub)
                    plt.xticks([10**-2, 10**-1])
                    #sub.legend([mock_list[i_cat]], [mock_label_list[i_cat]], loc='lower right', scatterpoints=1, 
                    #        prop={'size':20}, borderpad=0.1, handletextpad=0.5) 
                #sub.text(0.8, -850., mock_label_list[i_cat], fontsize=20, horizontalalignment='right') 
                    
                    # k_chi^2
                    #sub.annotate(r'$\mathtt{k}^\mathtt{NN}_{\chi^2}$', 
                    #        xy=(0.30, 0.0),
                    #        xytext=(0.30, -400), 
                    #        arrowprops=dict(facecolor='black', shrink=0.05)
                    #        )
            #if i_cat == 0:
            #    sub.set_ylabel(
            #            r"$\mathtt{\overline{P_{"+str(ell)+"}^{NN}} - \overline{P_{"+str(ell)+"}^{true}}(k)}$", 
            #            fontsize=24)
            if i_cat != 0 : 
                sub.set_yticklabels([])

            fig_data.append({
                'ell': ell,
                'name': cat , 
                'k': k_arr, 
                'true_Pk': tru_avg_pk, 
                'corr_Pk': upw_avg_pk,
                'resid': upw_avg_pk - tru_avg_pk,
                'cv': pk_err                
                })

    #yax = fig.get_yaxis()
    all_fig.set_xticklabels([])
    all_fig.set_yticklabels([])
    all_fig.set_ylabel(
            r"$\Delta \mathtt{P_l(k)}$", 
            fontsize=24, labelpad=50)
    #r"$\mathtt{\overline{P_l^{NN}} - \overline{P_l^{true}}(k)}$", 
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    
    if rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(rebin)+'x'
    fig_file = ''.join(['figure/fc_paper/',
        'mock_catalog_NN_true_P4k_resid', rebin_str, 
        '.png'])
    fig_data_file = 'figure/fc_paper/mock_catalog_NN_true_P4k_resid'+rebin_str+'.figdata.p' 
    pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()


# NN Method Residual plot
def fig_NN_p02k_mocks_resid(catalogs=['nseries'], n_mocks=[84], Ngrid=960, rebin=None): 
    ''' Figure that compares the fiber collision residuals Del P_l(k) 
    of the NN upweighting method for l = 0, 2 for mock catalogs  

    Parameters
    ----------
    catalogs : 
        List of strings that specify the mock catalogs to include. 
        (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    n_mocks : 
        
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14,8), dpi=1000)
    all_fig = fig.add_subplot(111)

    if len(catalogs) < 3: 
        for i in xrange(3 - len(catalogs)): 
            catalogs.append(catalogs[0])
            n_mocks.append(n_mocks[0])
    if len(catalogs) != len(n_mocks): 
        raise ValueError
    
    fig_data = [] 
    for i_ell, ell in enumerate([0, 2]):  # monopole or quadrupole
        mock_list, mock_label_list = [], [] 
        for i_cat, cat in enumerate(catalogs):
            sub = fig.add_subplot(2, 3, 3*i_ell + i_cat+1)

            if cat == 'qpm': 
                cat_label = qpm_label 
                cat_color = qpm_color
            elif cat == 'bigmd': 
                cat_label = bmd_label
                cat_color = bmd_color
            elif cat == 'nseries': 
                cat_label = nsr_label 
                cat_color = nsr_color 
            else: 
                raise NotImplementedError 

            tru_catdict = {'name': cat, 'n_mock': 1}
            tru_corrdict = {'name': 'true'}
            upw_corrdict = {'name': 'upweight'}
            tru_specdict = {
                    'P0': 20000,
                    'Lbox': 3600, 
                    'Ngrid': Ngrid[i_ell], 
                    'ell': ell 
                    }
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
            tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
            tru_avg_spec.read(rebin=rebin)
            
            upw_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', upw_cat_corr)
            upw_avg_spec.read(rebin=rebin)

            now_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', now_cat_corr)
            now_avg_spec.read(rebin=rebin)
            now_avg_pk = getattr(now_avg_spec, ''.join(['p', str(ell), 'k']))
        
            k_arr = tru_avg_spec.k
            tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
            upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
            
            if n_mocks[i_cat] > 1: 
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

            mock = sub.scatter(k_arr, upw_avg_pk - tru_avg_pk, lw=0, color=cat_color)
            # P(k) with no weights 
            #noweight_plot, = sub.plot(k_arr, now_avg_pk - tru_avg_pk, 
            #        lw=1.2, ls='--', color='k', dashes=[10,5])
            noweight_plot = sub.scatter(k_arr, now_avg_pk - tru_avg_pk, 
                    marker='x', s=16, lw=0.75, color='k') # , edgecolor='k') #lw=0.3, facecolors='none', 
            noweight_label = r"$\Delta \mathtt{P_l^{NoW}(k)}$"

            if ell == 2: 
                mock_list.append(mock)
                mock_label_list.append(plt_label)

            # Axes 
            sub.set_xlim([10**-2,10**0])
            sub.set_xscale('log')

            if ell == 0: 
                sub.set_ylim([-750.,750.])
                sub.set_xticklabels([])
                plt.sca(sub)
                #plt.yticks(list(np.arange(-1000., 1250., 500.)) )
                plt.yticks(list(np.arange(-750., 1000., 250.)) )
                if i_cat == len(catalogs)-1: 
                    cv_label = [r"$\sigma_\mathtt{l} \mathtt{(k)}$"]
                    sub.legend([samp_var], cv_label, loc='upper right', scatterpoints=1, 
                            prop={'size':25}, handletextpad=0.5, borderpad=0.5) 
                elif i_cat == 0: 
                    sub.text(1.5*10**-2, 500., r"$\mathtt{l = 0}$", fontsize=24)
                else: 
                    lgd = sub.legend(
                            [mock, noweight_plot], 
                            [r"$\mathtt{\overline{P_l^{NN}} - \overline{P_l^{true}}}$", noweight_label], 
                            loc='upper right', scatterpoints=1, 
                            prop={'size':20}, borderpad=0.5, handletextpad=-0.25, 
                            markerscale=3., scatteryoffsets=[0.75,0.75]) 
                    lgd.legendHandles[0].set_color('black')
                    # k_chi^2
                    #sub.annotate(r'$\mathtt{k}^\mathtt{NN}_{\chi^2}$', 
                    #        xy=(0.14, 0.0),
                    #        xytext=(0.14, 250), 
                    #        arrowprops=dict(facecolor='black', shrink=0.05)
                    #        )
            elif ell == 2: 
                #sub.set_ylim([-0.1,1.2])
                sub.set_ylim([-1000.,1000.])
                if i_cat == len(catalogs)-1: 
                    pass
                    #sub.legend([mock_list[i_cat]], [mock_label_list[i_cat]], 
                    #        loc='lower right', scatterpoints=1, 
                    #        prop={'size':20}, borderpad=0.1, handletextpad=0.5) 
                elif i_cat == 0: 
                    sub.text(1.5*10**-2, -750., r"$\mathtt{l = 2}$", fontsize=24)
                    plt.sca(sub)
                    plt.yticks(list(np.arange(-1000., 1000., 500.)) )
                    plt.xticks([10**-2, 10**-1])
                    #sub.legend([mock_list[i_cat]], [mock_label_list[i_cat]], 
                    #        loc='lower right', scatterpoints=1, 
                    #        prop={'size':20}, borderpad=0.1, handletextpad=0.5) 
                else: 
                    sub.set_xlabel('k (h/Mpc)', fontsize=24)
                    plt.sca(sub)
                    plt.xticks([10**-2, 10**-1])
                    #sub.legend([mock_list[i_cat]], [mock_label_list[i_cat]], loc='lower right', scatterpoints=1, 
                    #        prop={'size':20}, borderpad=0.1, handletextpad=0.5) 
                sub.text(0.8, -850., mock_label_list[i_cat], fontsize=20, horizontalalignment='right') 
                    
                    # k_chi^2
                    #sub.annotate(r'$\mathtt{k}^\mathtt{NN}_{\chi^2}$', 
                    #        xy=(0.30, 0.0),
                    #        xytext=(0.30, -400), 
                    #        arrowprops=dict(facecolor='black', shrink=0.05)
                    #        )
            #if i_cat == 0:
            #    sub.set_ylabel(
            #            r"$\mathtt{\overline{P_{"+str(ell)+"}^{NN}} - \overline{P_{"+str(ell)+"}^{true}}(k)}$", 
            #            fontsize=24)
            if i_cat != 0 : 
                sub.set_yticklabels([])

            fig_data.append({
                'ell': ell,
                'name': cat , 
                'k': k_arr, 
                'true_Pk': tru_avg_pk, 
                'corr_Pk': upw_avg_pk,
                'resid': upw_avg_pk - tru_avg_pk,
                'cv': pk_err                
                })

    #yax = fig.get_yaxis()
    all_fig.set_xticklabels([])
    all_fig.set_yticklabels([])
    all_fig.set_ylabel(
            r"$\Delta \mathtt{P_l(k)}$", 
            fontsize=24, labelpad=50)
    #r"$\mathtt{\overline{P_l^{NN}} - \overline{P_l^{true}}(k)}$", 
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    
    if rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(rebin)+'x'
    fig_file = ''.join(['figure/fc_paper/',
        'mock_catalog_NN_true_Plk_resid', rebin_str, 
        '.png'])
    fig_data_file = 'figure/fc_paper/mock_catalog_NN_true_Plk_resid'+rebin_str+'.figdata.p' 
    pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()

def fig_NN_p024k_mocks_resid(catalogs=['nseries'], n_mocks=[84], Ngrid=960, rebin=None): 
    ''' Figure that compares the fiber collision residuals Del P_l(k) 
    of the NN upweighting method for l = 0, 2 for mock catalogs  

    Parameters
    ----------
    catalogs : 
        List of strings that specify the mock catalogs to include. 
        (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    n_mocks : 
        
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14,12), dpi=1000)
    all_fig = fig.add_subplot(111)

    if len(catalogs) < 3: 
        for i in xrange(3 - len(catalogs)): 
            catalogs.append(catalogs[0])
            n_mocks.append(n_mocks[0])
    if len(catalogs) != len(n_mocks): 
        raise ValueError
    
    fig_data = [] 
    for i_ell, ell in enumerate([0, 2, 4]):  # monopole or quadrupole
        mock_list, mock_label_list = [], [] 
        for i_cat, cat in enumerate(catalogs):
            sub = fig.add_subplot(3, 3, 3*i_ell + i_cat+1)

            if cat == 'qpm': 
                cat_label = qpm_label 
                cat_color = qpm_color
            elif cat == 'bigmd': 
                cat_label = bmd_label
                cat_color = bmd_color
            elif cat == 'nseries': 
                cat_label = nsr_label 
                cat_color = nsr_color 
            else: 
                raise NotImplementedError 

            tru_catdict = {'name': cat, 'n_mock': 1}
            tru_corrdict = {'name': 'true'}
            upw_corrdict = {'name': 'upweight'}
            tru_specdict = {
                    'P0': 20000,
                    'Lbox': 3600, 
                    'Ngrid': Ngrid[i_ell], 
                    'ell': ell 
                    }
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
            tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
            tru_avg_spec.read(rebin=rebin)
            
            upw_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', upw_cat_corr)
            upw_avg_spec.read(rebin=rebin)

            now_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', now_cat_corr)
            now_avg_spec.read(rebin=rebin)
            now_avg_pk = getattr(now_avg_spec, ''.join(['p', str(ell), 'k']))
        
            k_arr = tru_avg_spec.k
            tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
            upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
            
            if n_mocks[i_cat] > 1: 
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
            elif ell == 4: 
                plt_label = cat_label

            mock = sub.scatter(k_arr, upw_avg_pk - tru_avg_pk, lw=0, color=cat_color)
            # P(k) with no weights 
            #noweight_plot, = sub.plot(k_arr, now_avg_pk - tru_avg_pk, 
            #        lw=1.2, ls='--', color='k', dashes=[10,5])
            noweight_plot = sub.scatter(k_arr, now_avg_pk - tru_avg_pk, 
                    marker='x', s=16, lw=0.75, color='k') # , edgecolor='k') #lw=0.3, facecolors='none', 
            noweight_label = r"$\Delta \mathtt{P_l^{NoW}(k)}$"

            if ell == 4: 
                mock_list.append(mock)
                mock_label_list.append(plt_label)

            # Axes 
            sub.set_xlim([10**-2,10**0])
            sub.set_xscale('log')

            if ell == 0: 
                sub.set_ylim([-750.,750.])
                sub.set_xticklabels([])
                plt.sca(sub)
                plt.yticks(list(np.arange(-750., 1000., 250.)) )
                if i_cat == len(catalogs)-1: 
                    cv_label = [r"$\sigma_\mathtt{l} \mathtt{(k)}$"]
                    sub.legend([samp_var], cv_label, loc='upper right', scatterpoints=1, 
                            prop={'size':25}, handletextpad=0.5, borderpad=0.5) 
                elif i_cat == 0: 
                    sub.text(1.5*10**-2, 500., r"$\ell\mathtt{= 0}$", fontsize=24)
                else: 
                    lgd = sub.legend(
                            [mock, noweight_plot], 
                            [r"$\mathtt{\overline{P_l^{NN}} - \overline{P_l^{true}}}$", noweight_label], 
                            loc='upper right', scatterpoints=1, 
                            prop={'size':20}, borderpad=0.5, handletextpad=-0.25, 
                            markerscale=3., scatteryoffsets=[0.75,0.75]) 
                    lgd.legendHandles[0].set_color('black')
            elif ell == 2: 
                sub.set_ylim([-1000.,1000.])
                if i_cat == len(catalogs)-1: 
                    pass
                elif i_cat == 0: 
                    sub.text(1.5*10**-2, -750., r"$\ell\mathtt{= 2}$", fontsize=24)
                    plt.yticks(list(np.arange(-500., 1000., 500.)) )
                sub.set_xticklabels([])
            elif ell == 4: 
                sub.set_ylim([-500.,500.])
                if i_cat == len(catalogs)-1: 
                    pass
                elif i_cat == 0: 
                    sub.text(1.5*10**-2, -400., r"$\ell\mathtt{= 4}$", fontsize=24)
                    plt.sca(sub)
                    #plt.yticks(list(np.arange(-500., 1000., 500.)) )
                    plt.xticks([10**-2, 10**-1])
                else: 
                    sub.set_xlabel('k (h/Mpc)', fontsize=24)
                    plt.sca(sub)
                    plt.xticks([10**-2, 10**-1])
                sub.text(0.8, -400., mock_label_list[i_cat], fontsize=20, horizontalalignment='right') 
                    
            if i_cat != 0 : 
                sub.set_yticklabels([])

            fig_data.append({
                'ell': ell,
                'name': cat , 
                'k': k_arr, 
                'true_Pk': tru_avg_pk, 
                'corr_Pk': upw_avg_pk,
                'resid': upw_avg_pk - tru_avg_pk,
                'cv': pk_err                
                })

    all_fig.set_xticklabels([])
    all_fig.set_yticklabels([])
    all_fig.set_ylabel(
            r"$\Delta \mathtt{P_l(k)}$", 
            fontsize=24, labelpad=50)
    #r"$\mathtt{\overline{P_l^{NN}} - \overline{P_l^{true}}(k)}$", 
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    
    if rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(rebin)+'x'
    fig_file = ''.join(['figure/fc_paper/',
        'mock_catalog_NN_true_P024k_resid', rebin_str, 
        '.png'])
    fig_data_file = 'figure/fc_paper/mock_catalog_NN_true_P024k_resid'+rebin_str+'.figdata.p' 
    pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()


# NN Method Normalized residual plot
def fig_NN_p0k_mock_norm_resid(catalogs=['nseries'], n_mocks=[84], Ngrid=960, rebin=None, chi2_use='nseries'): 
    ''' Figure that compares the normalized power spectrum multipole residuals 
    of the NN method for only l = 0.

    Parameters
    ----------
    catalogs : 
        List of strings that specify the mock catalogs to include. 
        (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    n_mocks : 
        
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(8,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    sub = plt.subplot(gs[0])#fig.add_subplot(121)
    chi_plot = plt.subplot(gs[1]) # fig.add_subplot(122)

    if len(catalogs) < 3: 
        for i in xrange(3 - len(catalogs)): 
            catalogs.append(catalogs[0])
            n_mocks.append(n_mocks[0])
    if len(catalogs) != len(n_mocks): 
        raise ValueError
    
    mock_list, mock_label_list = [], [] 
    for i_cat, cat in enumerate(catalogs):
        if cat == 'qpm': 
            cat_label = qpm_label 
            cat_color = qpm_color
        elif cat == 'bigmd': 
            cat_label = bmd_label
            cat_color = bmd_color
        elif cat == 'nseries': 
            cat_label = nsr_label 
            cat_color = nsr_color 
        else: 
            raise NotImplementedError 

        tru_catdict = {'name': cat, 'n_mock': 1}
        tru_corrdict = {'name': 'true'}
        upw_corrdict = {'name': 'upweight'}
        tru_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': 0}
        tru_cat_corr = {
                'catalog': tru_catdict, 
                'correction': tru_corrdict, 
                'spec': tru_specdict
                }
        upw_cat_corr = {
                'catalog': tru_catdict, 
                'correction': upw_corrdict, 
                'spec': tru_specdict
                }
        tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
        tru_avg_spec.read(rebin=rebin)
        
        upw_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', upw_cat_corr)
        upw_avg_spec.read(rebin=rebin)
        
        k_arr = tru_avg_spec.k
        tru_avg_pk = getattr(tru_avg_spec, ''.join(['p0k']))
        upw_avg_pk = getattr(upw_avg_spec, ''.join(['p0k']))

        if n_mocks[i_cat] > 1: 
            pk_err = tru_avg_spec.stddev(rebin=rebin)
            if cat == 'nseries':  
                #nseries_cv, = sub.plot(k_arr, np.abs(pk_err/tru_avg_pk), 
                #        lw=1, color='gray')
                nseries_cv = sub.fill_between(
                        k_arr, 
                        np.zeros(len(k_arr)), 
                        np.abs(pk_err/tru_avg_pk),
                        lw=0, 
                        color='grey',
                        alpha=0.6
                        )
        if cat == 'bigmd': 
            bigmd = sub.scatter(k_arr, np.abs(1.- upw_avg_pk/tru_avg_pk), s=6, lw=0, color=cat_color)
            #bigmd, = sub.plot(k_arr, np.abs(1.- upw_avg_pk/tru_avg_pk), lw=3, color=cat_color)
        elif cat == 'qpm':
            qpm, = sub.plot(k_arr, np.abs(1.- upw_avg_pk/tru_avg_pk), lw=3, ls='--', color=cat_color)
        elif cat == 'nseries':
            nseries, = sub.plot(k_arr, np.abs(1.- upw_avg_pk/tru_avg_pk), lw=3, color=cat_color)

        sub.plot((10**-3, 10**0), (0.0,0.0), 'k', lw=2, ls='--')
 
        for ell in [0, 2]:
            k_chi = chisquared_p02k_mock(ell, corr='NN', catalogs=[cat], n_mocks=[n_mocks[i_cat]], 
                    Ngrid=Ngrid, rebin='log', use=chi2_use)
            print ell, cat, k_chi
            if ell == 0: 
                if cat == 'nseries': 
                    ypos = 12.
                elif cat == 'qpm': 
                    ypos = 8.
                elif cat == 'bigmd': 
                    ypos = 4.
            elif ell == 2: 
                if cat == 'nseries': 
                    ypos = -12.
                elif cat == 'qpm': 
                    ypos = -4.
                elif cat == 'bigmd': 
                    ypos = -8.
            chi_plot.errorbar(np.repeat(k_chi,2), np.repeat(0.5 * ypos, 2), 
                    xerr=0.15 * k_chi, xuplims=True, color=cat_color, 
                    capsize=4., capthick=2, elinewidth=2) 

            #chi_plot.annotate('', #r'$\mathtt{k}^\mathtt{NN}_{\chi^2}$', 
            #        xy=(k_chi, 0.0),
            #        xytext=(k_chi, ypos), 
            #        arrowprops=dict(facecolor=cat_color, shrink=0.05, linewidth=0., alpha=0.8)#hatch='X')#linewidth=0.)
            #        )

    for ell in [0, 2]: 
        if ell == 0: 
            ypos = 6.
        elif ell == 2: 
            ypos = -14.
        chi_plot.text(1.2*10**-2, ypos, 
                r'$\mathtt{k}^\mathtt{NN; l='+str(ell)+'}_{\chi^2}$', fontsize=24)
    chi_plot.plot(np.arange(10**-2, 10**0, 0.01), np.zeros(len(np.arange(10**-2, 10**0, 0.01))), c='k', lw=2, ls='--') 
    # x-axis 
    sub.set_xlim([10**-2,10**0])
    sub.set_xscale('log')
    sub.set_xticklabels([])

    chi_plot.set_xlim([10**-2,10**0])
    chi_plot.set_xscale('log') 
    chi_plot.set_xlabel('k (h/Mpc)', fontsize=24)
    # y-axis    
    sub.set_ylim([-0.01,0.30])
    chi_plot.set_ylim([-20, 20]) 
    chi_plot.set_yticklabels([])
    chi_plot.set_yticks([])

    cv_labels = [r"$\sigma_\mathtt{l}\mathtt{(k)}/\mathtt{P_l(k)}$"]
    sub.legend([bigmd, nseries, qpm]+[nseries_cv], [bmd_label, nsr_label, qpm_label]+cv_labels, loc='upper left', scatterpoints=1, 
            prop={'size':20}, handletextpad=0.4, borderpad=1.)#, ncol=len(mock_list)-1)

    sub.set_ylabel(
            r"$\mathtt{1- \overline{P_0^{NN}}/\overline{P_0^{true}}(k)}$", 
            fontsize=24)

    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(rebin)+'x'
    fig_file = ''.join(['figure/fc_paper/', 
        'mock_catalog_NN_true_P0k_norm_resid', rebin_str, '.png'])
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()

# LOS Reconstruction Residual plot
def fig_dlospeak_p02k_mock_resid(catalogs=['nseries'], n_mocks=[84], Ngrid=[960], florian_offset=None, rebin=None): 
    ''' Figure that compares the P_l(k) residual of the line-of-sight reconstructed 
    power spectrum monopole and quadrupole of specified mock catalogs. 
    
    P_l(k)^dlospeak - P_l(k)^true l = 0, 2 for 

    Parameters
    ----------
    catalogs : list
        List of strings that specify the mock catalogs to include. 
        (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])

    n_mocks : list 
        List of ints of mocks to average the P_l(k) over. Same length as 
        the list of catalog names. 

    Ngrid : list 
        List of ints that specify the Ngrids of the power spectrum 
        monopole and quadrupole. Ngrid list has two elements. 
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14,8))
    all_fig = fig.add_subplot(111)

    if len(catalogs) < 3: 
        for i in xrange(3 - len(catalogs)): 
            catalogs.append(catalogs[0])
            n_mocks.append(n_mocks[0])
    if len(catalogs) != len(n_mocks): 
        raise ValueError

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

    fig_data = [] 
    for i_ell, ell in enumerate([0, 2]):  # monopole or quadrupole
        mock_list, mock_label_list = [], [] 
        for i_cat, cat in enumerate(catalogs):
            sub = fig.add_subplot(2, 3, 3*i_ell + i_cat+1)

            # dLOS peak correction parameters hardcoded for 
            # each of the mock catalogs
            if cat == 'qpm': 
                cat_label = qpm_label 
                cat_color = qpm_color
            elif cat == 'bigmd': 
                cat_label = bmd_label 
                cat_color = bmd_color
            elif cat == 'nseries': 
                cat_label = nsr_label
                cat_color = nsr_color
            else: 
                raise NotImplementedError 
            dlospeak_corrdict = dlospeak_dict[cat]  

            tru_catdict = {'name': cat, 'n_mock': 1}
            tru_corrdict = {'name': 'true'}
            tru_specdict = {
                    'P0': 20000,
                    'Lbox': 3600, 
                    'Ngrid': Ngrid[i_ell], 
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

            tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
            tru_avg_spec.read(rebin=rebin)
            
            dlospeak_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', dlospeak_cat_corr)
            dlospeak_avg_spec.read(rebin=rebin)
            
            k_arr = tru_avg_spec.k
            tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
            dlospeak_avg_pk = getattr(dlospeak_avg_spec, ''.join(['p', str(ell), 'k']))

            if n_mocks[i_cat] > 1: 
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

            mock = sub.scatter(k_arr, dlospeak_avg_pk - tru_avg_pk, lw=0, color=cat_color)
            if ell == 0: 
                florian_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', florian_cat_corr)
                florian_avg_spec.read(rebin=rebin)
                florian_avg_pk = getattr(florian_avg_spec, ''.join(['p', str(ell), 'k']))
                #florian, = sub.plot(k_arr, florian_avg_pk - tru_avg_pk, 
                #        lw=1., ls='--', color='k', dashes=[10,7])
                if florian_offset is None: 
                    florian_offset = 0 

                #florian = sub.scatter(k_arr, florian_avg_pk - tru_avg_pk - florian_offset, 
                #        marker='+', s=20, lw=0.7, color='k') # , edgecolor='k') #lw=0.3, facecolors='none', 

                hector_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', hector_cat_corr)
                hector_avg_spec.read(rebin=rebin)
                hector_avg_pk = getattr(hector_avg_spec, ''.join(['p', str(ell), 'k']))
                #hector = sub.scatter(k_arr, hector_avg_pk - tru_avg_pk, 
                #        s=8, lw=0.5, facecolors='none', edgecolor='k')#, dashes=[10,5])
                hector, = sub.plot(k_arr, hector_avg_pk - tru_avg_pk, 
                        lw=1., ls='--', color='k', dashes=[10,7])
            
            if ell == 2: 
                mock_list.append(mock)
                mock_label_list.append(plt_label)

            # Axes 
            sub.set_xlim([10**-2,10**0])
            sub.set_xscale('log')

            if ell == 0: 
                sub.set_ylim([-750.,750.])
                sub.set_xticklabels([])
                plt.sca(sub)
                #plt.yticks(list(np.arange(-1000., 1250., 500.)) )
                plt.yticks(list(np.arange(-750., 1000., 250.)) )
                if i_cat == len(catalogs)-1: 
                    cv_label = [r"$\sigma_\mathtt{l} \mathtt{(k)}$"]
                    sub.legend([samp_var], cv_label, loc='upper right', scatterpoints=1, 
                            prop={'size':25}, handletextpad=0.5, borderpad=0.5) 
                elif i_cat == 0: 
                    sub.text(1.5*10**-2, 500., r"$\ell = 0$", fontsize=24)
                else: 
                    leg = sub.legend([mock], ['$\mathtt{P_l^{LRec}-P_l^{true}}$'], 
                            loc='upper right', scatterpoints=1, 
                            prop={'size':20}, borderpad=0.5, handletextpad=-0.25, 
                            markerscale=3., scatteryoffsets=[0.75]) 
                    leg.legendHandles[0].set_color('black')
            elif ell == 2: 
                #sub.set_ylim([-0.1,1.2])
                sub.set_ylim([-1000.,1000.])
                if i_cat == len(catalogs)-1: 
                    #sub.legend([florian, hector], ['Beutler et al.(2014)', r'Gil-Mar\'{i}n et al.(2014)'], 
                    #        loc='upper right', scatterpoints=1, 
                    #        prop={'size':17}, handletextpad=0.0, borderpad=0.5, 
                    #        markerscale=3, scatteryoffsets=[0.5, 0.0])
                    sub.legend([hector], [r'Gil-Mar\'{i}n et al.(2014)'], 
                            loc='upper right', scatterpoints=1, 
                            prop={'size':17}, handletextpad=0.0, borderpad=0.5, 
                            markerscale=3, scatteryoffsets=[0.0])
                elif i_cat == 0: 
                    sub.text(1.5*10**-2, 650., r"$\ell = 2$", fontsize=24)
                    plt.sca(sub)
                    plt.yticks(list(np.arange(-1000., 1000., 500.)) )
                    plt.xticks([10**-2, 10**-1])
                    #sub.legend([mock_list[i_cat]], [mock_label_list[i_cat]], loc='lower right', scatterpoints=1, 
                    #        prop={'size':20}, borderpad=0.1, handletextpad=0.5) 
                else: 
                    sub.set_xlabel('k (h/Mpc)', fontsize=24)
                    plt.sca(sub)
                    plt.xticks([10**-2, 10**-1])
                    #sub.legend([mock_list[i_cat]], [mock_label_list[i_cat]], loc='lower right', scatterpoints=1, 
                    #        prop={'size':20}, borderpad=0.1, handletextpad=0.5) 
                    # k_chi^2                    
                    #sub.annotate(r'$\mathtt{k}^\mathtt{LRec}_{\chi^2}$', 
                    #        xy=(0.32, 0.0),
                    #        xytext=(0.32, 350), 
                    #        arrowprops=dict(facecolor='black', shrink=0.05)
                    #        )
                sub.text(0.8, -850., mock_label_list[i_cat], fontsize=20, horizontalalignment='right') 

            if i_cat != 0 : 
                sub.set_yticklabels([])

            fig_data.append({
                'ell': ell,
                'name': cat , 
                'k': k_arr,
                'true_Pk': tru_avg_pk, 
                'corr_Pk': dlospeak_avg_pk, 
                'resid': dlospeak_avg_pk - tru_avg_pk,
                'cv': pk_err
                })

    all_fig.set_xticklabels([])
    all_fig.set_yticklabels([])
    all_fig.set_ylabel(r"$\Delta \mathtt{P_l(k)}$", fontsize=24, labelpad=50)

    #r"$\mathtt{\overline{P_l^{LRec}} - \overline{P_l^{true}}(k)}$", 
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    
    if florian_offset is None: 
        florian_str = ''
    else: 
        florian_str = '.floran_offset'+str(round(florian_offset,1))
    if rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(rebin)+'x'
    fig_file = ''.join(['figure/fc_paper/', 
        'mock_catalog_dlospeak_true_Plk_resid', florian_str, rebin_str, '.png' 
        ])
    fig_data_file = ''.join(['figure/fc_paper/', 
        'mock_catalog_dlospeak_true_Plk_resid', florian_str, rebin_str, '.figdata.p'
        ])
    pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()

# LOS Reconstruction Normalized Residual plot
def fig_dlospeak_p0k_mock_norm_resid(catalogs=['nseries'], n_mocks=[84], Ngrid=960, rebin=None, chi2_use='nseries'): 
    ''' Figure that compares the normalized power spectrum multipole residuals 
    of the NN method for only l = 0.

    Parameters
    ----------
    catalogs : 
        List of strings that specify the mock catalogs to include. 
        (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    n_mocks : 
        
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(8,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    sub = plt.subplot(gs[0])
    chi_plot = plt.subplot(gs[1])

    if len(catalogs) < 3: 
        for i in xrange(3 - len(catalogs)): 
            catalogs.append(catalogs[0])
            n_mocks.append(n_mocks[0])
    if len(catalogs) != len(n_mocks): 
        raise ValueError
    
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

    mock_list, mock_label_list = [], [] 
    for i_cat, cat in enumerate(catalogs):
        if cat == 'qpm': 
            cat_label = qpm_label 
            cat_color = qpm_color
        elif cat == 'bigmd': 
            cat_label = bmd_label
            cat_color = bmd_color
        elif cat == 'nseries': 
            cat_label = nsr_label 
            cat_color = nsr_color 
        else: 
            raise NotImplementedError 

        tru_catdict = {'name': cat, 'n_mock': 1}
        tru_corrdict = {'name': 'true'}
        upw_corrdict = {'name': 'upweight'}
        tru_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': 0}
        tru_cat_corr = {
                'catalog': tru_catdict, 
                'correction': tru_corrdict, 
                'spec': tru_specdict
                }
        dlospeak_cat_corr = {
                'catalog': tru_catdict, 
                'correction': dlospeak_dict[cat], 
                'spec': tru_specdict
                }
        tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
        tru_avg_spec.read(rebin=rebin)
        
        dlospeak_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', dlospeak_cat_corr)
        dlospeak_avg_spec.read(rebin=rebin)
        
        k_arr = tru_avg_spec.k
        tru_avg_pk = getattr(tru_avg_spec, ''.join(['p0k']))
        dlospeak_avg_pk = getattr(dlospeak_avg_spec, ''.join(['p0k']))

        if n_mocks[i_cat] > 1: 
            pk_err = tru_avg_spec.stddev(rebin=rebin)
            if cat == 'nseries':  
                nseries_cv = sub.fill_between(
                        k_arr, 
                        np.zeros(len(k_arr)), 
                        np.abs(pk_err/tru_avg_pk),
                        lw=0, 
                        color='grey',
                        alpha=0.6
                        )
        if cat == 'bigmd': 
            bigmd = sub.scatter(k_arr, np.abs(1.- dlospeak_avg_pk/tru_avg_pk), s=6, lw=0, color=cat_color)
        elif cat == 'qpm':
            qpm, = sub.plot(k_arr, np.abs(1.- dlospeak_avg_pk/tru_avg_pk), lw=3, ls='--', color=cat_color)
        elif cat == 'nseries':
            nseries, = sub.plot(k_arr, np.abs(1.- dlospeak_avg_pk/tru_avg_pk), lw=3, color=cat_color)

        sub.plot((10**-3, 10**0), (0.0,0.0), 'k', lw=2, ls='--')
        
        # chi-squared panel 
        for ell in [0, 2]:
            k_chi = chisquared_p02k_mock(ell, corr='LR', catalogs=[cat], n_mocks=[n_mocks[i_cat]], 
                    Ngrid=Ngrid, rebin='log', use=chi2_use)
            print ell, cat, k_chi
            if ell == 0: 
                if cat == 'nseries': 
                    ypos = 8.
                elif cat == 'qpm': 
                    ypos = 12.
                elif cat == 'bigmd': 
                    ypos = 4.
            elif ell == 2: 
                if cat == 'nseries': 
                    ypos = -8.
                elif cat == 'qpm': 
                    ypos = -8.
                elif cat == 'bigmd': 
                    ypos = -8.

            chi_plot.errorbar(np.repeat(k_chi,2), np.repeat(0.5 * ypos, 2), 
                    xerr=0.15 * k_chi, xuplims=True, color=cat_color, 
                    capsize=4., capthick=2, elinewidth=2) 

            #chi_plot.annotate('', #r'$\mathtt{k}^\mathtt{NN}_{\chi^2}$', 
            #        xy=(k_chi, 0.0),
            #        xytext=(k_chi, ypos), 
            #        arrowprops=dict(facecolor=cat_color, shrink=0.05, linewidth=0., alpha=0.8)#hatch='X')#linewidth=0.)
            #        )
    for ell in [0, 2]: 
        k_chi = chisquared_p02k_mock(ell, corr='NN', catalogs=catalogs, n_mocks=n_mocks, 
                Ngrid=Ngrid, rebin='log', use=chi2_use)

        if ell == 0: 
            factor = 8 
        elif ell == 2: 
            factor = -8 

        chi_plot.errorbar(np.repeat(k_chi,2), np.repeat(0.5 * factor, 2), 
                xerr=0.15 * k_chi, xuplims=True, color='k', 
                capsize=4., capthick=2, elinewidth=2) 
        if ell == 0: 
            chi_plot.text(0.85 * k_chi, 0.7 * factor, r'$\mathtt{NN}$') 
        elif ell == 2:  
            chi_plot.text(0.85 * k_chi, 1.1 * factor, r'$\mathtt{NN}$') 

        #chi_plot.annotate(r'$\mathtt{NN}$', #r'$\mathtt{k}^\mathtt{NN}_{\chi^2}$', 
        #        xy=(k_chi, 0.0),
        #        xytext=(k_chi*0.9, factor), 
        #        arrowprops=dict(facecolor='k', shrink=0.05, headwidth=8, width=3, linewidth=0., alpha=1.)
        #        )
        if ell == 0: 
            gm_k_chi = chisquared_p02k_mock(ell, corr='GM', catalogs=catalogs, n_mocks=n_mocks, 
                Ngrid=Ngrid, rebin='log', use=chi2_use)
            print 'Gil-Marin k_chi^2  = ', gm_k_chi

            chi_plot.errorbar(np.repeat(gm_k_chi,2), np.repeat(0.5 * factor, 2), 
                    xerr=0.15 * gm_k_chi, xuplims=True, color='gray', 
                    capsize=4., capthick=2, elinewidth=2) 
            chi_plot.text(0.8 * gm_k_chi, 0.7 * factor, r'$\mathtt{GM+}$') 

            #chi_plot.annotate(r'$\mathtt{GM+}$',  
            #        xy=(gm_k_chi, 0.0),
            #        xytext=(gm_k_chi*0.85, factor), 
            #        arrowprops=dict(facecolor='k', shrink=0.05, headwidth=8, width=3, linewidth=0., alpha=0.5)
            #        )
            '''
            fb_k_chi = chisquared_p02k_mock(ell, corr='FB', catalogs=catalogs, n_mocks=n_mocks, 
                Ngrid=Ngrid, rebin='log', use='qpm')
            chi_plot.annotate(r'$\mathtt{B}$', #r'$\mathtt{k}^\mathtt{NN}_{\chi^2}$', 
                    xy=(fb_k_chi, 0.0),
                    xytext=(fb_k_chi*0.9, factor), 
                    arrowprops=dict(facecolor='k', shrink=0.05, headwidth=8, width=3, linewidth=0., alpha=0.5)
                    )
            '''
        if ell == 0: 
            ypos = 6.
        elif ell == 2: 
            ypos = -14.
        chi_plot.text(1.2*10**-2, ypos, 
                r'$\mathtt{k}^{\mathtt{LRec}; \mathtt{l='+str(ell)+'}}_{\chi^2}$', fontsize=24)
    chi_plot.plot(np.arange(10**-2, 10**0, 0.01), np.zeros(len(np.arange(10**-2, 10**0, 0.01))), c='k', lw=2, ls='--') 
    # x-axis 
    sub.set_xlim([10**-2,10**0])
    sub.set_xscale('log')
    sub.set_xticklabels([])

    chi_plot.set_xlim([10**-2,10**0])
    chi_plot.set_xscale('log') 
    chi_plot.set_xlabel('k (h/Mpc)', fontsize=24)
    # y-axis    
    sub.set_ylim([-0.01,0.30])
    chi_plot.set_ylim([-20, 20]) 
    chi_plot.set_yticklabels([])
    chi_plot.set_yticks([])

    cv_labels = [r"$\sigma_\mathtt{l}\mathtt{(k)}/\mathtt{P_l(k)}$"]
    sub.legend([bigmd, nseries, qpm]+[nseries_cv], [bmd_label, nsr_label, qpm_label]+cv_labels, loc='upper left', scatterpoints=1, 
            prop={'size':20}, handletextpad=0.4, borderpad=1.)#, ncol=len(mock_list)-1)

    sub.set_ylabel(
            r"$\mathtt{1- \overline{P_0^{LRec}}/\overline{P_0^{true}}(k)}$", 
            fontsize=24)

    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    if rebin is None: 
        rebin_str = ''
    else:
        rebin_str = '.rebin'+str(rebin)+'x'
    fig_file = ''.join(['figure/fc_paper/', 
        'mock_catalog_dlospeak_true_P0k_norm_resid', rebin_str, '.png']) 
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()

def fig_dlospeakknown_p02k_mock_resid(n_mocks, Ngrid=960): 
    ''' Figure that compares the P_l(k) residual of the dLOS Peak Known corrected 
    power spectrum monopole and quadrupole of specified mock catalogs. 
    
    P_l(k)^dlospeakknown - P_l(k)^true l = 0, 2 for 

    Parameters
    ----------
    n_mocks : int 
        int of mocks to average the P_l(k) over. Same length as 
        the list of catalog names. 

    Ngrid : int 
        int that specify the Ngrids of the power spectrum 
        monopole and quadrupole. Ngrid list has two elements. 
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14,8))
    all_fig = fig.add_subplot(111)

    dlospeak_dict = {
            'nseries': {
                'name': 'dlospeakknown', 
                'fit': 'gauss', 
                'sigma': 3.9, 
                'fpeak': 0.68
                }
            }

    fig_data = [] 

    sub = fig.add_subplot(111)
    ell = 2 

    # dLOS peak correction parameters hardcoded for 
    # each of the mock catalogs
    cat_label = nsr_label
    cat_color = nsr_color
    dlospeak_corrdict = dlospeak_dict['nseries']  

    tru_catdict = {'name': 'nseries', 'n_mock': 1}
    tru_corrdict = {'name': 'true'}
    tru_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell}
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
    tru_avg_spec = AvgSpec(84, 'pk', tru_cat_corr)
    tru_avg_spec.read()

    dlospeak_avg_spec = AvgSpec(n_mocks, 'pk', dlospeak_cat_corr)
    dlospeak_avg_spec.read()
        
    k_arr = tru_avg_spec.k
    tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
    dlospeak_avg_pk = getattr(dlospeak_avg_spec, ''.join(['p', str(ell), 'k']))

    if n_mocks > 1: 
        pk_err = tru_avg_spec.stddev()

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

    mock = sub.scatter(k_arr, dlospeak_avg_pk - tru_avg_pk, lw=0, color=cat_color)

    # Axes 
    sub.set_xlim([10**-2,10**0])
    sub.set_xscale('log')

    if ell == 0: 
        sub.set_ylim([-750.,750.])
        sub.set_xticklabels([])
        plt.sca(sub)
        #plt.yticks(list(np.arange(-1000., 1250., 500.)) )
        plt.yticks(list(np.arange(-750., 1000., 250.)) )
        cv_label = [r"$\sigma_\mathtt{l} \mathtt{(k)}$"]
        sub.legend([samp_var], cv_label, loc='upper right', scatterpoints=1, 
                prop={'size':25}, handletextpad=0.5, borderpad=0.5) 
        sub.text(1.5*10**-2, 500., r"$\mathtt{l = 0}$", fontsize=24)
    elif ell == 2: 
        #sub.set_ylim([-0.1,1.2])
        sub.set_ylim([-1000.,1000.])
        sub.text(1.5*10**-2, 600., r"$\mathtt{l = 2}$", fontsize=24)
        plt.sca(sub)
        plt.yticks(list(np.arange(-1000., 1000., 500.)) )
        plt.xticks([10**-2, 10**-1])
        sub.set_xlabel('k (h/Mpc)', fontsize=24)
        plt.sca(sub)
        plt.xticks([10**-2, 10**-1])

    fig_data.append({
        'ell': ell,
        'name': 'nseries', 
        'k': k_arr,
        'true_Pk': tru_avg_pk, 
        'corr_Pk': dlospeak_avg_pk, 
        'resid': dlospeak_avg_pk - tru_avg_pk,
        'cv': pk_err
        })

    all_fig.set_xticklabels([])
    all_fig.set_yticklabels([])
    all_fig.set_ylabel(
            r"$\mathtt{\overline{P_l^{NN}} - \overline{P_l^{true}}(k)}$", 
            fontsize=24, labelpad=50)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)

    fig_file = 'figure/fc_paper/mock_catalog_dlospeakknown_true_Plk_resid.png' 
    fig_data_file = 'figure/fc_paper/mock_catalog_dlospeakknown_true_Plk_resid.figdata.p' 
    pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()

def fig_tophatconv_delP(n_mocks=[84, 84], Ngrid=[960, 480], figdata=False, k_rebin=None):
    ''' Figure that compares Del P calculated from the tophat convolution to the
    empirical Del P from the Nseries mock catalog powerspectrum     
    
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

        # combined sigma_Pk (THIS IS WRONG SINCE THE ERRORS ARE CORRELATED)
        #comb_pk_err = np.sqrt(true_pk_err**2 + upw_pk_err**2) 
        #samp_var = sub.fill_between( k_arr, delP_emp - comb_pk_err, delP_emp + comb_pk_err, lw=0, color='grey', alpha=0.6)
        samp_var = sub.fill_between( 
                k_arr, 
                delP_emp - delP_stddev, 
                delP_emp + delP_stddev, 
                lw=0, color='grey', alpha=0.6
                )
        emp_plot = sub.plot(k_arr, delP_emp, lw=2, c='k', ls='--')

        # Del P Uncorrelated
        k_uncorr, delP_uncorr = four_corr.DelPuncorr(ell, fs=fs, rc=rc, k_arr=k_arr)

        # Del P Correlated
        k_corr, delP_corr = four_corr.DelPcorr_pkmu(ell, fs=fs, rc=rc, 
                fold=fold, rebin=rebin, dqdmu=True, dmudq=False)
        corr_interp = interp1d(k_corr, delP_corr, kind='cubic')

        k_range = np.where((k_arr > k_corr[0]) & (k_arr < k_corr[-1]))

        delP_tophat = delP_uncorr[k_range] + corr_interp(k_arr[k_range])
        tophat_plot = sub.scatter(k_arr[k_range], delP_tophat, lw=0, s=5, c=cat_color)

        # x-axis 
        sub.set_xscale('log')
        sub.set_xlim([10**-2,10**0])
        sub.set_xlabel('k (h/Mpc)', fontsize=24)
        # y-axis 
        if ell == 0: 
            sub.text(1.5*10**-2, 125., r"$\mathtt{l = 0}$", fontsize=24)
            sub.set_ylim([-1100., 300.])
            sub.set_ylabel(r"$\mathtt{P_l(k)}$ Residuals", fontsize=24)
            
            sub.legend(emp_plot, [r"$\mathtt{P^{NN}_l(k) - P^{true}_l(k)}$"], 
                    loc='lower right', scatterpoints=1, 
                    prop={'size':20}, handletextpad=0.2, borderpad=0.5) 
        elif ell == 2: 
            sub.text(1.5*10**-2, 900., r"$\mathtt{l = 2}$", fontsize=24)
            sub.set_ylim([-1200., 1200.])
            #sub.set_ylabel(r"$\mathtt{P_2(k)}$ Residuals", fontsize=24)#, labelpad=20)
            sub.yaxis.tick_right()
            #sub.yaxis.set_label_position('right')
            samp_leg = sub.legend(
                    [samp_var], 
                    [r'$\sigma_{\Delta \mathtt{P_l}}$'], 
                    loc='upper right', scatterpoints=1, 
                    prop={'size':20}, handletextpad=0.2, borderpad=0.5) 
            # [r'$\sqrt{\sigma_\mathtt{l;\,true}^2 + \sigma_\mathtt{l;\,NN}^2}$'], 
            plt.gca().add_artist(samp_leg)
            sub.legend([tophat_plot], [r"Eq.$26+34$"], #r"Nseries Box $\mathtt{\Delta P(k)}$"], 
                    loc='lower right', scatterpoints=1, 
                    prop={'size':20}, handletextpad=-0.2, borderpad=0.5, 
                    markerscale=7, scatteryoffsets=[0.5]) 

        fig_data.append({
            'ell': ell,
            'name': 'nseries', 
            'k_emp': k_arr,
            'k_tophat': k_arr[k_range],
            'delP_emp': delP_emp, 
            'delP_tophat': delP_tophat,
            'delP_stddev': delP_stddev
            })
        #'comb_sigma_pk': comb_pk_err

    if k_rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(k_rebin)+'x'
    fig.subplots_adjust(wspace=0.15, hspace=0.15)

    fig_data_file = ''.join(['figure/fc_paper/mock_catalog_tophatconv_upw_delPlk', rebin_str, '.figdata.p'])
    pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig_file = ''.join(['figure/fc_paper/mock_catalog_tophatconv_upw_delPlk', rebin_str, '.png'])
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()
    return None

def fig_tophatconv_delP_untrusted(ktrust, n_mocks=[84, 84], Ngrid=[960, 480], order=10, k_rebin=None, figdata=False): 
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
        
        # sigma_Pk of P^true(k) and P^upw(k) 
        #true_pk_err = true_avg_spec.stddev()
        #upw_pk_err = upw_avg_spec.stddev()
        ## combined sigma_Pk 
        #comb_pk_err = np.sqrt(true_pk_err**2 + upw_pk_err**2) 
        
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
        tophat_plot = sub.scatter(k_arr[lp_lim], delpcorr_untrust[lp_lim], 
                lw=0, s=8, c=cat_color)
    
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
            sub.text(1.5*10**-2, 125., r"$\mathtt{l = 0}$", fontsize=24)
            #sub.set_ylim([-600., 0.])
            sub.set_ylim([-1100., 300.])
            sub.set_ylabel(
                    r"$\Delta \mathtt{P^{corr}_l(k)} |_\mathtt{q=k_{trust}}^{\mathtt{q=}\infty}$", 
                    fontsize=24)
            sub.legend(
                    [emp_plot, poly2_plot, tophat_plot], 
                    [r"Eq.38 [Mocks]", r"Eq.37 for $\mathtt{l' \leq 2}$", 
                        r"Eq.37 for $\mathtt{l' \leq "+str(order)+"}$"], 
                    loc='lower left', scatterpoints=1, 
                    prop={'size':20}, handletextpad=0.0, borderpad=0.5) 
            
        elif ell == 2: 
            sub.text(1.5*10**-2, 900., r"$\mathtt{l = 2}$", fontsize=24)
            #sub.set_ylim([-500., 500.])
            sub.set_ylim([-1200., 1200.])
            #sub.set_ylabel(r"$\mathtt{P_2(k)}$ Residuals", fontsize=24)#, labelpad=20)
            sub.yaxis.tick_right()
            #sub.yaxis.set_label_position('right')
            samp_leg = sub.legend(
                    [samp_var], 
                    [r'$\sigma_{\Delta \mathtt{P_l}}$'], 
                    loc='upper right', scatterpoints=1, 
                    prop={'size':20}, handletextpad=0.2, borderpad=0.5) 
        # [r'$\sqrt{\sigma_\mathtt{true}^2 + \sigma_\mathtt{NN}^2}$'], 
        #fig_data.append({
        #    'ell': ell,
        #    'name': 'nseries', 
        #    'k_emp': k_arr,
        #    'k_tophat': k_arr[k_range],
        #    'delP_emp': delP_emp, 
        #    'delP_tophat': delP_tophat,
        #    'comb_sigma_pk': comb_pk_err
        #    })
    if k_rebin is None: 
        rebin_str = ''
    else: 
        rebin_str = '.rebin'+str(k_rebin)+'x'

    fig.subplots_adjust(wspace=0.15, hspace=0.15)

    fig_file = 'figure/fc_paper/mock_catalog_tophatconv_delPlk_untrusted_'+str(order)+rebin_str+'.png' 
    #fig_data_file = 'figure/fc_paper/mock_catalog_tophatconv_upw_delPlk.figdata.p' 
    #pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close()
    return None

def fig_mock_dlos(catalogs=['nseries'], n_mocks=[84], peak_range=[-15.0, 15.0]): 
    ''' dLOS distribution plot for fc_paper
    '''
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=[14, 5])
    sub = fig.add_subplot(111)

    catalogs = ['cmass'] + catalogs
    n_mocks = [1] + n_mocks
    
    fig_data = [] 
    for i_cat, cat in enumerate(catalogs): 
        cat_dlos = [] 
        for i_mock in xrange(1, n_mocks[i_cat]+1): 
            cat_corr = {
                    'catalog': {'name': cat, 'n_mock': i_mock}, 
                    'correction': {'name': 'upweight'}
                    }
            deelos = Dlos(cat_corr)
            deelos.read()
            cat_dlos += list(deelos.dlos)

        cat_deelos = Dlos(cat_corr)
        cat_deelos.dlos = np.array(cat_dlos)

        # Freedman Diaconis binsize
        fd_binsize = cat_deelos.fd_binsize(peak_range = peak_range)
        
        if cat == 'cmass': 
            dlos_mid, dlos_dist = cat_deelos.dlos_dist(binsize = 0.5, normed=True)
        else:
            dlos_mid, dlos_dist = cat_deelos.dlos_dist(binsize = fd_binsize, normed=True)

        fig_data.append(
                {
                    'name': cat,
                    'n_mocks': n_mocks[i_cat], 
                    'dlos': np.array(cat_dlos), 
                    'dlos_mid': dlos_mid, 
                    'dlos_hist': dlos_dist
                })

        if cat == 'qpm': 
            cat_label = qpm_label 
            cat_color = qpm_color
        elif cat == 'bigmd': 
            cat_label = bmd_label
            cat_color = bmd_color
        elif cat == 'nseries': 
            cat_label = nsr_label
            cat_color = nsr_color
        elif cat == 'cmass': 
            cat_label = cms_label
            cat_color = cms_color
        else: 
            raise NotImplementedError 

        sub.plot(dlos_mid, dlos_dist, lw=4, color=cat_color, label=cat_label)
    
    sub.set_xlim([-45.0, 45.0])
    sub.set_xlabel('$\mathtt{d_{LOS}}$ (Mpc)', fontsize=30)
    sub.set_ylabel('$\mathtt{p(d_{LOS})}$', fontsize=30)
    sub.legend(loc='upper left', borderpad=1.0)
    
    fig_file = 'figure/fc_paper/mock_catalog_dlos.png'
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    figdata_file = 'figure/fc_paper/mock_catalog_dlos.figdata.p'
    pickle.dump(fig_data, open(figdata_file, 'wb'))
    plt.close()

def fig_2pcf_tophat(n_mocks=5, scale='small', fs=0.6, rc=0.43): 
    '''
    Plot 1 - (1+xi^fc)/(1+xi^true) from CUTE 2PCF code
    '''
    prettyplot()
    pretty_colors = prettycolors()

    if scale == 'large': 
        contour_range = np.arange(-0.05, 0.05, 0.005)
        n_rp = 40 
        n_pi = 40
    elif scale == 'small': 
        contour_range = np.arange(-0.1, 0.11, 0.01)
    elif scale == 'smaller': 
        contour_range = np.arange(-0.5, 0.5, 0.05)
    elif scale == 'verysmall': 
        contour_range = 20
        n_rp = 20 
        n_pi = 20
    elif scale == '5x5':    # 5 Mpc x 5 Mpc
        #contour_range = np.arange(-0.1, 0.11, 0.01)
        contour_range = 20
        n_rp = 100 
        n_pi = 100
    
    fig = plt.figure(figsize=(24,10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.25, 1]) 
    sub = fig.add_subplot(gs[0])#121)

    corrections = ['true', 'upweighted']
    for corr in corrections:    # for each correction
        for i_mock in n_mocks:  # for each mocks
            corr_file = ''.join([
                '/mount/riachuelo1/hahn/2pcf/corr/',
                'corr_2pcf_CutskyN', 
                str(i_mock), '.', 
                corr, 
                '.cute2pcf.', 
                scale, 
                'scale.dat'
                ])
            print corr_file 

            tpcf = np.loadtxt(corr_file, unpack=True)
    
            if i_mock == n_mocks[0]: 
                rp_bins = tpcf[0].reshape(n_rp, n_pi)[0]
                pi_bins = tpcf[1].reshape(n_rp, n_pi)[:,0]
                r_p, pi = np.meshgrid(rp_bins, pi_bins)

                twop_corr = tpcf[2].reshape(n_rp, n_pi)

            else: 
                twop_corr += tpcf[2].reshape(n_rp, n_pi)

        twop_corr /= np.float(len(n_mocks))  # average 2PCF

        if corr == 'true': 
            true_corr = twop_corr.T
            continue

        # contour of 1 - (1 + xi^fc)/(1+xi^true)
        residual_2pcf = 1.0 - (1.0 + twop_corr.T)/(1.0 + true_corr)
        #print np.max(residual_2pcf)
        cont = sub.contourf(
                r_p, pi, 
                residual_2pcf, 
                contour_range, 
                cmap=plt.cm.afmhot
                )
        plt.colorbar(cont)

        sub.vlines(rc, 0.0, np.max(r_p), lw=8, linestyle='--', color='k')

        sub.set_ylabel('$\pi \; \mathtt{(Mpc/h)}$', fontsize=40)
        sub.set_xlabel('$\mathtt{r_{p} \; (Mpc/h)}$', fontsize=40)
        sub.set_xlim([np.min(rp_bins), np.max(rp_bins)])
        sub.set_ylim([np.min(pi_bins), np.max(pi_bins)])
        
        #sub.set_title(r"$1 - (1 + \xi^\mathtt{"+corr.upper()+r"})/(1+ \xi^\mathtt{TRUE})$", fontsize=40)
    
    sub_proj = fig.add_subplot(gs[1])#122)

    proj_resid = np.sum(residual_2pcf, axis=0)
    sub_proj.plot(r_p[0], proj_resid/np.float(n_pi), 
            lw=4, c=pretty_colors[7], label=r'Nseries') 
    # top hat function for comparison  
    rbin = np.linspace(0.0, 50.0, 10000)
    inhat = np.where(rbin < rc)
    tophat = np.repeat(0., len(rbin))
    tophat[inhat] = fs
    sub_proj.plot(rbin, tophat, lw=3, ls='--', c='k', 
            label=r'$\mathtt{f_s W_{fc}(r_p)}$')
    
    # x-axis
    sub_proj.set_xlim([0, 0.5*np.max(r_p[0])])
    sub_proj.set_xlabel('$\mathtt{r_{p} \; (Mpc/h)}$', fontsize=40)
    # y-axis
    sub_proj.set_ylim([-0.1, 1.0])
    sub_proj.set_ylabel(
            r"$1 - (1 + \xi^\mathtt{NN})/(1+ \xi^\mathtt{true})$", 
            fontsize=30)

    sub_proj.legend(loc='upper right', fontsize=40)

    fig_name = ''.join([
            'figure/fc_paper/', 
            '2pcf_Nseries_', corr, '_', str(len(n_mocks)), 'mocks.', scale, '.tophat.png'
            ])
    fig.savefig(fig_name, bbox_inches="tight", dpi=150)
    plt.close()

def table_dlospeak_fit(catalogs, fit='gauss'): 
    ''' Construct LaTex table from best fit fpeak and sigma from peak of the 
    dlos distribution from specified catalogs.  

    Parameters
    ----------
    catalogs : list 
        List of strings that specify the names of the catalogs. 

    fit : string
        Stirng that specifies the functional form of the dlos Peak fit. 
        Default is 'gauss'. 
    '''
    header_str = ''.join([r'\begin{tabular}{ccc} \hline \hline', '\n', 
        r'Catalog &$\sigma$ ($\mathrm{Mpc}$) & $f_{\mathrm{peak}}$\\ \hline', '\n'])

    row_str = ''
    bestfit_file = lambda c, f: ''.join(['dlos/dlos_fit', '.', c, '.', f, '_fit',  '.p'])
    for cat in catalogs: 
        save_dict = pickle.load(open(bestfit_file(cat, fit), 'rb'))
        fpeak = save_dict['fpeak']
        sigma = save_dict['sigma']

        if cat == 'nseries': 
            cat_name = 'Nseries'
        elif cat == 'qpm': 
            cat_name = 'QPM'
        elif cat == 'bigmd': 
            cat_name = 'BigMultiDark'
        elif cat == 'cmass': 
            cat_name = 'CMASS'

        row_str += '&'.join([cat_name, str(round(sigma, 2)), str(round(fpeak, 2))])

        if cat == catalogs[-1]: 
            row_str += r'\\ \hline'+'\n'
        else: 
            row_str += r'\\'+'\n'

    tail_str = '\end{tabular} \par'

    table_str = ''.join([header_str, row_str, tail_str])
        
    return table_str 

# Chi^2 
def fig_mock_covariance(ell, catalogs=['nseries'], n_mocks=[84], Ngrid=960, rebin=10): 
    ''' Plot the covariance matrix for mock catalog power spectrum .
    '''
    if len(catalogs) != len(n_mocks): 
        raise ValueError
     
    #fig_data = [] 
    mock_list, mock_label_list = [], [] 
    for i_cat, cat in enumerate(catalogs):
        catdict = {'name': cat, 'n_mock': 1}
        specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell}
        tru_cat_corr = {
                'catalog': {'name': catdict['name'], 'n_mock': 1}, 
                'correction': {'name': 'true'}, 
                'spec': specdict
                }
        upw_cat_corr = {
                'catalog': {'name': catdict['name'], 'n_mock': 1}, 
                'correction': {'name': 'upweight'}, 
                'spec': specdict
                }

        # figure to check that the rebinning of the power spectrum is  
        # being performed properly 
        test_fig = plt.figure(figsize=[10,10])
        test_sub = test_fig.add_subplot(111)
        test_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr) 
        test_avg_spec.read(rebin=None)
        test_sub.scatter(
                test_avg_spec.k, 
                getattr(test_avg_spec, ''.join(['p', str(ell), 'k'])), 
                s=10, lw=0)
        # test different rebinning sizes
        for rb in [5, 10, 20]: 
            test_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr) 
            test_avg_spec.read(rebin=rb)
            test_sub.plot(
                    test_avg_spec.k, 
                    getattr(test_avg_spec, ''.join(['p', str(ell), 'k'])), 
                    label='rebin='+str(rb), lw=2)
        test_sub.legend(loc='lower left')
        test_sub.set_xlim([10**-3, 10**0])
        test_sub.set_yscale('log')
        test_sub.set_xscale('log')
        test_file = ''.join([
            'figure/fc_paper/', 
            'test_p', str(ell), 'k.', 
            cat, '.k_rebinning.png'
            ])
        test_fig.savefig(test_file, bbox_inches='tight', dpi=150)
        
        # re-binned true + upweight power spectrum
        tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
        tru_avg_spec.read(rebin=rebin)
        upw_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', upw_cat_corr)
        upw_avg_spec.read(rebin=rebin)
    
        k_arr = tru_avg_spec.k
        tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
        upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
        counts = tru_avg_spec.count
         
        if n_mocks[i_cat] <= 1: 
            continue 

        full_cov = tru_avg_spec.covariance(rebin=rebin)

        inear = (np.abs(k_arr - 0.5)).argmin() + 2
        prettyplot()
        fig = plt.figure()
        sub = fig.add_subplot(111)
        colormap = plt.cm.jet
        norm = mpl.colors.SymLogNorm(1000., vmin=-1.*10**7, vmax=10**9.5, clip=True)
        cont = sub.pcolormesh(
                k_arr[:inear], k_arr[:inear], full_cov[:inear, :inear], 
                norm=norm, cmap=colormap)
        
        sub.set_xlim([k_arr[0], 0.5])
        sub.set_ylim([k_arr[0], 0.5])
        sub.set_yscale('log')
        sub.set_ylabel('k', fontsize=30)
        sub.set_xscale('log')
        sub.set_xlabel('k', fontsize=30)
        plt.colorbar(cont)
        
        rebin_str = ''
        if isinstance(rebin, int): 
            rebin_str = '.rebin'+str(rebin)
        elif isinstance(rebin, str): 
            rebin_str = '.rebin'+rebin

        fig_file = ''.join([
            'figure/fc_paper/', 
            'p', str(ell), 'k.', 
            cat, '.covariance', 
            rebin_str, '.png'
            ])
        fig.savefig(fig_file, bbox_inches='tight', dpi=150)

    return None 

def chisquared_p02k_mock(ell, chi_goal=1., corr='NN', catalogs=['nseries'], n_mocks=[84], Ngrid=960, rebin=10, use='nseries'): 
    ''' Calculate the k value where the reduced chi^2 value is 1.  
    '''
    if len(catalogs) != len(n_mocks): 
        raise ValueError
    
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
    
    k_chi, k_w_chi = [], [] 
    chi_w, chi = [], [] 
    for i_cat, cat in enumerate(catalogs):
        tru_catdict = {'name': cat, 'n_mock': 1}
        tru_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell}
        if corr == 'NN': 
            corr_corrdict = {'name': 'upweight'} 
        elif corr == 'LR': 
            corr_corrdict = dlospeak_dict[cat]
        elif corr == 'LK': 
            corr_corrdict = dlospeak_dict[cat]
            corr_corrdict['name'] = 'dlospeakknown'
        elif corr == 'FB': 
            corr_corrdict = {'name': 'floriansn'} 
        elif corr == 'GM': 
            corr_corrdict = {'name': 'hectorsn'} 

        tru_cat_corr = {
                'catalog': {'name': tru_catdict['name'], 'n_mock': 1}, 
                'correction': {'name': 'true'}, 
                'spec': tru_specdict
                }
        upw_cat_corr = {
                'catalog': {'name': tru_catdict['name'], 'n_mock': 1}, 
                'correction': corr_corrdict, 
                'spec': tru_specdict
                }

        n_mock_cat = n_mocks[i_cat]

        tru_avg_spec = AvgSpec(n_mock_cat, 'pk', tru_cat_corr)
        tru_avg_spec.read(rebin=rebin)

        upw_avg_spec = AvgSpec(n_mock_cat, 'pk', upw_cat_corr)
        upw_avg_spec.read(rebin=rebin)
    
        k_arr = tru_avg_spec.k
        tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
        upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
        counts = tru_avg_spec.count
        
        # covariance matrix
        if n_mocks[i_cat] > 1: 
            full_cov = tru_avg_spec.covariance(rebin=rebin)
        else: 
            # for bigMD use either Nseries or QPM, whichever is specified by the 
            # 'use' kwarg
            if use == 'nseries':
                borrow_cat_corr = {
                        'catalog': {'name': 'nseries', 'n_mock': 1}, 
                        'correction': {'name': 'true'}, 
                        'spec': tru_specdict
                        }
                borrow_n_mock = 84
            elif use == 'qpm':
                borrow_cat_corr = {
                    'catalog': {'name': 'qpm', 'n_mock': 1}, 
                    'correction': {'name': 'true'}, 
                    'spec': tru_specdict
                    }
                borrow_n_mock = 100 

            borrow_avg_spec = AvgSpec(borrow_n_mock, 'pk', borrow_cat_corr)
            borrow_avg_spec.read(rebin=rebin)
            full_cov = borrow_avg_spec.covariance(rebin=rebin)
            n_mock_cat = borrow_n_mock 
    
        print k_arr
        print upw_avg_pk - tru_avg_pk 
    
        # calculate chi^2(k) 
        i_start = (np.abs(k_arr - 0.05)).argmin()
        once = None
        for ik in range(i_start, len(k_arr)): 
            cov = full_cov[:ik, :ik]
            #cov = np.diag(np.diag(cov))
            P_resid = upw_avg_pk[:ik] - tru_avg_pk[:ik]
            if corr == 'FB':
                P_resid -= 250.
            #print cov.shape
            #print np.linalg.cond(cov)

            if np.linalg.cond(cov) < 10**20:    # otherwise the cov matrix can't be inverted
                inv_chi2 = np.dot(P_resid, np.dot(np.linalg.inv(cov), P_resid))     #/np.float(len(P_resid))
                pinv_chi2 = np.dot(P_resid, np.dot(np.linalg.pinv(cov), P_resid))   #/np.float(len(P_resid))
                if pinv_chi2 > chi_goal: 
                    if once is None: 
                        k_chi.append(0.5*(k_arr[ik-2] + k_arr[ik-1]))
                        k_w_chi.append(
                                (k_arr[ik-2] * counts[ik-2] + k_arr[ik-1] * counts[ik-1])/(counts[ik-2] + counts[ik-1]))
                        prev_inv_chi2 = np.dot(P_resid[:-1], np.dot(np.linalg.inv(cov[:-1, :-1]), P_resid[:-1])) #/np.float(len(P_resid)-1)
                        chi.append(inv_chi2)
                        chi_w.append(0.5 * (inv_chi2 + prev_inv_chi2))

                        #print cat, cov.shape
                        #print 'k_avg = ',  0.5*(k_arr[ik-2] + k_arr[ik-1])
                        #print 'Chi^2 =', inv_chi2, ' ; Chi^2 =', pinv_chi2, ' ', np.linalg.cond(cov)
                        once = 'twice'
                    else: 
                        continue
            else: 
                continue

    #print 'k_chi = ', np.mean(k_chi)
    #print 'weighted k_chi = ', np.mean(k_w_chi)
    #print 'chi^2 = ', np.mean(chi)
    #print 'weighted chi^2 = ', np.mean(chi_w)
    #print 'k_chi', np.mean(k_chi), '; k_w_chi ', k_w_chi
    return np.mean(k_w_chi) 

def chisquared_p02k_mock_SNoffset(ell, chi_goal=1., corr='NN', catalogs=['nseries'], n_mocks=[84], Ngrid=960, 
        rebin=10, use='nseries'): 
    ''' Calculate the k value where the reduced chi^2 value is 1.  
    '''
    if len(catalogs) != len(n_mocks): 
        raise ValueError

    offsets = np.arange(-550., -150., 50.)
    k_offsets = [] 
    for offset in offsets: 
        k_chi, k_w_chi = [], [] 
        chi_w, chi = [], [] 
        for i_cat, cat in enumerate(catalogs):
            tru_catdict = {'name': cat, 'n_mock': 1}
            tru_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell}
            corr_corrdict = {'name': 'floriansn'} 

            tru_cat_corr = {
                    'catalog': {'name': tru_catdict['name'], 'n_mock': 1}, 
                    'correction': {'name': 'true'}, 
                    'spec': tru_specdict
                    }
            upw_cat_corr = {
                    'catalog': {'name': tru_catdict['name'], 'n_mock': 1}, 
                    'correction': corr_corrdict, 
                    'spec': tru_specdict
                    }

            n_mock_cat = n_mocks[i_cat]

            tru_avg_spec = AvgSpec(n_mock_cat, 'pk', tru_cat_corr)
            tru_avg_spec.read(rebin=rebin)

            upw_avg_spec = AvgSpec(n_mock_cat, 'pk', upw_cat_corr)
            upw_avg_spec.read(rebin=rebin)
        
            k_arr = tru_avg_spec.k
            tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
            upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
            counts = tru_avg_spec.count
            
            # covariance matrix
            if n_mocks[i_cat] > 1: 
                full_cov = tru_avg_spec.covariance(rebin=rebin)
            else: 
                # for bigMD use either Nseries or QPM, whichever is specified by the 
                # 'use' kwarg
                if use == 'nseries':
                    borrow_cat_corr = {
                            'catalog': {'name': 'nseries', 'n_mock': 1}, 
                            'correction': {'name': 'true'}, 
                            'spec': tru_specdict
                            }
                    borrow_n_mock = 84
                elif use == 'qpm':
                    borrow_cat_corr = {
                        'catalog': {'name': 'qpm', 'n_mock': 1}, 
                        'correction': {'name': 'true'}, 
                        'spec': tru_specdict
                        }
                    borrow_n_mock = 100 

                borrow_avg_spec = AvgSpec(borrow_n_mock, 'pk', borrow_cat_corr)
                borrow_avg_spec.read(rebin=rebin)
                full_cov = borrow_avg_spec.covariance(rebin=rebin)
                n_mock_cat = borrow_n_mock 
            
            # calculate chi^2(k) 
            i_start = (np.abs(k_arr - 0.05)).argmin()
            once = None
            for ik in range(i_start, len(k_arr)): 
                cov = full_cov[:ik, :ik]
                #cov = np.diag(np.diag(cov))
                P_resid = upw_avg_pk[:ik] - tru_avg_pk[:ik] + offset
                #print cov.shape
                #print np.linalg.cond(cov)

                if np.linalg.cond(cov) < 10**20:    # otherwise the cov matrix can't be inverted
                    inv_chi2 = np.dot(P_resid, np.dot(np.linalg.inv(cov), P_resid))     #/np.float(len(P_resid))
                    pinv_chi2 = np.dot(P_resid, np.dot(np.linalg.pinv(cov), P_resid))   #/np.float(len(P_resid))
                    if pinv_chi2 > chi_goal: 
                        if once is None: 
                            k_chi.append(0.5*(k_arr[ik-2] + k_arr[ik-1]))
                            k_w_chi.append(
                                    (k_arr[ik-2] * counts[ik-2] + k_arr[ik-1] * counts[ik-1])/(counts[ik-2] + counts[ik-1]))
                            prev_inv_chi2 = np.dot(P_resid[:-1], np.dot(np.linalg.inv(cov[:-1, :-1]), P_resid[:-1])) #/np.float(len(P_resid)-1)
                            chi.append(inv_chi2)
                            chi_w.append(0.5 * (inv_chi2 + prev_inv_chi2))

                            #print cat, cov.shape
                            #print 'k_avg = ',  0.5*(k_arr[ik-2] + k_arr[ik-1])
                            #print 'Chi^2 =', inv_chi2, ' ; Chi^2 =', pinv_chi2, ' ', np.linalg.cond(cov)
                            once = 'twice'
                        else: 
                            continue
                else: 
                    continue
        k_offsets.append(np.mean(k_w_chi))
    
    print k_offsets
    print max(k_offsets) 
    print k_offsets.index(max(k_offsets))
    print offsets[k_offsets.index(max(k_offsets))]
    #print 'k_chi = ', np.mean(k_chi)
    #print 'weighted k_chi = ', np.mean(k_w_chi)
    #print 'chi^2 = ', np.mean(chi)
    #print 'weighted chi^2 = ', np.mean(chi_w)
    #print 'k_chi', np.mean(k_chi), '; k_w_chi ', k_w_chi
    return np.mean(k_w_chi) 

def fig_chisquared_NN_p02k_mock(ell, catalogs=['nseries'], n_mocks=[84], Ngrid=960, rebin=10, use='nseries', reduced=False): 
    ''' Calculate the k value where chi^2 value is 1.  
    '''
    if len(catalogs) != len(n_mocks): 
        raise ValueError
    
    prettyplot()
    fig = plt.figure()
    sub = fig.add_subplot(111)
     
    fig_data = [] 
    mock_list, mock_label_list = [], [] 
    for i_cat, cat in enumerate(catalogs):
        if cat == 'qpm': 
            cat_label = qpm_label 
            cat_color = qpm_color
        elif cat == 'bigmd': 
            cat_label = bmd_label
            cat_color = bmd_color
        elif cat == 'nseries': 
            cat_label = nsr_label 
            cat_color = nsr_color 
        else: 
            raise NotImplementedError 

        tru_catdict = {'name': cat, 'n_mock': 1}
        tru_corrdict = {'name': 'true'}
        upw_corrdict = {'name': 'upweight'}
        tru_specdict = {
                'P0': 20000,
                'Lbox': 3600, 
                'Ngrid': Ngrid, 
                'ell': ell 
                }
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
        tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
        tru_avg_spec.read(rebin=rebin)
        
        upw_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', upw_cat_corr)
        upw_avg_spec.read(rebin=rebin)
    
        k_arr = tru_avg_spec.k
        tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
        upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
        counts = tru_avg_spec.count
          
        if n_mocks[i_cat] > 1: 
            full_cov = tru_avg_spec.covariance(rebin=rebin)
        else: 
            if use == 'nseries':
                nseries_cat_corr = {
                        'catalog': {'name': 'nseries', 'n_mock': 1}, 
                        'correction': tru_corrdict, 
                        'spec': tru_specdict
                        }
                nseries_avg_spec = AvgSpec(84, 'pk', nseries_cat_corr)
                nseries_avg_spec.read(rebin=rebin)
                full_cov = nseries_avg_spec.covariance(rebin=rebin)
            elif use == 'qpm':
                qpm_cat_corr = {
                    'catalog': {'name': 'qpm', 'n_mock': 1}, 
                    'correction': tru_corrdict, 
                    'spec': tru_specdict
                    }
                qpm_avg_spec = AvgSpec(100, 'pk', qpm_cat_corr)
                qpm_avg_spec.read(rebin=rebin)
                full_cov = qpm_avg_spec.covariance(rebin=rebin)
        
        chi2 = [] 
        for ik in range(1, len(k_arr)): 
            cov = full_cov[:ik, :ik]
            #cov = np.diag(np.diag(cov))
            P_resid = upw_avg_pk[:ik] - tru_avg_pk[:ik]

            if np.linalg.cond(cov) < 10**20: 
                inv_chi2 = np.dot(P_resid, np.dot(np.linalg.inv(cov), P_resid))/np.float(len(P_resid))
                pinv_chi2 = np.dot(P_resid, np.dot(np.linalg.pinv(cov), P_resid))/np.float(len(P_resid))
                
                if round(inv_chi2,2) == round(pinv_chi2,2): 
                    chi2.append(inv_chi2)
                else: 
                    print inv_chi2, pinv_chi2
                    raise ValueError
            else: 
                continue
            
        sub.plot(k_arr[1:], chi2, color=cat_color, label=cat_label)

        #fig_data.append({
        #    'ell': ell,
        #    'name': cat , 
        #    'k': k_arr, 
        #    'true_Pk': tru_avg_pk, 
        #    'upw_Pk': upw_avg_pk,
        #    'resid': upw_avg_pk - tru_avg_pk,
        #    'cv': pk_err                
        #    })

    sub.plot(np.linspace(0.001, 1.0, 1000), np.repeat(1.0, len(np.linspace(0.001, 1.0, 1000))), c='k', lw=4, ls='--')
    sub.set_ylabel(r"Reduced $\chi^\mathtt{2}\mathtt{(k'<k)}$", fontsize=25)
    sub.set_yscale('log')
    sub.set_xscale('log') 
    sub.set_xlim([10**-3, 10**0])
    sub.set_xlabel(r"$\mathtt{k}\; h/\mathtt{Mpc}$", fontsize=25)
    sub.legend(loc='upper left', scatterpoints=1, prop={'size':20}, borderpad=0.1, handletextpad=0.5) 
    
    #pickle.dump(fig_data, open(fig_data_file, 'wb'))
    fig_file = ''.join([
        'figure/fc_paper/', 
        'reduced_chi2.p', str(ell), 'k.NNmethod.png'
        ])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    return None




"""
    def fig_dlospeak_p0k_p2k_mock_catalogs(catalogs=['nseries'], n_mocks=[84], Ngrid=[960]): 
        ''' Figure that compares the ratio of the line-of-sight reconstructed power spectrum 
        monopole and quadrupole over the true power spectrum monopole and quadrupole of 
        specified mock catalogs. 
        
        P_l(k)^dlospeak / P_l(k)^true l = 0, 2 for 

        Parameters
        ----------
        catalogs : list
            List of strings that specify the mock catalogs to include. 
            (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])

        n_mocks : list 
            List of ints of mocks to average the P_l(k) over. Same length as 
            the list of catalog names. 

        Ngrid : list 
            List of ints that specify the Ngrids of the power spectrum 
            monopole and quadrupole. Ngrid list has two elements. 
        '''
        prettyplot()
        pretty_colors = prettycolors()
        
        fig = plt.figure(1, figsize=(14,8))

        if len(catalogs) < 3: 
            for i in xrange(3 - len(catalogs)): 
                catalogs.append(catalogs[0])
                n_mocks.append(n_mocks[0])
        if len(catalogs) != len(n_mocks): 
            raise ValueError
        
        fig_data = [] 
        for i_ell, ell in enumerate([0, 2]):  # monopole or quadrupole
            mock_list, mock_label_list = [], [] 
            for i_cat, cat in enumerate(catalogs):
                sub = fig.add_subplot(2, 3, 3*i_ell + i_cat+1)

                tru_catdict = {'name': cat, 'n_mock': 1}
                tru_corrdict = {'name': 'true'}
                # dLOS peak correction parameters hardcoded for 
                # each of the mock catalogs
                if cat == 'qpm': 
                    cat_label = qpm_label 
                    cat_color = qpm_color
                    dlospeak_corrdict = {
                            'name': 'dlospeak', 
                            'fit': 'gauss', 
                            'sigma': 4.4, 
                            'fpeak': 0.62
                            }
                elif cat == 'bigmd': 
                    cat_label = bmd_label 
                    cat_color = bmd_color
                    dlospeak_corrdict = {
                            'name': 'dlospeak', 
                            'fit': 'gauss', 
                            'sigma': 5.5, 
                            'fpeak': 0.6
                            }
                elif cat == 'nseries': 
                    cat_label = nsr_label
                    cat_color = nsr_color
                    dlospeak_corrdict = {
                            'name': 'dlospeak', 
                            'fit': 'gauss', 
                            'sigma': 3.9, 
                            'fpeak': 0.68
                            }
                else: 
                    raise NotImplementedError 

                tru_specdict = {
                        'P0': 20000,
                        'Lbox': 3600, 
                        'Ngrid': Ngrid[i_ell], 
                        'ell': ell 
                        }
                tru_cat_corr = {
                        'catalog': {'name': tru_catdict['name'], 'n_mock': 1}, 
                        'correction': tru_corrdict, 
                        'spec': tru_specdict
                        }
                dlospeak_cat_corr = {
                        'catalog': {'name': tru_catdict['name'], 'n_mock': 1}, 
                        'correction': dlospeak_corrdict, 
                        'spec': tru_specdict
                        }
                tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
                tru_avg_spec.read()
                
                dlospeak_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', dlospeak_cat_corr)
                dlospeak_avg_spec.read()
            
                k_arr = tru_avg_spec.k
                tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
                dlospeak_avg_pk = getattr(dlospeak_avg_spec, ''.join(['p', str(ell), 'k']))

                if n_mocks[i_cat] > 1: 
                    pk_err = tru_avg_spec.stddev()

                    if ell == 0: 
                        mono_cv = sub.plot(k_arr, np.abs(pk_err/tru_avg_pk), 
                                lw=1, ls='--', color='gray')
                    else: 
                        quad_cv = sub.plot(k_arr, np.abs(pk_err/tru_avg_pk), 
                                lw=1, ls='--', color='gray')
                if ell == 0: 
                    plt_label = None
                elif ell == 2: 
                    plt_label = cat_label

                mock = sub.scatter(k_arr, np.abs(1. - dlospeak_avg_pk/tru_avg_pk), lw=0, color=cat_color)
                if ell == 2: 
                    mock_list.append(mock)
                    mock_label_list.append(plt_label)
                sub.plot((10**-3, 10**0), (0.0,0.0), 'k--')

                # Axes 
                sub.set_xlim([10**-3,10**0])
                sub.set_xscale('log')

                if ell == 0: 
                    sub.set_xticklabels([])
                    sub.set_ylim([-0.05,0.25])
                    if i_cat == len(catalogs)-1: 
                        mono_label = [r"$|\frac{\sigma_\mathtt{l}\mathtt{(k)}}{\mathtt{P_l}(k)}|$"]
                        sub.legend(mono_cv, mono_label, loc='upper left', scatterpoints=1, 
                                prop={'size':25}, handletextpad=0.2, borderpad=0.5) 
                elif ell == 2: 
                    sub.set_ylim([-0.1,1.0])
                    if i_cat == len(catalogs)-1: 
                        sub.legend(mock_list, mock_label_list, loc='upper left', scatterpoints=1, 
                                prop={'size':20}, borderpad=0.1, handletextpad=0.2) 
                    elif i_cat == 0: 
                        plt.sca(sub)
                        plt.yticks(list(np.arange(0.0, 1.0, 0.2)) )
                        plt.xticks([10**-3, 10**-2, 10**-1])
                    else: 
                        sub.set_xlabel('k (h/Mpc)', fontsize=24)
                        plt.sca(sub)
                        plt.xticks([10**-3, 10**-2, 10**-1])
                    #sub.yaxis.tick_right()
                    #sub.yaxis.set_label_position('right')
                if i_cat == 0:
                    sub.set_ylabel(
                            r"$\mathtt{|1-\overline{P_{"+str(ell)+"}^{LRec}}/\overline{P_{"+str(ell)+"}^{true}}|}$", 
                            fontsize=24)
                else: 
                    sub.set_yticklabels([])

                fig_data.append({
                    'ell': ell,
                    'name': cat , 
                    'k': k_arr,
                    'true_Pk': tru_avg_pk, 
                    'upw_Pk': dlospeak_avg_pk, 
                    'norm_resid': np.abs(1.- dlospeak_avg_pk/tru_avg_pk),
                    'norm_cv': np.abs(pk_err/tru_avg_pk)
                    })

        fig.subplots_adjust(wspace=0.0, hspace=0.0)

        fig_file = 'figure/fc_paper/mock_catalog_dlospeak_true_Plk.png' 
        fig_data_file = 'figure/fc_paper/mock_catalog_dlospeak_true_Plk.figdata.p' 
        pickle.dump(fig_data, open(fig_data_file, 'wb'))
        fig.savefig(fig_file, bbox_inches="tight", dpi=150)
        plt.close()

    def fig_tophatconv_p0k_p2k_mock_catalogs(n_mocks=[84, 84], Ngrid=[960, 480]): 
        ''' Figure that compares the ratio of the tophat convolved power spectrum 
        monopole and quadrupole over the upweighted power spectrum monopole and quadrupole of 
        specified mock catalogs. 
        
        P_l(k)^tophatconv / P_l(k)^upweight l = 0, 2 for 

        Parameters
        ----------
        n_mocks : list 
            Length 2 list of ints that specify the number of mocks to average 
            the P_l(k) over. 

        Ngrid : list 
            Length 2 list of ints that specify the Ngrids of the power spectrum 
            monopole and quadrupole, respectively
        '''
        prettyplot()
        pretty_colors = prettycolors()
        
        fig = plt.figure(1, figsize=(14,8))

        fig_data = [] 
        for i_ell, ell in enumerate([0, 2]):  # monopole or quadrupole
            sub = fig.add_subplot(1, 2, i_ell+1)

            # TophatConv correction parameters hardcoded for 
            # each of the mock catalogs
            cat_label = nsr_label
            cat_color = nsr_color
            tophat_conv_corrdict = {
                    'name': 'tophat_conv', 
                    'fs': 0.6, 
                    'rc': 0.43, 
                    'fold': 10, 
                    'rebin': 20
                    }

            upw_catdict = {'name': 'nseries', 'n_mock': 1}
            upw_specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid[i_ell], 'ell': ell}
            upw_cat_corr = {
                    'catalog': upw_catdict, 
                    'correction': {'name': 'upweight'},  
                    'spec': upw_specdict
                    }
            true_cat_corr = {
                    'catalog': upw_catdict, 
                    'correction': {'name': 'true'},  
                    'spec': upw_specdict
                    }
            tophat_conv_cat_corr = {
                    'catalog': upw_catdict, 
                    'correction': tophat_conv_corrdict, 
                    'spec': upw_specdict
                    }
            upw_avg_spec = AvgSpec(n_mocks[i_ell], 'pk', upw_cat_corr)
            upw_avg_spec.read()
            
            true_avg_spec = AvgSpec(n_mocks[i_ell], 'pk', true_cat_corr)
            true_avg_spec.read()
            
            tophat_conv_avg_spec = AvgSpec(n_mocks[i_ell], 'pk', tophat_conv_cat_corr)
            tophat_conv_avg_spec.read()
        
            k_arr = upw_avg_spec.k
            upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
            true_avg_pk = getattr(true_avg_spec, ''.join(['p', str(ell), 'k']))
            tophat_conv_avg_pk = getattr(tophat_conv_avg_spec, ''.join(['p', str(ell), 'k']))
            
            # Del P / P (normalized sample variance of upweighted power spectrum) 
            pk_err = upw_avg_spec.stddev()

            #if ell == 0: 
            #    mono_cv = sub.plot(k_arr, np.abs(pk_err/upw_avg_pk), lw=1, ls='--', color='gray')
            #else: 
            #    quad_cv = sub.plot(k_arr, np.abs(pk_err/upw_avg_pk), lw=1, ls='--', color='gray')
            if ell == 0: 
                plt_label = None
            elif ell == 2: 
                plt_label = cat_label

            mock = sub.scatter(
                    k_arr, 
                    np.abs((upw_avg_pk - tophat_conv_avg_pk)/pk_err),
                    lw=0, color=cat_color
                    )
            sub.plot(
                    k_arr, 
                    np.abs((upw_avg_pk - true_avg_pk)/pk_err),
                    lw=2, color='k', ls='--'
                    )
            sub.plot((10**-3, 10**0), (0.0,0.0), 'k--')

            # x-axis 
            sub.set_xscale('log')
            sub.set_xlim([10**-3,10**0])
            sub.set_xlabel('k (h/Mpc)', fontsize=24)
            # y-axis 
            if ell == 0: 
                #sub.set_xticklabels([])
                #sub.set_ylim([-0.05,0.25])
                #mono_label = [r"$\mathtt{|\frac{\Delta P^{upw}_l}{P^{upw}_l}(k)|}$"]
                #sub.legend(mono_cv, mono_label, loc='upper right', scatterpoints=1, 
                #        prop={'size':25}, handletextpad=0.2, borderpad=0.5) 
                sub.set_ylabel(
                        r"$\mathtt{|1-\overline{P_{"+str(ell)+"}^{tophat}}/\overline{P_{"+str(ell)+"}^{upw}}|}$", 
                        fontsize=24)

            elif ell == 2: 
                #sub.set_ylim([-0.1,1.0])
                #sub.yaxis.tick_right()
                sub.set_ylabel(
                        r"$\mathtt{|1-\overline{P_{"+str(ell)+"}^{tophat}}/\overline{P_{"+str(ell)+"}^{upw}}|}$", 
                        fontsize=24, labelpad=20)
                sub.yaxis.set_label_position('right')

            fig_data.append({
                'ell': ell,
                'name': 'nseries', 
                'k': k_arr,
                'upw_Pk': upw_avg_pk, 
                'tophat_conv_Pk': tophat_conv_avg_pk, 
                'norm_resid': np.abs(1.- tophat_conv_avg_pk/upw_avg_pk),
                'norm_upw_cv': np.abs(pk_err/upw_avg_pk)
                })

        #fig.subplots_adjust(wspace=0.0, hspace=0.0)

        fig_file = 'figure/fc_paper/mock_catalog_tophatconv_upw_Plk.png' 
        fig_data_file = 'figure/fc_paper/mock_catalog_tophatconv_upw_Plk.figdata.p' 
        pickle.dump(fig_data, open(fig_data_file, 'wb'))
        fig.savefig(fig_file, bbox_inches="tight", dpi=150)
        plt.close()
"""


if __name__ == "__main__": 
    #fig_chisquared_NN_p02k_mock(0, catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], Ngrid=960, rebin='log', use='nseries')
    #fig_chisquared_NN_p02k_mock(2, catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[70, 70, 1], Ngrid=960, rebin='log', use='nseries')
    #chisquared_p02k_mock_SNoffset(0, corr='NN', catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], Ngrid=960, rebin='log', use='nseries')
    #for ell in [0, 2]: 
    #    print 'l = ', ell 
        #print 'NN', chisquared_p02k_mock(ell, corr='NN', catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], Ngrid=960, rebin='log', use='nseries')
        #print 'NN', chisquared_p02k_mock(ell, corr='NN', catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], Ngrid=960, rebin='log', use='qpm')
        #print 'LOS Recon', chisquared_p02k_mock(ell, corr='LR', catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], Ngrid=960, rebin='log', use='nseries')
        #print 'LOS Recon', chisquared_p02k_mock(ell, corr='LR', catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], Ngrid=960, rebin='log', use='qpm')
    #for ell in [0]:
    #    print 'Florian', chisquared_p02k_mock(ell, corr='FB', catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], Ngrid=960, rebin='log', use='qpm')
        #print 'Hector ', chisquared_p02k_mock(ell, corr='GM', catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], Ngrid=960, rebin='log', use='qpm')
    #print table_dlospeak_fit(['nseries', 'qpm', 'bigmd', 'cmass'], fit='gauss')
    #fig_2pcf_tophat(n_mocks=range(1,6), scale='5x5')
    #fig_p0k_p2k_mock_catalogs(catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84,100,1], Ngrid=[960, 480], rebin=6)
    #fig_nn_upw_p0k_p2k_mock_catalogs(catalogs=['nseries', 'qpm', 'bigmd'], 
    #        n_mocks=[84,100,1], Ngrid=[960, 480])
    #fig_NN_p02k_mocks(catalogs=['nseries','qpm', 'bigmd'], n_mocks=[84,100,1], Ngrid=[960, 480])
    #fig_dlospeakknown_p02k_mock_resid(10, Ngrid=960)
    #fig_NN_p02k_mocks_resid(catalogs=['nseries', 'qpm', 'bigmd'], 
    #        n_mocks=[84,100,1], Ngrid=[960, 960], rebin=6)
    #fig_NN_p024k_mocks_resid(catalogs=['nseries', 'qpm', 'bigmd'], 
    #        n_mocks=[84,100,1], Ngrid=[960, 960, 960], rebin=6)
    #fig_NN_p0k_mock_norm_resid(catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], 
    #        Ngrid=960, rebin=6, chi2_use='nseries')
    fig_dlospeak_p02k_mock_resid(catalogs=['nseries', 'qpm', 'bigmd'], 
            n_mocks=[84,100,1], Ngrid=[960, 960], florian_offset=250., rebin=6)
    #fig_dlospeak_p0k_mock_norm_resid(catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1], 
    #        Ngrid=960, rebin=6, chi2_use='nseries')
    #fig_tophatconv_delP(n_mocks=[84, 84], Ngrid=[960, 960], figdata=False, k_rebin=6)
    #fig_tophatconv_delP_untrusted(0.3, n_mocks=[84, 84], Ngrid=[960, 960], order=18, k_rebin='log')
    #fig_tophatconv_p0k_p2k_mock_catalogs(n_mocks=[84, 84], Ngrid=[960, 480])
    #fibcol_fraction(catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[5,5,1])
    #fig_mock_catalogs(catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[5,5,1])
    #fig_mock_dlos(catalogs=['nseries', 'qpm', 'bigmd'], n_mocks=[84, 100, 1])
