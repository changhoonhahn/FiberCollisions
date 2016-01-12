'''

Figures for the Fiber Collisions paper

'''
import numpy as np 
import scipy as sp
import os.path
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages

# --- Local --- 
from dlos.dlos import Dlos
from corr_spec.corr_corrdata import CorrCorrData
from corr_spec.corr_average import CorrAvgSpec as AvgSpec

# --- plotting ---
from defutility.plotting import prettyplot 
from defutility.plotting import prettycolors 

def fig_mock_catalogs(catalogs=['nseries'], n_mocks=[84]): 
    '''
    Figure that compares the redshift distribution of mock catalogs 

    Parameters
    ----------
    catalogs : 
        List of strings that specify the mock catalogs to include. (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(5,4))
    sub = fig.add_subplot(111)
    
    # CMASS DR12v4
    z_mid, z_low, z_high, cmass_ngal_eff = np.loadtxt(
            '/mount/riachuelo1/hahn/data/nbar-cmass-dr12v4-N-Reid.dat', 
            unpack=True, 
            usecols=[0,1,2,6]
            )

    for i_cat, cat in enumerate(catalogs):
        for i_mock in range(1, n_mocks[i_cat]+1): 
            if cat == 'qpm': 
                cat_label = 'QPM' 
                cat_color = pretty_colors[3]
            elif cat == 'tilingmock': 
                cat_label = 'Tiling Mock' 
                cat_color = pretty_colors[5]
            elif cat == 'nseries': 
                cat_label = 'Nseries'
                cat_color = pretty_colors[7]
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

        
    sub.step(
            z_mid, 
            cmass_ngal_eff / np.sum(cmass_ngal_eff), 
            lw=4, 
            color=pretty_colors[0], 
            label='BOSS CMASS', 
            where='mid') 

    sub.set_xlabel(r"Redshift $(z)$", fontsize=24) 
    sub.set_xlim([0.3, 0.85])
    sub.set_ylim([0.0, 0.04])
    sub.legend(loc='upper left', scatterpoints=1, prop={'size':20}) 
    
    fig.savefig(
            ''.join([
                'figure/fc_paper/', 
                'mock_catalog_z_dist.png'
                ]), 
            bbox_inches="tight")
    #fig.clear() 

def fig_p0k_p2k_mock_catalogs(catalogs=['nseries'], n_mocks=[84], Ngrid=[960]): 
    '''
    Figure that compares the redshift distribution of mock catalogs 

    Parameters
    ----------
    catalogs : 
        List of strings that specify the mock catalogs to include. (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14,8))
    mono = fig.add_subplot(121)
    quad = fig.add_subplot(122)

    for i_cat, cat in enumerate(catalogs):
        if cat == 'qpm': 
            cat_label = 'QPM' 
            cat_color = pretty_colors[3]
        elif cat == 'tilingmock': 
            cat_label = 'Tiling Mock' 
            cat_color = pretty_colors[7]
        elif cat == 'nseries': 
            cat_label = 'Nseries'
            cat_color = pretty_colors[1]
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
            avg_spec.read()
        
            k_arr = avg_spec.k
            avg_pk = getattr(avg_spec, ''.join(['p', str(ell), 'k']))
            pk_err = avg_spec.stddev()
            pk_err_lower = pk_err 
            pk_err_lower[np.abs(pk_err) > np.abs(avg_pk)] = avg_pk[np.abs(pk_err) > np.abs(avg_pk)] * 0.99
            
            if ell == 0: 
                mono.fill_between(k_arr, avg_pk - pk_err_lower, avg_pk + pk_err, 
                     color=cat_color, alpha=0.8)
                mono.plot(k_arr, np.zeros(len(avg_pk)), lw=4, 
                        color=cat_color, label=None)
            elif ell == 2: 
                quad.fill_between(k_arr, np.abs(avg_pk - pk_err_lower), np.abs(avg_pk + pk_err),
                     color=cat_color, alpha=0.8)
                quad.plot(k_arr, np.zeros(len(avg_pk)), lw=4, 
                        color=cat_color, label=cat_label)
    
    # Plot CMASS P0(k) and P2(k)
    #sub.step(
    #        z_mid, 
    #        cmass_ngal_eff / np.sum(cmass_ngal_eff), 
    #        lw=4, 
    #        color=pretty_colors[0], 
    #        label='BOSS CMASS', 
    #        where='mid') 

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
    quad.legend(loc='upper left', scatterpoints=1, prop={'size':20}) 
    
    plt.sca(mono)
    plt.xticks([10**-3, 10**-2, 10**-1])
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig.savefig(
            ''.join([
                'figure/fc_paper/', 
                'mock_catalog_Plk.png'
                ]), 
            bbox_inches="tight")
    plt.close()

def fig_nn_upw_p0k_p2k_mock_catalogs(catalogs=['nseries'], n_mocks=[84], Ngrid=960): 
    '''
    Figure that compares P_l(k)^nn / P_l(k)^true for l = 0, 2 for mock catalogs  

    Parameters
    ----------
    catalogs : 
        List of strings that specify the mock catalogs to include. 
        (e.g. ['nseries', 'qpm', 'bigmultidark', 'tilingmock'])
    n_mocks : 
        
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14,8))

    if len(catalogs) < 3: 
        for i in xrange(3 - len(catalogs)): 
            catalogs.append(catalogs[0])
            n_mocks.append(n_mocks[0])
    
    for i_ell, ell in enumerate([0, 2]):  # monopole or quadrupole
        #sub = fig.add_subplot(1, 2, i_ell+1)
        for i_cat, cat in enumerate(catalogs):
            sub = fig.add_subplot(2, 3, 3*i_ell + i_cat+1)

            if cat == 'qpm': 
                cat_label = 'QPM' 
                cat_color = pretty_colors[3]
            elif cat == 'tilingmock': 
                cat_label = 'Tiling Mock' 
                cat_color = pretty_colors[5]
            elif cat == 'nseries': 
                cat_label = 'Nseries'
                cat_color = pretty_colors[7]
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
            tru_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', tru_cat_corr)
            tru_avg_spec.read()
            
            upw_avg_spec = AvgSpec(n_mocks[i_cat], 'pk', upw_cat_corr)
            upw_avg_spec.read()
        
            k_arr = tru_avg_spec.k
            tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
            upw_avg_pk = getattr(upw_avg_spec, ''.join(['p', str(ell), 'k']))
            pk_err = tru_avg_spec.stddev()
            
            if ell == 0: 
                plt_label = None
                delP_label = r"$\mathtt{1 - \frac{\Delta P_0}{P_0}(k)}$"
            elif ell == 2: 
                plt_label = cat_label
                delP_label = r"$\mathtt{1 + \frac{\Delta P_2}{P_2}(k)}$"

            if ell == 0: 
                sub.plot(
                        k_arr, 1. - (pk_err/tru_avg_pk), 
                        lw=1, ls='--', color='gray', label=delP_label)
            else: 
                sub.plot(
                        k_arr, 1. + pk_err/tru_avg_pk, 
                        lw=1, ls='--', color='gray', label=delP_label)

            sub.scatter(k_arr, upw_avg_pk/tru_avg_pk, lw=0, color=cat_color, label=plt_label)
            sub.plot((10**-3, 10**0), (1.0,1.0), 'k--')

            # Axes 
            sub.set_xlim([10**-3,10**0])
            sub.set_xscale('log')

            if ell == 0: 
                sub.set_xticklabels([])
                sub.set_ylim([0.65,1.025])
                #plt.sca(sub)
                #plt.yticks(list(np.arange(0.7, 1.1, 0.1)) )
                #plt.xticks([10**-3, 10**-2, 10**-1])
                if i_cat == len(catalogs)-1: 
                    sub.legend(loc='lower left', scatterpoints=1, prop={'size':20}) 
            elif ell == 2: 
                sub.set_ylim([0.9,2.0])
                #sub.set_yscale('log')

                #sub.legend(loc='upper left', scatterpoints=1, prop={'size':20}) 
                #plt.sca(sub)
                #plt.yticks(list(np.arange(1.0, 2.2, 0.2)) )
                
                if i_cat == len(catalogs)-1: 
                    sub.legend(loc='upper left', scatterpoints=1, prop={'size':20}) 
                elif i_cat == 0: 
                    plt.sca(sub)
                    plt.yticks(list(np.arange(1.0, 2.0, 0.2)) )
                    plt.xticks([10**-3, 10**-2, 10**-1])
                else: 
                    sub.set_xlabel('k (h/Mpc)', fontsize=24)
                    plt.sca(sub)
                    plt.xticks([10**-3, 10**-2, 10**-1])
                #sub.yaxis.tick_right()
                #sub.yaxis.set_label_position('right')

            if i_cat == 0:
                sub.set_ylabel(r"$\mathtt{\overline{P_{"+str(ell)+"}^{NN}}/\overline{P_{"+str(ell)+"}^{true}}(k)}$", fontsize=24)
                #sub.set_xlabel('k (h/Mpc)', fontsize=24)
            else: 
                sub.set_yticklabels([])

    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig.savefig(
            ''.join([
                'figure/fc_paper/', 
                'mock_catalog_NN_true_Plk.png'
                ]), 
            bbox_inches="tight")
    plt.close()

def fig_mock_dlos(catalogs=['nseries'], n_mocks=[84], peak_range=[-15.0, 15.0]): 
    '''
    dLOS distribution plot for fc_paper
    '''
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=[14, 5])
    sub = fig.add_subplot(111)

    catalogs = ['cmass'] + catalogs
    n_mocks = [1] + n_mocks

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
        
        dlos_mid, dlos_dist = cat_deelos.dlos_dist(binsize = fd_binsize, normed=True)

        if cat == 'qpm': 
            cat_label = 'QPM' 
            cat_color = pretty_colors[3]
        elif cat == 'tilingmock': 
            cat_label = 'Tiling Mock' 
            cat_color = pretty_colors[5]
        elif cat == 'nseries': 
            cat_label = 'Nseries'
            cat_color = pretty_colors[7]
        elif cat == 'cmass': 
            cat_label = 'CMASS'
            cat_color = pretty_colors[0]
        else: 
            raise NotImplementedError 

        sub.plot(dlos_mid, dlos_dist, lw=4, color=cat_color, label=cat_label)
    
    sub.set_xlim([-45.0, 45.0])
    sub.set_xlabel('$\mathtt{d_{LOS}}$ (Mpc)')
    sub.set_ylabel('$\mathtt{p(d_{LOS})}$')
    sub.legend(loc='upper left')

    fig.savefig(
            ''.join([
                'figure/fc_paper/', 
                'mock_catalog_dlos.png'
                ]), 
            bbox_inches="tight")
    plt.close()

def fig_2pcf_tophat(n_mocks=5, scale='small'): 
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
        contour_range = np.arange(-0.1, 0.11, 0.01)
        #contour_range = 20
        n_rp = 100 
        n_pi = 100
    
    fig = plt.figure(figsize=(15,10))
    sub = fig.add_subplot(111)

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

        sub.vlines(0.43, 0.0, np.max(r_p), lw=8, linestyle='--', color='k')

        sub.set_xlabel('$\mathtt{r_{p} \; (Mpc/h)}$', fontsize=40)
        sub.set_ylabel('$\pi \; \mathtt{(Mpc/h)}$', fontsize=40)
        sub.set_xlim([np.min(rp_bins), np.max(rp_bins)])
        sub.set_ylim([np.min(pi_bins), np.max(pi_bins)])
        
        #sub.set_title(r"$1 - (1 + \xi^\mathtt{"+corr.upper()+r"})/(1+ \xi^\mathtt{TRUE})$", fontsize=40)
    
        fig_name = ''.join([
                'figure/fc_paper/', 
                '2pcf_Nseries_', corr, '_', str(len(n_mocks)), 'mocks.', scale, '.tophat.png'
                ])
        fig.savefig(fig_name, bbox_inches="tight")
        plt.close()


if __name__ == "__main__": 
    fig_2pcf_tophat(n_mocks=range(1,6), scale='5x5')
    #fig_nn_upw_p0k_p2k_mock_catalogs(catalogs=['nseries', 'qpm'], n_mocks=[10,10], Ngrid=[960, 480])
    #fig_p0k_p2k_mock_catalogs(catalogs=['nseries', 'qpm'], n_mocks=[5,5], Ngrid=[960, 480])
    #fig_mock_catalogs(catalogs=['nseries'], n_mocks=[5])
    #fig_mock_dlos(catalogs=['nseries'], n_mocks=[5])
