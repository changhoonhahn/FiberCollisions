'''

Plot figures for paper
FiberCollisions Project

Author(s): ChangHoon Hahn


'''


import numpy as np 
import os.path
import subprocess
import cosmolopy as cosmos
from matplotlib.collections import LineCollection

# --- Local --- 
import fibcol_data as fc_data
import fibcol_nbar as fc_nbar
import fibcol_spec as fc_spec
import plot_fibcol as fc_plot
import fibcol_utility as fc_util
import fibercollisions as fc
import galaxy_environment as genv



def plot_z_dist_fcpaper(cat_corrs):
    '''
    Plot z distribution for ldg, qpm, and tm
    '''
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    for icatcorr, cat_corr in enumerate(cat_corrs): 
        catalog = cat_corr['catalog']
    
        if catalog['name'].lower() == 'cmass': 
            # compute nbar_ngal using nbar file 
            x_mid, x_low, x_high, corr_ngal  = np.loadtxt('/mount/riachuelo1/hahn/data/nbar-cmass-dr12v4-N-Reid.dat', unpack=True, usecols=[0,1,2,6])
        else: 
            corr_ngal_file = fc_nbar.get_nbar_ngal_file('random', **cat_corr)
            if os.path.isfile(corr_ngal_file) == False: 
                print 'Constructing ', corr_ngal_file
                fc_nbar.write_nbar_ngal('random', **cat_corr)

            x_mid, x_low, x_high, corr_ngal = np.loadtxt(corr_ngal_file, unpack=True, usecols=[0,1,2,3]) 
        total_ngal = np.sum(corr_ngal) 

        norm_corr_ngal = np.array([np.float(corr_ngal[i])/np.float(total_ngal) for i in range(len(corr_ngal))])

        if catalog['name'].lower() == 'lasdamasgeo': 
            cat_label = 'Las Damas'
            cat_color = pretty_colors[1]
        elif catalog['name'].lower() == 'qpm': 
            cat_label = 'QPM' 
            cat_color = pretty_colors[3]
        elif catalog['name'].lower() == 'tilingmock': 
            cat_label = 'Tiling Mock' 
            cat_color = pretty_colors[5]
        elif catalog['name'].lower() == 'cmass': 
            cat_label = 'BOSS CMASS'
            cat_color = pretty_colors[0]
        else: 
            raise NameError('asdf') 

        sub.step(x_mid, norm_corr_ngal, lw=4, color=cat_color, label=cat_label, where='mid') 

    sub.set_xlabel(r"Redshift $(z)$", fontsize=24) 
    sub.set_xlim([0.0, 0.8])
    sub.set_ylim([0.0, 0.04])
    sub.legend(loc='upper left', scatterpoints=1, prop={'size':20}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+'fcpaper_z_dist.png', bbox_inches="tight")
    fig.clear() 

def plot_pk_fcpaper(catalogs): 
    '''
    Comparison of true P(k)avg for all mocks
    '''
    ldg_nmock = 1 
    qpm_nmock = 5
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(7, 8))
    sub = fig.add_subplot(111)

    for i_cat, catalog in enumerate(catalogs): 
        catalog = {'name': catalog}    # catalog dict

        # default power/bispectrum box settings 
        if catalog['name'].lower() in ('lasdamasgeo', 'qpm', 'cmass'): 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
        elif catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 
        else: 
            raise NameError('not coded yet') 

        corr_methods = [{'name': 'true'}]

        for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
            # for LasDamasGeo 
            if catalog['name'].lower() == 'lasdamasgeo': 
                n_file = 0
                for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                    for letter in ['a', 'b', 'c', 'd']: 
                        # set catalog correction dictionary for specific file 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_catalog['letter'] = letter
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                        # import power spectrum 
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        try: 
                            avg_k 
                        except NameError: 
                            avg_k = power_i.k
                            sum_Pk = power_i.Pk
                        else: 
                            sum_Pk = sum_Pk + power_i.Pk

                        n_file = n_file+1

            # Tiling Mocks and CMASS
            elif catalog['name'].lower() in ('tilingmock', 'cmass'): 
                i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
                # read tiling mock P(k)
                power = fc_spec.Spec('power', **i_cat_corr)
                power.readfile()

                avg_k = power.k
                sum_Pk = power.Pk
                n_file = 1          # only one file 

            elif catalog['name'].lower() == 'qpm': 
                n_file = 0  
                for i_mock in range(1, qpm_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                        
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_Pk = power_i.Pk
                    else: 
                        sum_Pk = sum_Pk + power_i.Pk

                    n_file = n_file+1
            else: 
                raise NameError('not yet coded!')
            print n_file 
            avg_Pk = np.array([sum_Pk[i]/np.float(n_file) for i in range(len(sum_Pk))])    # average P(k)
        
            Pk_var = np.zeros(len(avg_k)) 
            if catalog['name'].lower() == 'lasdamasgeo': 
                for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                    for letter in ['a', 'b', 'c', 'd']: 
                        # set catalog correction dictionary for specific file 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_catalog['letter'] = letter
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        Pk_i = power_i.Pk
                        
                        # calculate sample variance
                        Pk_var = Pk_var + (avg_Pk - Pk_i)**2

            elif catalog['name'].lower() == 'qpm': 
                for i_mock in range(1, qpm_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                    
                    Pk_i = power_i.Pk
                    
                    Pk_var = Pk_var + (avg_Pk - Pk_i)**2

            if catalog['name'].lower() == 'lasdamasgeo': 
                cat_label = 'Las Damas'
                cat_color = pretty_colors[1]
            elif catalog['name'].lower() == 'qpm': 
                cat_label = 'QPM'
                cat_color = pretty_colors[3]
            elif catalog['name'].lower() == 'tilingmock': 
                cat_label = 'Tiling Mock'
                cat_color = pretty_colors[5]
            elif catalog['name'].lower() == 'cmass': 
                cat_label = 'BOSS DR12'
                cat_color = pretty_colors[0]
            
            if catalog['name'].lower() in ('lasdamasgeo', 'qpm'): 
                Pk_var = np.sqrt(Pk_var/np.float(n_file))
                sub.fill_between(avg_k, avg_Pk-Pk_var, avg_Pk+Pk_var, facecolor=cat_color, edgecolor=cat_color)
                sub.plot(avg_k, [avg_Pk[i]/10.**10 for i in range(len(avg_Pk))], lw=2, color=cat_color, label=cat_label)
                #sub.errorbar(avg_k, avg_Pk, yerr=Pk_var, color=pretty_colors[i_cat])
            else: 
                sub.plot(avg_k, avg_Pk, lw=3, color=cat_color, label=cat_label)

            del avg_k
            del avg_Pk

    sub.set_xscale('log')
    sub.set_yscale('log')
    sub.set_xlim([10**-3,10**0])
    sub.set_ylim([10**3, 10**5.5])

    sub.set_ylabel(r"$\mathtt{\overline{P(k)}}$", fontsize=20)
    sub.legend(loc='lower left', scatterpoints=1, prop={'size':20})
    sub.set_xlabel('k (h/Mpc)', fontsize=24)
    
    fig_name = ''.join(['fcpaper_pk_comp.png'])     
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]), bbox_inches="tight")
    fig.clear()

def plot_pk_upw_fcpaper(catalogs): 
    ''' Comparison of P(k)avg for NN-upweight method with sample variance of P(k) mocks. 

    Parameters
    ----------
    * catalogs : list of catalogs 

    '''

    ldg_nmock = 40 
    qpm_nmock = 100
    nseries_nmock = 84 
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(5*len(catalogs), 5))

    for i_cat, catalog in enumerate(catalogs): 
        sub = fig.add_subplot(1, len(catalogs), i_cat+1) 

        catalog = {'name': catalog}    # catalog dict

        # default power/bispectrum box settings 
        if catalog['name'].lower() in ('lasdamasgeo', 'qpm', 'nseries'): 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
        elif catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 
        else: 
            raise NameError('not coded yet') 

        corr_methods = [{'name': 'true'}, {'name': 'upweight'}]

        for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
            # for LasDamasGeo 
            if catalog['name'].lower() == 'lasdamasgeo': 
                n_file = 0
                for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                    for letter in ['a', 'b', 'c', 'd']: 
                        # set catalog correction dictionary for specific file 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_catalog['letter'] = letter
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                        # import power spectrum 
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        try: 
                            avg_k 
                        except NameError: 
                            avg_k = power_i.k
                            sum_Pk = power_i.Pk
                        else: 
                            sum_Pk = sum_Pk + power_i.Pk

                        n_file = n_file+1

            # Tiling Mocks
            elif catalog['name'].lower() == 'tilingmock': 
                i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
                # read tiling mock P(k)
                power = fc_spec.Spec('power', **i_cat_corr)
                power.readfile()

                avg_k = power.k
                sum_Pk = power.Pk
                n_file = 1          # only one file 

            elif catalog['name'].lower() == 'qpm': 
                n_file = 0  
                for i_mock in range(1, qpm_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                        
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_Pk = power_i.Pk
                    else: 
                        sum_Pk = sum_Pk + power_i.Pk

                    n_file = n_file+1

            elif catalog['name'].lower() == 'nseries': 
                n_file = 0  
                for i_mock in range(1, nseries_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                        
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_Pk = power_i.Pk
                    else: 
                        sum_Pk = sum_Pk + power_i.Pk

                    n_file = n_file+1
            else: 
                raise NameError('not yet coded!')

            avg_Pk = [sum_Pk[i]/np.float(n_file) for i in range(len(sum_Pk))]    # average P(k)
        
            # P(k) residual comparison
            if correction['name'].lower() == 'true':
                # Should be first and the demoninator otherwise code crashes 
                avg_Pk_true = avg_Pk

                Pk_var = np.zeros(len(avg_k)) 
                if catalog['name'].lower() == 'lasdamasgeo': 
                    for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                        for letter in ['a', 'b', 'c', 'd']: 
                            # set catalog correction dictionary for specific file 
                            i_catalog = catalog.copy() 
                            i_catalog['n_mock'] = i_mock
                            i_catalog['letter'] = letter
                            i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                            power_i = fc_spec.Spec('power', **i_cat_corr)
                            power_i.readfile()
                            
                            Pk_i = power_i.Pk
                            
                            # calculate sample variance
                            Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2

                elif catalog['name'].lower() == 'qpm': 
                    for i_mock in range(1, qpm_nmock+1): 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        Pk_i = power_i.Pk
                        
                        Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2
                
                elif catalog['name'].lower() == 'nseries': 
                    for i_mock in range(1, nseries_nmock+1): 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        Pk_i = power_i.Pk
                        
                        Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2
                
                if catalog['name'].lower() != 'tilingmock': 
                    Pk_var = np.sqrt(Pk_var/np.float(n_file))
                    print Pk_var/avg_Pk_true
                    sub.plot(avg_k, Pk_var/avg_Pk_true+1.0, lw=2, ls='--', color=pretty_colors[0], label=r'$\mathtt{\Delta P(k)/\overline{P_{true}(k)}}$')
                else: 
                    sub.plot(avg_k, Pk_var/avg_Pk_true, lw=4, color=pretty_colors[0], label=r'$\mathtt{\Delta P/P}$')
            else: 
                if correction['name'].lower() == 'upweight': 
                    corr_label = 'NN-upweight'
                    corr_color = 1
                elif correction['name'].lower() == 'peaknbar':
                    corr_label = 'd_{LOS}-peak'
                    corr_color = 3 

                # set residual label
                resid_label = ''.join([r"$\mathtt{", corr_label,"}$"])
                
                # plot residual 
                sub.scatter(avg_k, fc_plot.residual(avg_Pk, avg_Pk_true), color=pretty_colors[corr_color], label=resid_label)
                #print resid_label 
                #print residual(avg_Pk, avg_Pk_true)

            sub.plot((10**-3, 10**0), (1.0,1.0), 'k--')          # draw horizontal line? 
            
            del avg_k
            del avg_Pk

        # set axes
        sub.text(0.002, 0.95, ''.join([str(n_file), ' ', catalog['name'].upper()]))              # number of mocks + Catalog name 
        sub.set_xscale('log')
        sub.set_xlim([10**-3,10**0])
        sub.set_ylim([0.85, 1.15]) 

        if i_cat == 0 : 
            sub.set_ylabel(r"$\mathtt{\overline{P_{\rm{corr}}(k)}/\overline{P_{\rm{true}}(k)}}$", fontsize=20)
        elif i_cat == len(catalogs)-1: 
            sub.set_yticklabels([])
            sub.legend(loc='lower right', scatterpoints=1, prop={'size':18})
        else: 
            sub.set_yticklabels([])
            sub.set_xlabel('k (h/Mpc)', fontsize=24)
            #sub.axes.get_yaxis().set_ticks([])
    
    fig_name = ''.join(['fcpaper_pk_upw_comp.png'])     
    fig.savefig(''.join(['figure/', fig_name]), bbox_inches="tight")
    fig.clear()

def plot_pk_peakonly_fcpaper(catalogs): 
    '''
    Comparison of P(k)avg for peaknbar fc correction method with sample variance of P(k) mocks. 
    '''
    ldg_nmock = 40 
    qpm_nmock = 100
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(15, 5))

    for i_cat, catalog in enumerate(catalogs): 
        sub = fig.add_subplot(1, 3, i_cat+1) 

        catalog = {'name': catalog}    # catalog dict

        # default power/bispectrum box settings 
        if catalog['name'].lower() in ('lasdamasgeo', 'qpm'): 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
        elif catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 
        else: 
            raise NameError('not coded yet') 

        # dlos distribution sigma values 
        if catalog['name'].lower() == 'lasdamasgeo': 
            sigma = 6.5 
        elif catalog['name'].lower() == 'qpm': 
            sigma = 4.4 
        elif catalog['name'].lower() == 'tilingmock': 
            sigma = 4.8 
        else: 
            raise NameError('not coded yet') 

        corr_methods = [{'name': 'true'}, {'name': 'upweight'}, {'name': 'peaknbar', 'sigma':sigma, 'fpeak':1.0, 'fit': 'gauss'}]

        for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
            # for LasDamasGeo 
            if catalog['name'].lower() == 'lasdamasgeo': 
                n_file = 0
                for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                    for letter in ['a', 'b', 'c', 'd']: 
                        # set catalog correction dictionary for specific file 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_catalog['letter'] = letter
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                        # import power spectrum 
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        try: 
                            avg_k 
                        except NameError: 
                            avg_k = power_i.k
                            sum_Pk = power_i.Pk
                        else: 
                            sum_Pk = sum_Pk + power_i.Pk

                        n_file = n_file+1

            # Tiling Mocks
            elif catalog['name'].lower() == 'tilingmock': 
                i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
                # read tiling mock P(k)
                power = fc_spec.Spec('power', **i_cat_corr)
                power.readfile()

                avg_k = power.k
                sum_Pk = power.Pk
                n_file = 1          # only one file 

            elif catalog['name'].lower() == 'qpm': 
                n_file = 0  
                for i_mock in range(1, qpm_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                        
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_Pk = power_i.Pk
                    else: 
                        sum_Pk = sum_Pk + power_i.Pk

                    n_file = n_file+1
            else: 
                raise NameError('not yet coded!')

            avg_Pk = [sum_Pk[i]/np.float(n_file) for i in range(len(sum_Pk))]    # average P(k)
        
            # P(k) residual comparison
            if correction['name'].lower() == 'true':
                # Should be first and the demoninator otherwise code crashes 
                avg_Pk_true = avg_Pk

                Pk_var = np.zeros(len(avg_k)) 
                if catalog['name'].lower() == 'lasdamasgeo': 
                    for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                        for letter in ['a', 'b', 'c', 'd']: 
                            # set catalog correction dictionary for specific file 
                            i_catalog = catalog.copy() 
                            i_catalog['n_mock'] = i_mock
                            i_catalog['letter'] = letter
                            i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                            power_i = fc_spec.Spec('power', **i_cat_corr)
                            power_i.readfile()
                            
                            Pk_i = power_i.Pk
                            
                            # calculate sample variance
                            Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2

                elif catalog['name'].lower() == 'qpm': 
                    for i_mock in range(1, qpm_nmock+1): 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        Pk_i = power_i.Pk
                        
                        Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2
                
                if catalog['name'].lower() != 'tilingmock': 
                    Pk_var = np.sqrt(Pk_var/np.float(n_file))
                    print Pk_var/avg_Pk_true
                    sub.plot(avg_k, Pk_var/avg_Pk_true+1.0, lw=2, ls='--', color=pretty_colors[0], label=r'$\mathtt{\Delta P(k)/\overline{P_{true}(k)}}$')
                else: 
                    sub.plot(avg_k, Pk_var/avg_Pk_true, lw=4, color=pretty_colors[0], label=r'$\mathtt{\Delta P/P}$')
            else: 
                if correction['name'].lower() == 'upweight': 
                    corr_label = 'NN'
                    corr_color = 1
                elif correction['name'].lower() == 'peaknbar':
                    corr_label = 'd_{LOS}-peak'
                    corr_color = 3 

                # set residual label
                resid_label = ''.join([r"$\mathtt{", "\overline{P(k)_{", corr_label,"}}/\overline{P(k)_{true}}}$"])
                
                # plot residual 
                sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=pretty_colors[corr_color], label=resid_label)
                #print resid_label 
                #print residual(avg_Pk, avg_Pk_true)

            sub.plot((10**-3, 10**0), (1.0,1.0), 'k--')          # draw horizontal line? 
            
            del avg_k
            del avg_Pk

        # set axes
        sub.text(0.01, 0.9, ''.join([str(n_file), ' ', catalog['name'].upper()]))              # number of mocks + Catalog name 
        sub.set_xscale('log')
        sub.set_xlim([10**-3,10**0])
        sub.set_ylim([0.85, 1.15]) 

        if i_cat == 0 : 
            sub.set_ylabel(r"$\mathtt{\overline{P_{\rm{corr}}(k)}/\overline{P_{\rm{true}}(k)}}$", fontsize=20)
        elif i_cat == 2: 
            sub.set_yticklabels([])
            sub.legend(loc='lower right', scatterpoints=1, prop={'size':18})
        else: 
            sub.set_yticklabels([])
            sub.set_xlabel('k (h/Mpc)', fontsize=24)
            #sub.axes.get_yaxis().set_ticks([])
    
    fig_name = ''.join(['fcpaper_pk_peakonly_comp.png'])     
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]), bbox_inches="tight")
    fig.clear()

def plot_pk_shotnoiseonly_fcpaper(catalogs): 
    '''
    Comparison of P(k)avg for peaknbar fc correction method with sample variance of P(k) mocks. 
    '''
    ldg_nmock = 40 
    qpm_nmock = 100 
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(15, 5))

    for i_cat, catalogue in enumerate(catalogs): 
        sub = fig.add_subplot(1, 3, i_cat+1) 

        catalog = {'name': catalogue}    # catalog dict

        # default power/bispectrum box settings 
        if catalog['name'].lower() in ('lasdamasgeo', 'qpm'): 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 'quad': False} 
        elif catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 
        else: 
            raise NameError('not coded yet') 
    
        corr_methods = [{'name': 'true'}, {'name': 'upweight'}, {'name': 'shotnoise'}, {'name':'floriansn'}, {'name':'hectorsn'}]

        for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
            # for LasDamasGeo 
            if catalog['name'].lower() == 'lasdamasgeo': 
                n_file = 0
                for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                    for letter in ['a', 'b', 'c', 'd']: 
                        # set catalog correction dictionary for specific file 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_catalog['letter'] = letter
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        try: 
                            avg_k 
                        except NameError: 
                            avg_k = power_i.k
                            sum_Pk = power_i.Pk
                        else: 
                            sum_Pk = sum_Pk + power_i.Pk

                        n_file = n_file+1

            # Tiling Mocks
            elif catalog['name'].lower() == 'tilingmock': 
                i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
                # read tiling mock P(k)
                power = fc_spec.Spec('power', **i_cat_corr)
                power.readfile()

                avg_k = power.k
                sum_Pk = power.Pk
                n_file = 1          # only one file 

            elif catalog['name'].lower() == 'qpm': 
                n_file = 0  
                for i_mock in range(1, qpm_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                        
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_Pk = power_i.Pk
                    else: 
                        sum_Pk = sum_Pk + power_i.Pk

                    n_file = n_file+1
            else: 
                raise NameError('not yet coded!')

            avg_Pk = [sum_Pk[i]/np.float(n_file) for i in range(len(sum_Pk))]    # average P(k)
        
            # P(k) residual comparison
            if correction['name'].lower() == 'true':
                # Should be first and the demoninator otherwise code crashes 
                avg_Pk_true = avg_Pk

                Pk_var = np.zeros(len(avg_k)) 
                if catalog['name'].lower() == 'lasdamasgeo': 
                    for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                        for letter in ['a', 'b', 'c', 'd']: 
                            # set catalog correction dictionary for specific file 
                            i_catalog = catalog.copy() 
                            i_catalog['n_mock'] = i_mock
                            i_catalog['letter'] = letter
                            i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                            power_i = fc_spec.Spec('power', **i_cat_corr)
                            power_i.readfile()
                            
                            Pk_i = power_i.Pk
                            
                            # calculate sample variance
                            Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2

                elif catalog['name'].lower() == 'qpm': 
                    for i_mock in range(1, qpm_nmock+1): 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        Pk_i = power_i.Pk
                        
                        Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2
                
                if catalog['name'].lower() != 'tilingmock': 
                    Pk_var = np.sqrt(Pk_var/np.float(n_file))
                    print Pk_var/avg_Pk_true
                    if i_cat == 2: 
                        sub.plot(avg_k, Pk_var/avg_Pk_true+1.0, lw=2, ls='--', color=pretty_colors[0], label=r'$\mathtt{\Delta P(k)/\overline{P_{true}(k)}}$')
                    else: 
                        sub.plot(avg_k, Pk_var/avg_Pk_true+1.0, lw=2, ls='--', color=pretty_colors[0])
                else: 
                    sub.plot(avg_k, Pk_var/avg_Pk_true, lw=4, color=pretty_colors[0], label=r'$\mathtt{\Delta P(k)/\overline{P_{true}(k)}}$')
            else: 
                if correction['name'].lower() == 'upweight': 
                    corr_color = 1
                elif correction['name'].lower() == 'shotnoise': 
                    corr_color = 5 
                elif correction['name'].lower() == 'floriansn': 
                    corr_color = 7 
                elif correction['name'].lower() == 'hectorsn':
                    corr_color = 9

                # plot residual 
                if i_cat == 1: 
                    if correction['name'].lower() in ('upweight', 'shotnoise'): 
                        if correction['name'].lower() == 'upweight': 
                            corr_label = 'Eq.10'
                        elif correction['name'].lower() == 'shotnoise': 
                            corr_label = 'Eq.15'
                        sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=pretty_colors[corr_color], label=corr_label)
                    else: 
                        sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=pretty_colors[corr_color])

                    sub.legend(loc='lower left', scatterpoints=1, prop={'size':18})
                elif i_cat == 2: 
                    if correction['name'].lower() in ('floriansn', 'hectorsn'):
                        if correction['name'].lower() == 'floriansn': 
                            corr_label = 'Beutler+2014'
                        elif correction['name'].lower() == 'hectorsn':
                            corr_label = 'Gil-Marin+2014'
                        sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=pretty_colors[corr_color], label=corr_label)
                    else: 
                        sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=pretty_colors[corr_color])
                    sub.legend(loc='lower left', scatterpoints=1, prop={'size':18})
                else:
                    sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=pretty_colors[corr_color])

                sub.plot((10**-3, 10**0), (1.0,1.0), 'k--')          # draw horizontal line? 

                #print resid_label 
                #print residual(avg_Pk, avg_Pk_true)

            if catalog['name'].lower() == 'lasdamasgeo': 
                cat_label = 'Las Damas'
            elif catalog['name'].lower() == 'qpm': 
                cat_label = 'QPM' 
            elif catalog['name'].lower() == 'tilingmock': 
                cat_label = 'Tiling Mock' 
            else: 
                raise NameError('asdf') 
            sub.set_title(cat_label) 

            del avg_k
            del avg_Pk

        # set axes
        ylimit = [0.85,1.15] 
        
        #sub.text(1.5*10**-3, np.mean(ylimit), ''.join([str(n_file), ' ', catalog['name'].upper()]))              # number of mocks + Catalog name 
        sub.set_xscale('log')
        sub.set_xlim([10**-3,10**0])
        sub.set_ylim([0.85, 1.15]) 

        if i_cat == 0 : 
            sub.set_ylabel(r"$\mathtt{\overline{P_{\rm{corr}}(k)}/\overline{P_{\rm{true}}(k)}}$", fontsize=20)
        elif i_cat == 2: 
            sub.set_yticklabels([])
        else: 
            sub.set_yticklabels([])
            sub.set_xlabel('k (h/Mpc)', fontsize=24)
            #sub.axes.get_yaxis().set_ticks([])
    
    fig_name = ''.join(['fcpaper_pk_shotnoiseonly_comp.png'])     
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]), bbox_inches="tight")
    fig.clear()

def plot_pk_mpfit_peakshotnoise_fcpaper(catalogs): 
    '''
    Comparison of P(k)avg for peak+shotnoise fc correction method using MPFit sigma and fpeak with sample variance of P(k) mocks. 
    '''
    ldg_nmock = 39 
    qpm_nmock = 100 
    nseries_nmock = 19 
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(15, 5))

    for i_cat, catalogue in enumerate(catalogs): 
        sub = fig.add_subplot(1, 3, i_cat+1) 

        catalog = {'name': catalogue}    # catalog dict
        '''
        # default power/bispectrum box settings 
        if catalog['name'].lower() in ('lasdamasgeo', 'qpm', 'nseries'): 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':960} 
        elif catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':960} 
        else: 
            raise NameError('not coded yet') 
        '''
            
        spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':960} 

        # Gaussian fit
        if catalog['name'].lower() == 'lasdamasgeo': 
            corr = {'name': 'peakshot', 'sigma': 6.5, 'fpeak': 0.76, 'fit':'gauss'}
        elif catalog['name'].lower() == 'qpm': 
            corr = {'name': 'peakshot', 'sigma': 4.4, 'fpeak':0.65, 'fit':'gauss'}
        elif catalog['name'].lower() == 'nseries': 
            corr = {'name': 'peakshot', 'sigma': 4.0, 'fpeak':0.68, 'fit':'gauss'}
        elif catalog['name'].lower() == 'tilingmock': 
            corr = {'name': 'peakshot', 'sigma': 4.8, 'fpeak':0.62, 'fit':'gauss'}
        else: 
            raise NameError('not coded yet') 

        corr_methods = [{'name': 'true'}, corr]

        for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
            # for LasDamasGeo 
            if catalog['name'].lower() == 'lasdamasgeo': 
                n_file = 0
                for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                    for letter in ['a', 'b', 'c', 'd']: 
                        # set catalog correction dictionary for specific file 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_catalog['letter'] = letter
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                        if correction['name'] != 'Upweight Igal_Irand': 
                            # import power spectrum 
                            power_i = fc_spec.Spec('power', **i_cat_corr)
                        else: 
                            i_cat_corr = {'catalog': i_catalog, 'correction': {'name':'upweight'}, 'spec':spec}
                            power_i = fc_spec.Spec('power', Igal_Irand==True, **i_cat_corr)
                        power_i.readfile()
                        
                        try: 
                            avg_k 
                        except NameError: 
                            avg_k = power_i.k
                            sum_Pk = power_i.Pk
                        else: 
                            sum_Pk = sum_Pk + power_i.Pk

                        n_file = n_file+1

            # Tiling Mocks
            elif catalog['name'].lower() == 'tilingmock': 
                i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
                # read tiling mock P(k)
                power = fc_spec.Spec('power', **i_cat_corr)
                power.readfile()

                avg_k = power.k
                sum_Pk = power.Pk
                n_file = 1          # only one file 

            elif catalog['name'].lower() == 'qpm': 
                n_file = 0  
                for i_mock in range(1, qpm_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                        
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_Pk = power_i.Pk
                    else: 
                        sum_Pk = sum_Pk + power_i.Pk

                    n_file = n_file+1

            elif catalog['name'].lower() == 'nseries': 
                n_file = 0  
                for i_mock in range(1, nseries_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                        
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_Pk = power_i.Pk
                    else: 
                        sum_Pk = sum_Pk + power_i.Pk

                    n_file = n_file+1
            else: 
                raise NameError('not yet coded!')

            avg_Pk = [sum_Pk[i]/np.float(n_file) for i in range(len(sum_Pk))]    # average P(k)
        
            # P(k) residual comparison
            if correction['name'].lower() == 'true':
                # Should be first and the demoninator otherwise code crashes 
                avg_Pk_true = avg_Pk

                Pk_var = np.zeros(len(avg_k)) 
                if catalog['name'].lower() == 'lasdamasgeo': 
                    for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                        for letter in ['a', 'b', 'c', 'd']: 
                            # set catalog correction dictionary for specific file 
                            i_catalog = catalog.copy() 
                            i_catalog['n_mock'] = i_mock
                            i_catalog['letter'] = letter
                            i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                            power_i = fc_spec.Spec('power', **i_cat_corr)
                            power_i.readfile()
                            
                            Pk_i = power_i.Pk
                            
                            # calculate sample variance
                            Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2

                elif catalog['name'].lower() == 'qpm': 
                    for i_mock in range(1, qpm_nmock+1): 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        Pk_i = power_i.Pk
                        
                        Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2
                
                elif catalog['name'].lower() == 'nseries': 
                    for i_mock in range(1, nseries_nmock+1): 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        Pk_i = power_i.Pk
                        
                        Pk_var = Pk_var + (avg_Pk_true - Pk_i)**2
                
                if catalog['name'].lower() != 'tilingmock': 
                    Pk_var = np.sqrt(Pk_var/np.float(n_file))
                    #print Pk_var/avg_Pk_true
                    sub.plot(avg_k, Pk_var/avg_Pk_true+1.0, lw=2, ls='--', 
                            color=pretty_colors[0], label=r'$\Delta P/P$')
                else: 
                    sub.plot(avg_k, Pk_var/avg_Pk_true, lw=4, color=pretty_colors[0], label='Sample Variance')
            else: 
                corr_label = 'Hahn et al.'
                corr_color = 7  

                # set residual label
                #resid_label = ''.join([r"$\mathtt{", "\overline{P(k)_{", corr_label,"}}/\overline{P(k)_{True}}}$"])
                
                # plot residual 
                sub.scatter(avg_k, fc_plot.residual(avg_Pk, avg_Pk_true), color=pretty_colors[corr_color], label=corr_label)
                #print resid_label 
                print catalog['name'].lower() 
                print np.min(fc_plot.residual(avg_Pk, avg_Pk_true)), \
                        np.max(fc_plot.residual(avg_Pk, avg_Pk_true))
    

            sub.plot((10**-3, 10**0), (1.0, 1.0), 'k--') 
            del avg_k
            del avg_Pk

        # set axes
        ylimit = [0.95,1.15] 
        
        #sub.text(1.5*10**-3, np.mean(ylimit), ''.join([str(n_file), ' ', catalog['name'].upper()]))              # number of mocks + Catalog name 
        sub.set_xscale('log')
        sub.set_xlim([10**-3,10**0])
        sub.set_ylim([0.85, 1.15]) 
        
        if catalog['name'].lower() == 'lasdamasgeo': 
            cat_label = 'Las Damas'
        elif catalog['name'].lower() == 'qpm': 
            cat_label = 'QPM' 
        elif catalog['name'].lower() == 'nseries': 
            cat_label = 'N Series' 
        elif catalog['name'].lower() == 'tilingmock': 
            cat_label = 'Tiling Mock' 
        else: 
            raise NameError('asdf') 
        sub.set_title(cat_label) 

        if i_cat == 0 : 
            sub.set_ylabel(r"$\mathtt{\overline{P_{\rm{corr}}(k)}/\overline{P_{\rm{true}}(k)}}$", fontsize=20)
        elif i_cat == 2: 
            sub.set_yticklabels([])
            sub.legend(loc='lower right', scatterpoints=1, prop={'size':18})
        else: 
            sub.set_yticklabels([])
            sub.set_xlabel('k (h/Mpc)', fontsize=24)
            #sub.axes.get_yaxis().set_ticks([])
    
    fig_name = ''.join(['fcpaper_pk_peakshotnoise_mpfit_comp.png'])     
    fig.savefig(''.join(['figure/', fig_name]), bbox_inches="tight")
    fig.clear()

def plot_p2k_fcpaper(catalogs): 
    ''' Comparison of P2(k)avg (powerspectrum quadrupole) for peak+shotnoise fc correction method 
    with sample variance of P(k) mocks. 

    '''
    ldg_nmock = 40 
    qpm_nmock = 100 
    nseries_nmock = 9 
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(15, 5))

    for i_cat, catalogue in enumerate(catalogs): 
        sub = fig.add_subplot(1, 3, i_cat+1) 

        catalog = {'name': catalogue}    # catalog dict
        '''
        # default power/bispectrum box settings 
        if catalog['name'].lower() in ('lasdamasgeo', 'qpm', 'nseries'): 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':960} 
        '''
        if catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 
        else: 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 

        # Gaussian fit
        if catalog['name'].lower() == 'lasdamasgeo': 
            corr = {'name': 'peakshot', 'sigma': 6.5, 'fpeak': 0.76, 'fit':'gauss'}
        elif catalog['name'].lower() == 'qpm': 
            corr = {'name': 'peakshot', 'sigma': 4.4, 'fpeak':0.65, 'fit':'gauss'}
        elif catalog['name'].lower() == 'nseries': 
            corr = {'name': 'peakshot', 'sigma': 4.0, 'fpeak':0.7, 'fit':'gauss'}
        elif catalog['name'].lower() == 'tilingmock': 
            corr = {'name': 'peakshot', 'sigma': 4.8, 'fpeak':0.63, 'fit':'gauss'}
        else: 
            raise NameError('not coded yet') 

        spec['quad'] = True 

        corr_methods = [{'name': 'true'}, {'name': 'upweight'}, corr]

        for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
            # for LasDamasGeo 
            if catalog['name'].lower() == 'lasdamasgeo': 
                n_file = 0
                for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                    for letter in ['a', 'b', 'c', 'd']: 
                        # set catalog correction dictionary for specific file 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_catalog['letter'] = letter
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                        if correction['name'] != 'Upweight Igal_Irand': 
                            # import power spectrum 
                            power_i = fc_spec.Spec('power', **i_cat_corr)
                        else: 
                            i_cat_corr = {'catalog': i_catalog, 'correction': {'name':'upweight'}, 'spec':spec}
                            power_i = fc_spec.Spec('power', Igal_Irand==True, **i_cat_corr)
                        power_i.readfile()
                        
                        try: 
                            avg_k 
                        except NameError: 
                            avg_k = power_i.k
                            sum_P2k = power_i.P2k
                        else: 
                            sum_P2k = sum_P2k + power_i.P2k

                        n_file = n_file+1

            # Tiling Mocks
            elif catalog['name'].lower() == 'tilingmock': 
                i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
                # read tiling mock P(k)
                power = fc_spec.Spec('power', **i_cat_corr)
                power.readfile()

                avg_k = power.k
                sum_P2k = power.P2k
                n_file = 1          # only one file 

            elif catalog['name'].lower() == 'qpm': 
                n_file = 0  
                for i_mock in range(1, qpm_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                        
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_P2k = power_i.P2k
                    else: 
                        sum_P2k = sum_P2k + power_i.P2k

                    n_file = n_file+1

            elif catalog['name'].lower() == 'nseries': 
                n_file = 0  
                for i_mock in range(1, nseries_nmock+1): 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                        
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    power_i.readfile()
                        
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_P2k = power_i.P2k
                    else: 
                        sum_P2k = sum_P2k + power_i.P2k

                    n_file = n_file+1
            else: 
                raise NameError('not yet coded!')

            avg_P2k = [sum_P2k[i]/np.float(n_file) for i in range(len(sum_P2k))]    # average P(k)
        
            # P(k) residual comparison
            if correction['name'].lower() == 'true':
                # Should be first and the demoninator otherwise code crashes 
                avg_P2k_true = avg_P2k

                Pk_var = np.zeros(len(avg_k)) 
                if catalog['name'].lower() == 'lasdamasgeo': 
                    for i_mock in range(1, ldg_nmock+1):                       # compute average[P(k)] for each correction method
                        for letter in ['a', 'b', 'c', 'd']: 
                            # set catalog correction dictionary for specific file 
                            i_catalog = catalog.copy() 
                            i_catalog['n_mock'] = i_mock
                            i_catalog['letter'] = letter
                            i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                            power_i = fc_spec.Spec('power', **i_cat_corr)
                            power_i.readfile()
                            
                            Pk_i = power_i.P2k
                            
                            # calculate sample variance
                            Pk_var = Pk_var + (avg_P2k_true - Pk_i)**2

                elif catalog['name'].lower() == 'qpm': 
                    for i_mock in range(1, qpm_nmock+1): 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        Pk_i = power_i.P2k
                        
                        Pk_var = Pk_var + (avg_P2k_true - Pk_i)**2
                
                elif catalog['name'].lower() == 'nseries': 
                    for i_mock in range(1, nseries_nmock+1): 
                        i_catalog = catalog.copy() 
                        i_catalog['n_mock'] = i_mock
                        i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                            
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                        power_i.readfile()
                        
                        Pk_i = power_i.P2k
                        
                        Pk_var = Pk_var + (avg_P2k_true - Pk_i)**2
                
                if catalog['name'].lower() != 'tilingmock': 

                    Pk_var = np.sqrt(Pk_var/np.float(n_file))
                    #print Pk_var/avg_Pk_true
                    sub.plot(avg_k, Pk_var/avg_P2k_true+1.0, lw=2, ls='--', 
                            color=pretty_colors[0], label=r'$\Delta P/P$')
                else: 
                    sub.plot(avg_k, Pk_var/avg_P2k_true, lw=4, color=pretty_colors[0], label='Sample Variance')
            else: 
                if correction['name'].lower() == 'upweight':
                    corr_label = 'NN-upweight'
                    corr_color = 1
                else: 
                    corr_label = 'Hahn et al.'
                    corr_color = 7  

                # set residual label
                #resid_label = ''.join([r"$\mathtt{", "\overline{P(k)_{", corr_label,"}}/\overline{P(k)_{True}}}$"])
                
                # plot residual 
                sub.scatter(avg_k, fc_plot.residual(avg_P2k, avg_P2k_true), 
                        color=pretty_colors[corr_color], label=corr_label)
                #print resid_label 
                print catalog['name'].lower() 
                print np.min(fc_plot.residual(avg_P2k, avg_P2k_true)), \
                        np.max(fc_plot.residual(avg_P2k, avg_P2k_true))
    

            sub.plot((10**-3, 10**0), (1.0, 1.0), 'k--') 
            del avg_k
            del avg_P2k

        # set axes
        ylimit = [0.95,1.15] 
        
        #sub.text(1.5*10**-3, np.mean(ylimit), ''.join([str(n_file), ' ', catalog['name'].upper()]))              # number of mocks + Catalog name 
        sub.set_xscale('log')
        sub.set_xlim([2*10**-2,10**0])
        sub.set_ylim([0.5, 1.8]) 
        
        if catalog['name'].lower() == 'lasdamasgeo': 
            cat_label = 'Las Damas'
        elif catalog['name'].lower() == 'qpm': 
            cat_label = 'QPM' 
        elif catalog['name'].lower() == 'nseries': 
            cat_label = 'N Series' 
        elif catalog['name'].lower() == 'tilingmock': 
            cat_label = 'Tiling Mock' 
        else: 
            raise NameError('asdf') 
        sub.set_title(cat_label) 

        if i_cat == 0 : 
            sub.set_ylabel(r"$\mathtt{\overline{P_2^{\rm{corr}}(k)}/\overline{P_2^{\rm{true}}(k)}}$", fontsize=20)
        elif i_cat == 1: 
            sub.set_yticklabels([])
            sub.legend(loc='lower left', scatterpoints=1, prop={'size':18})
        else: 
            sub.set_yticklabels([])
            sub.set_xlabel('k (h/Mpc)', fontsize=24)
            #sub.axes.get_yaxis().set_ticks([])
    
    fig_name = ''.join(['fcpaper_p2k.png'])     
    fig.savefig(''.join(['figure/', fig_name]), bbox_inches="tight")
    fig.clear()

def plot_chi2_fcpaper(catalogs): 
    '''
    chi-squared comparison for mock catalogues, correction scheme specified 
    '''
    qpm_n_mock = 100 
    ldg_n_mock = 40 
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(15, 5))

    for i_cat, catalogue in enumerate(catalogs): 
        sub = fig.add_subplot(1, 3, i_cat+1)

        catalog = {'name': catalogue}    # catalog dict
        # hardcoded default power/bispectrum box settings 
        if catalog['name'].lower() == 'lasdamasgeo': 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
            corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma':6.5, 'fpeak':0.76, 'fit':'gauss'}] 
            n_mock = 40 
        elif catalog['name'].lower() == 'qpm': 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
            corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}]
            n_mock = 100
        elif catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 
            corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma':4.8, 'fpeak':0.62, 'fit':'gauss'}]
            n_mock = 1 
        else: 
            raise NameError('not coded yet') 

        for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
            Pks = []
            
            if catalog['name'].lower() == 'tilingmock': 
                n_mock = 1
            elif catalog['name'].lower() == 'qpm': 
                n_mock = qpm_n_mock
            elif catalog['name'].lower() == 'lasdamasgeo': 
                n_mock = ldg_n_mock 
            else: 
                raise NameError('asdfklasdf') 

            cat_corr = {'catalog': catalog, 'correction': correction, 'spec': spec}  

            avg_Pk_data = fc.avgPk(n_mock, **cat_corr) 
            avg_k = avg_Pk_data[0]
            avg_Pk = avg_Pk_data[1]
            
            # P(k) Comparison
            if correction['name'].lower() == 'true':        
                avg_Pk_true = avg_Pk

                # get delta P
                delta_P = fc.deltaP(n_mock, **cat_corr) 
                delP = delta_P[1]
            else: 
                #print 'DELTA P', delP[-5:-1]
                #print 'Pk_corr', correction['name'].lower(),' - Pk_true', avg_Pk[-5:-1]-avg_Pk_true[-5:-1]

                chi2_total = (avg_Pk - avg_Pk_true)**2/delP**2
               
                # chi2(kmax) 
                chi2 = np.zeros(len(avg_k))
                for i_k, kval in enumerate(avg_k): 
                    chi2[i_k] = np.sum(chi2_total[0:i_k])/np.float(len(chi2_total[0:i_k]))

                if correction['name'].lower() == 'floriansn': 
                    dot_size = 8 
                else: 
                    dot_size = 24 

                if i_cat == 1: 
                    if correction['name'].lower() == 'upweight': 
                        corr_label = "Eq. 10"
                    elif correction['name'].lower() == 'floriansn': 
                        corr_label = "Beutler+2014" 
                    else:
                        corr_label = ''
                    sub.scatter(avg_k, chi2, color=pretty_colors[i_corr], s=dot_size, label=corr_label)
                elif i_cat == 2: 
                    if correction['name'].lower() in ('peakshot', 'vlospeakshot'): 
                        corr_label = "Hahn et al."
                    elif correction['name'].lower() == 'hectorsn': 
                        corr_label = "Gil-Marin+2014"
                    else:
                        corr_label = ''

                    sub.scatter(avg_k, chi2, color=pretty_colors[i_corr], s=dot_size, label=corr_label)
                else: 
                    sub.scatter(avg_k, chi2, color=pretty_colors[i_corr], s=dot_size)
                
                print avg_k[-5:-1]
                print chi2[-5:-1], catalog['name'], correction['name']

        # set axes
        #if catalog['name'].lower() == 'tilingmock': ylimit = [0.0,1.0] 
        #elif catalog['name'].lower() == 'qpm': ylimit = [0.0, 0.3] 
        #elif catalog['name'].lower() == 'lasdamasgeo': ylimit = [0.0, 0.2] 
        ylimit = [10**-4, 10**2] 
        if catalog['name'].lower() == 'lasdamasgeo': 
            cat_label = 'Las Damas'
        elif catalog['name'].lower() == 'qpm': 
            cat_label = 'QPM' 
        elif catalog['name'].lower() == 'tilingmock': 
            cat_label = 'Tiling Mock' 
        else: 
            raise NameError('asdf') 
        ytext = 1.15
        #ylabel = r"$\mathtt{\frac{1}{N[k<k_{max}]}\sum\limits_{k< k_{max}} (\overline{P(k)}-\overline{P(k)_{\rm{True}}})^2/(\Delta P)^2}$"
        ylabel = r"$\chi^2 (\mathtt{k_{max}})$"

        fig_name = ''.join(['fcpaper_pk_chisquared_comparison.png'])     

        #sub.text(1.5*10**-3, np.mean(ylimit), ''.join([str(n_file), ' ', catalog_name.upper()]))              # number of mocks + Catalog name 
        sub.set_title(cat_label) 
        sub.set_xscale('log')
        sub.set_yscale('log')
        sub.set_xlim([10**-3,10**0])
        sub.set_ylim(ylimit)
        if i_cat == 0: 
            sub.set_ylabel(ylabel, fontsize=20)
        elif i_cat == 2: 
            sub.set_yticklabels([])
            sub.legend(loc='upper left', scatterpoints=1, prop={'size':18})
        else: 
            sub.legend(loc='upper left', scatterpoints=1, prop={'size':18})
            sub.set_yticklabels([])
            sub.set_xlabel(r'$\mathtt{k_{max}}$ (h/Mpc)', fontsize=24)
    
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]), bbox_inches="tight")
    fig.clear()

def plot_fcpaper(): 
    # ------------------------------------------------------------
    # figures for FCPAPER
    plot_pk_peakonly_fcpaper(['lasdamasgeo', 'qpm', 'tilingmock'])
    plot_pk_shotnoiseonly_fcpaper(['lasdamasgeo', 'qpm', 'tilingmock'])
    plot_pk_mpfit_peakshotnoise_fcpaper(['lasdamasgeo', 'qpm', 'tilingmock'])
    plot_chi2_fcpaper(['lasdamasgeo', 'qpm', 'tilingmock'])
    # ---------------------------------------------------------------

if __name__=='__main__': 
    catalogs = ['lasdamasgeo', 'nseries', 'tilingmock']
    #plot_pk_upw_fcpaper(catalogs)
    #plot_pk_mpfit_peakshotnoise_fcpaper(catalogs)
    plot_p2k_fcpaper(['lasdamasgeo', 'nseries'])

