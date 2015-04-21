import numpy as np 
import os.path
import subprocess
import cosmolopy as cosmos
import fibcol_data as fc_data
import fibcol_nbar as fc_nbar
import fibcol_spec as fc_spec
import fibcol_utility as fc_util
import fibercollisions as fc
import galaxy_environment as genv
from matplotlib.collections import LineCollection

# Plotting -----------------------------------------------------------------
def plot_pk_fibcol_comp(catalog_name, n_mock, corr_methods, resid='False', quad=False): 
    '''
    Comparison of P(k)avg for different fibercollision correciton methods
    Make sure 'true' is first in the correction method list
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(7, 8))
    sub = fig.add_subplot(111)

    catalog = {'name': catalog_name}    # catalog dict
    # hardcoded default power/bispectrum box settings 
    if catalog_name.lower() in ('lasdamasgeo', 'qpm'): 
        spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':960, 'quad': quad} 
    elif catalog_name.lower() == 'tilingmock': 
        spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':960} 
    else: 
        raise NameError('not coded yet') 

    for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
        # for LasDamasGeo 
        if catalog_name.lower() == 'lasdamasgeo': 
            n_file = 0
            for i_mock in range(1,n_mock+1):                       # compute average[P(k)] for each correction method
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
        elif catalog_name.lower() == 'tilingmock': 
            i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
            # read tiling mock P(k)
            power = fc_spec.Spec('power', **i_cat_corr)
            power.readfile()

            avg_k = power.k
            sum_Pk = power.Pk
            n_file = 1          # only one file 

        elif catalog_name.lower() == 'qpm': 
            n_file = 0  
            for i_mock in range(1, n_mock+1): 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock
                i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                    
                power_i = fc_spec.Spec('power', **i_cat_corr)
                #print power_i.file_name
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
        
        # P(k) Comparison
        if resid == 'False':  
            if correction['name'].lower() == 'true':
                sub.plot(avg_k, avg_Pk, color=pretty_colors[i_corr], label=r"$\mathtt{\overline{P(k)_{"+correction['name']+"}}}$")

            elif (correction['name'].lower() == 'peak') or (correction['name'].lower() == 'peaknbar') or \
                    (correction['name'].lower() == 'peaktest') or (correction['name'].lower() == 'peakshot') or (correction['name'].lower() == 'allpeakshot'):
                # scatter plot
                sub.scatter(avg_k, avg_Pk, color=pretty_colors[i_corr], label=r"$\mathtt{\overline{P(k)_{"+correction['name']+','+\
                        ''.join([str(correction['sigma']), str(correction['fpeak'])])+"}}}$")
            else: 
                sub.scatter(avg_k, avg_Pk, color=pretty_colors[i_corr], label=r"$\mathtt{\overline{P(k)_{"+correction['name']+"}}}$")

        # P(k) residual comparison
        else:               
            if correction['name'].lower() == 'true':
                # Should be first and the demoninator otherwise code crashes 
                avg_Pk_true = avg_Pk
            else: 
                if correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'allpeakshot', 'vlospeakshot'):
                    if correction['fit'].lower() in ('gauss', 'expon'):
                        # set residual label
                        resid_label = ''.join([r"$\mathtt{", "\overline{P(k)_{", correction['fit'], ", ", correction['name'], ",", 
                            ''.join([str(correction['sigma']), str(correction['fpeak'])]), "}}/\overline{P(k)_{True}}}$"])
                    elif correction['fit'].lower() in ('true'): 
                        # set residual label
                        resid_label = ''.join([r"$\mathtt{", "\overline{P(k)_{", correction['fit'], ", ", correction['name'], ",", 
                            str(correction['fpeak']), "}}/\overline{P(k)_{True}}}$"])
                    else: 
                        raise NameError('asdflkjasdf')
                    
                    # plot residual 
                    sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=pretty_colors[i_corr], label=resid_label)
                    print resid_label 

                    print avg_k[(avg_k > 0.15) & (avg_k < 0.2)]
                    print residual(avg_Pk, avg_Pk_true)[(avg_k > 0.15) & (avg_k < 0.2)]
                    #print avg_k[-5:-1]
                    #print residual(avg_Pk, avg_Pk_true)[-5:-1]
                else:
                    resid_label = ''.join([r"$\mathtt{", "\overline{P(k)_{",correction['name'], "}}/\overline{P(k)_{True}}}$"])

                    sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=pretty_colors[i_corr], label=resid_label)
                    print resid_label 
                    print avg_k[(avg_k > 0.15) & (avg_k < 0.2)]
                    print residual(avg_Pk, avg_Pk_true)[(avg_k > 0.15) & (avg_k < 0.2)]
                    #print avg_k[-5:-1]
                    #print residual(avg_Pk, avg_Pk_true)[-5:-1]

                #chi2 = np.sum((residual(avg_Pk, avg_Pk_true)-1.0)**2) 
                #print chi2
        try: 
            corr_str 
        except NameError: 
            corr_str = ''.join(['_', correction['name']]) 
        else: 
            if correction['name'] in ('peak', 'peaknbar', 'peaktest', 'peakshot'): 
                if correction['fit'].lower() in ('expon', 'gauss'): 
                    corr_str = ''.join([corr_str, '_', correction['fit'], '_', correction['name'], 
                        '_sigma', str(correction['sigma']), '_fpeak', str(correction['fpeak'])]) 
                elif correction['fit'].lower() in ('true'): 
                    corr_str = ''.join([corr_str, '_', correction['fit'], '_', correction['name'], 
                        '_fpeak', str(correction['fpeak'])]) 
                else: 
                    raise NameError('asdlfkjasdf')
            else: 
                corr_str = ''.join([corr_str, '_', correction['name']]) 

        del avg_k
        del avg_Pk

    # set axes
    if resid == 'False':
        #        if catalog_name.lower() == 'lasdamasgeo':
        ylimit = [10**2,10**5.5]
        ytext = 10**5.2

        ylabel = 'P(k)'
        sub.set_yscale('log')
        resid_str = ''
        fig_name = ''.join(['powerspec_', catalog_name.lower(), 
            '_fibcoll_', corr_str, '_comparison.png'])
    else:
        ylimit = [0.85,1.15] 
        ytext = 1.15
        ylabel = r"$\mathtt{\overline{P(k)}/\overline{P(k)_{\rm{True}}}}$"
        resid_str = '_residual'

    if quad == True: 
        quad_tag = 'quad_'
    else: 
        quad_tag = ''

    fig_name = ''.join(['powerspec_', quad_tag, catalog_name.lower(), 
        '_fibcoll_', corr_str, resid_str, '_comparison.png'])     

    sub.text(1.5*10**-3, np.mean(ylimit), ''.join([str(n_file), ' ', catalog_name.upper()]))              # number of mocks + Catalog name 
    sub.set_xscale('log')
    sub.set_xlim([10**-3,10**0])
    sub.set_ylim(ylimit)
    sub.set_xlabel('k', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)

    if resid == 'False': 
        sub.legend(loc='lower left', scatterpoints=1, prop={'size':14})
    elif resid == 'True': 
        sub.legend(loc='upper left', scatterpoints=1, prop={'size':14})
    
    try:
        fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]), bbox_inches="tight")
    except IOError: 
        fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', 'powerspec_', catalog_name.lower(), '_fibcoll_toolong', resid_str, '_comparison.png']), bbox_inches="tight")

    fig.clear()

def plot_p2k_fibcol_comp(catalog_name, n_mock, corr_methods, resid='False'): 
    '''
    Comparison of quadrupole P2(k)avg for different fibercollision correciton methods
    Make sure 'true' is first in the correction method list
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(7, 8))
    sub = fig.add_subplot(111)

    catalog = {'name': catalog_name}    # catalog dict
    # hardcoded default power/bispectrum box settings 
    if catalog_name.lower() in ('lasdamasgeo', 'qpm'): 
        spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 'quad': True} 
    elif catalog_name.lower() == 'tilingmock': 
        spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360, 'quad': True} 
    else: 
        raise NameError('not coded yet') 

    for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
        # for LasDamasGeo 
        if catalog_name.lower() == 'lasdamasgeo': 
            n_file = 0
            for i_mock in range(1,n_mock+1):                       # compute average[P(k)] for each correction method
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
                        sum_P2k = power_i.P2k
                    else: 
                        sum_P2k = sum_P2k + power_i.P2k

                    n_file = n_file+1

        # Tiling Mocks
        elif catalog_name.lower() == 'tilingmock': 
            i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
            # read tiling mock P(k)
            power = fc_spec.Spec('power', **i_cat_corr)
            power.readfile()

            avg_k = power.k
            sum_P2k = power.P2k
            n_file = 1          # only one file 

        elif catalog_name.lower() == 'qpm': 
            n_file = 0  
            for i_mock in range(1, n_mock+1): 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock
                i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                    
                power_i = fc_spec.Spec('power', **i_cat_corr)
                power_i.readfile()
                    
                try: 
                    avg_k 
                except NameError: 
                    avg_k = power_i.k
                    sum_P2k = np.abs(power_i.P2k)
                else: 
                    sum_P2k = sum_P2k + np.abs(power_i.P2k)

                n_file = n_file+1

        else: 
            raise NameError('not yet coded!')


        avg_P2k = [sum_P2k[i]/np.float(n_file) for i in range(len(sum_P2k))]    # average P(k)
        
        # P2(k) Comparison
        if resid == 'False':  
            if correction['name'].lower() == 'true':
                # plot for true  Pk
                sub.plot(avg_k, avg_P2k, color=pretty_colors[i_corr], label=r"$\mathtt{\overline{P_{2}(k)_{"+correction['name']+"}}}$")

            elif correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'allpeakshot', 'vlospeakshot'):
                
                if correction['fit'].lower() in ('expon', 'gauss'): 
                    resid_label = r"$\mathtt{\overline{P_2(k)_{"+correction['name']+','+''.join([str(correction['sigma']), str(correction['fpeak'])])+"}}}$"
                elif correction['fit'].lower() in ('true'): 
                    resid_label = r"$\mathtt{\overline{P_2(k)_{"+correction['name']+','+str(correction['fpeak'])+"}}}$"
            
                sub.scatter(avg_k, avg_P2k, color=pretty_colors[i_corr], label=resid_label)
            else: 
                sub.scatter(avg_k, avg_P2k, color=pretty_colors[i_corr], label=r"$\mathtt{\overline{P_2(k)_{"+correction['name']+"}}}$")

        # P(k) residual comparison
        else:               
            if correction['name'].lower() == 'true':    # Should be first and the demoninator otherwise code crashes 
                avg_P2k_true = avg_P2k
            else: 
                if correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'allpeakshot', 'vlospeakshot'):
                    # set residual label
                    if correction['fit'].lower() in ('expon', 'gauss'): 
                        resid_label = ''.join([r"$\mathtt{", "\overline{P_2(k)_{", correction['fit'], ", ", correction['name'], ",", 
                            ''.join([str(correction['sigma']), str(correction['fpeak'])]), "}}/\overline{P_2(k)_{True}}}$"])
                    elif correction['fit'].lower() in ('true'): 
                        resid_label = ''.join([r"$\mathtt{", "\overline{P_2(k)_{", correction['fit'], ", ", correction['name'], ",", 
                            str(correction['fpeak']), "}}/\overline{P_2(k)_{True}}}$"])
                    else:
                        raise NameError('asdflkj')
                    
                    # plot residual 
                    sub.scatter(avg_k, residual(avg_P2k, avg_P2k_true), color=pretty_colors[i_corr], label=resid_label)
                    print avg_P2k[-5:-1]
                    print avg_P2k_true[-5:-1]
                    print resid_label 
                    print residual(avg_P2k, avg_P2k_true)[(avg_k > 0.15) & (avg_k < 0.2)]
                    print avg_k[(avg_k > 0.15) & (avg_k < 0.2)]
                else:
                    resid_label = ''.join([r"$\mathtt{", "\overline{P_2(k)_{",correction['name'], "}}/\overline{P_2(k)_{True}}}$"])

                    sub.scatter(avg_k, residual(avg_P2k, avg_P2k_true), color=pretty_colors[i_corr], label=resid_label)
                    print resid_label 
                    print residual(avg_P2k, avg_P2k_true)[(avg_k > 0.15) & (avg_k < 0.2)]
                    print avg_k[(avg_k > 0.15) & (avg_k < 0.2)]

        try: 
            corr_str 
        except NameError: 
            corr_str = ''.join(['_', correction['name']]) 
        else: 
            if correction['name'].lower in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'vlospeakshot'): 
                corr_str = ''.join([corr_str, '_', correction['fit'], '_', correction['name'], 
                    '_sigma', str(correction['sigma']), '_fpeak', str(correction['fpeak'])]) 
            else: 
                corr_str = ''.join([corr_str, '_', correction['name']]) 

        del avg_k
        del avg_P2k

    # set axes
    if resid == 'False':
        ylimit = [10**2,10**5.5]
        ytext = 10**5.2

        ylabel = 'P(k)'
        sub.set_yscale('log')
        resid_str = ''
    else:
        ylimit = [0.0, 2.0] 
        ytext = 1.15
        ylabel = r"$\mathtt{\overline{P_2(k)}/\overline{P_2(k)_{\rm{True}}}}$"
        resid_str = '_residual'


    fig_name = ''.join(['p2k_', catalog_name.lower(), 
        '_fibcoll_', corr_str, resid_str, '_comparison.png'])     

    sub.text(1.5*10**-3, np.mean(ylimit), ''.join([str(n_file), ' ', catalog_name.upper()]))              # number of mocks + Catalog name 
    sub.set_xscale('log')
    sub.set_xlim([10**-3,10**0])
    sub.set_ylim(ylimit)
    sub.set_xlabel('k', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)

    if resid == 'False': 
        sub.legend(loc='lower left', scatterpoints=1, prop={'size':14})
    elif resid == 'True': 
        sub.legend(loc='upper left', scatterpoints=1, prop={'size':14})
    
    try:
        fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]), bbox_inches="tight")
    except IOError: 
        fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', 'p2k_', catalog_name.lower(), '_fibcoll_toolong', resid_str, '_comparison.png']), bbox_inches="tight")

    fig.clear()

def plot_bispec_fibcolcorr_comparison(BorQ='B', x_axis='triangles', triangle='all', 
        catalog='PTHalo', filespec={'version':'v11p0', 'n_mock':range(1,11), 'Nrandom':100}, 
        correction_method=['true', 'delta', 'peak'], corrspec={'sigma':[0, 0, 6.3], 'fpeak':[1.0, 1.0, 1.0]}, resid='False'):
    '''
    Comparison of B(k)avg for different fibercollision correciton methods
    Make sure 'true' is first in the correction method list
    '''
    prettyplot()
    
    fig = plt.figure(1, figsize=(7, 8))
    sub = fig.add_subplot(111)

    corr_color = ['grey', 'blue', 'red']
    for i_corr, correction in enumerate(correction_method):     # loop through correction methods
        for i_mock in filespec['n_mock']:                       # compute average[B(k)] for each correction method
            i_filespec = filespec.copy()
            i_filespec['n_mock'] = i_mock 
            i_corrspec = corrspec.copy() 
            for key in corrspec.keys(): 
                i_corrspec[key] = corrspec[key][i_corr]

            bispec_i = spec(spectrum='bispec', catalog=catalog, filespec=i_filespec, fibcolcorr=correction, corrspec=i_corrspec) 
            bispec_i.readfile()
            if i_mock == filespec['n_mock'][0]: 
                if x_axis == 'triangles':                       # specify x-axis (triangle index, avgk, or max k) 
                    avg_k = bispec_i.i_triangle                 # x-axis is triangles
                elif x_axis == 'avg_k': 
                    avg_k = bispec_i.avgk                       # x-axis is avg(k1,k2,k3)
                elif x_axis == 'max_k': 
                    avg_k = bispec_i.kmax                       # x-axis is max(k1,k2,k3) 
                avg_k1 = bispec_i.k1
                avg_k2 = bispec_i.k2
                avg_k3 = bispec_i.k3

                if BorQ == 'B': 
                    sum_Bk = bispec_i.Bk
                elif BorQ == 'Q': 
                    sum_Bk = bispec_i.Q
            else: 
                if BorQ == 'B': 
                    sum_Bk = sum_Bk + bispec_i.Bk
                elif BorQ == 'Q': 
                    sum_Bk = sum_Bk + bispec_i.Q
        avg_Bk = [ sum_Bk[i]/float(len(filespec['n_mock'])) for i in range(len(sum_Bk))]    # average B(k)

        if triangle == 'all': 
            tri_index = np.array([True for i in range(len(avg_Bk))], dtype=bool)
        else: 
            tri_index = classify_triangles(avg_k1, avg_k2, avg_k3, triangle=triangle)
            print tri_index

        if resid == 'False':                                    # B(k) Comparison
            if correction == 'true':
                sub.scatter(avg_k[tri_index], avg_Bk[tri_index], 
                        color=corr_color[i_corr], label=r"$\mathtt{\overline{"+BorQ+"(k)_{"+correction+"}}}$")
            else:
                sub.scatter(avg_k[tri_index], avg_Bk[tri_index], 
                        color=corr_color[i_corr], label=r"$\mathtt{\overline{"+BorQ+"(k)_{"+correction+"}}}$")
        else:                                                   # B(k) residual comparison
            if correction == 'true':
                avg_Bk_true = avg_Bk
            else:
                sub.scatter(np.array(avg_k)[tri_index], residual(np.array(avg_Bk)[tri_index], np.array(avg_Bk_true)[tri_index]), 
                        color=corr_color[i_corr], 
                        label=''.join([r"$\mathtt{", "\overline{", BorQ, "(k)_{",correction, "}}/\overline{", BorQ, "(k)_{True}}}$"]))
        if correction == 'peak': 
            sigma_flag = ''.join(['_sigma', str(corrspec['sigma'][i_corr]), '_fpeak', str(corrspec['fpeak'][i_corr])]) 

    if resid=='False':
        if catalog.lower() == 'tilingmock':
            ylimit = [10**2,10**5.5]
            ytext = 10**5.2
        elif catalog.lower() == 'pthalo':
            ylimit = [10**2,10**5.5]
            ytext = 10**5.2
        ylabel = '$\mathtt{'+BorQ+'(k)}$'
        sub.set_yscale('log')
        fig_name = ''.join(['bispec_', BorQ, 'k_', x_axis, '_', catalog.lower(), 
            '_fibcoll_', '_'.join(correction_method), sigma_flag, '_', triangle, 'triangles_comparison.png'])
    else:
        ylimit = [0.9,1.2]
        ytext = 1.15
        ylabel = r"$\mathtt{\overline{"+BorQ+"(k)}/\overline{"+BorQ+"(k)_{\rm{True}}}}$"
        fig_name = ''.join(['bispec_', BorQ, 'k_', x_axis, '_', catalog.lower(), 
            '_fibcoll_', '_'.join(correction_method), sigma_flag, '_', triangle, 'triangles_residual_comparison.png'])     
        sub.set_ylim(ylimit)
    sub.text(1.5*10**-3, ytext, 
            ''.join([str(len(filespec['n_mock'])), ' ', catalog.upper()]))              # number of mocks + Catalog name 
    if x_axis == 'triangles': 
        sub.set_xlabel('Triangles', fontsize=20)
        sub.set_xlim([0, 7500])
    elif x_axis == 'avg_k': 
        sub.set_xlabel('avg(k1, k2, k3)', fontsize=20)
        sub.set_xscale('log')
        sub.set_xlim([10**-3,10**0])
    elif x_axis == 'max_k': 
        sub.set_xlabel('max(k1, k2, k3)', fontsize=20)
        sub.set_xscale('log')
        sub.set_xlim([10**-3,10**0])
    sub.set_ylabel(ylabel, fontsize=20)
    sub.legend(loc='lower left', prop={'size':14}, scatterpoints=1)
    sub.grid(True)
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]))
    fig.clear()

def plot_peakcorrection_check(catalog_name, n_mock, corr_methods):
    '''
    Check dLOS generated during peak correction 
    '''
    # set up figure
    prettyplot() 
    pretty_colors = prettycolors()  
    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    if catalog_name.lower() == 'lasdamasgeo': 
        # For Las Damas Geo 
        for i_correct, correction in enumerate(corr_methods): 
            # Loop through correction methods
            for i_mock in range(1,n_mock+1): 
                # Through n_mock
                for letter in ['a', 'b', 'c', 'd']: 
                    # Through letters

                    i_cat = {'name': catalog_name, 'n_mock': i_mock, 'letter': letter}
                    i_corr = correction.copy() 
                    i_cat_corr = {'catalog': i_cat, 'correction': i_corr}
                    i_gal_data = fc_data.galaxy_data('data', readdata=False, **i_cat_corr)
                    i_dlos_file = i_gal_data.file_name+'.sanitycheck'
                    i_dlos = np.loadtxt(i_dlos_file) 

                    try: 
                        combined_dlos 
                    except NameError: 
                        combined_dlos = i_dlos
                    else: 
                        combined_dlos = np.concatenate([combined_dlos, i_dlos]) 

            # calculate histogram
            x_min = -1000.0
            x_max = 1000.0
            n_bins = int((x_max-x_min)/0.5) 
            dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max], normed=True) 
            xlow = mpc_binedges[:-1]
            xhigh = mpc_binedges[1:] 
            xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
            
            # plot LOS displacement 
            sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[i_correct]) 

            del combined_dlos

    sub.set_xlim([-30, 30])
    sub.set_xlabel(r"$d_{LOS}$ Mpc") 
    plt.show() 

def plot_nbar_comparison(n_mock, cat_corr_list, resid='False'): 
    '''
    comparison of nbar(z) values 
    '''

    for cat_corr in cat_corr_list: 
        pass

def plot_pk_chi2_comp(catalog_name, n_mock, corr_methods): 
    '''
    Comparison of P(k)avg for different fibercollision correciton methods
    Make sure 'true' is first in the correction method list
    '''
    qpm_n_mock = 45
    ldg_n_mock = 40 
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(7, 8))
    sub = fig.add_subplot(111)

    catalog = {'name': catalog_name}    # catalog dict
    # hardcoded default power/bispectrum box settings 
    if catalog_name.lower() in ('lasdamasgeo', 'qpm'): 
        spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
    elif catalog_name.lower() == 'tilingmock': 
        spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 
    else: 
        raise NameError('not coded yet') 

    for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
        Pks = []
        # Tiling Mocks
        if catalog_name.lower() == 'tilingmock': 
            i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
            # read tiling mock P(k)
            power = fc_spec.Spec('power', **i_cat_corr)
            power.readfile()

            avg_k = power.k
            sum_Pk = power.Pk
            n_file = 1          # only one file 

        elif catalog_name.lower() == 'qpm': 
            n_file = 0
            for i_mock in range(1, n_mock+1): 
                i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
                i_cat_corr['catalog']['n_mock'] = i_mock 
                # read tiling mock P(k)
                power = fc_spec.Spec('power', **i_cat_corr)
                power.readfile()
                    
                if i_mock == 1: 
                    avg_k = power.k
                    sum_Pk = power.Pk
                else: 
                    sum_Pk = sum_Pk+power.Pk

                n_file = n_file+1

        elif catalog_name.lower() == 'lasdamasgeo': 
            n_file = 0
            for i_mock in range(1, n_mock+1): 
                for letter in ['a', 'b', 'c', 'd']: 
                    i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
                    i_cat_corr['catalog']['n_mock'] = i_mock 
                    i_cat_corr['catalog']['letter'] = letter
                    # read tiling mock P(k)
                    power = fc_spec.Spec('power', **i_cat_corr)
                    power.readfile()
                    
                    if (i_mock == 1) & (letter == 'a'): 
                        avg_k = power.k
                        sum_Pk = power.Pk
                    else: 
                        sum_Pk = sum_Pk+power.Pk

                    n_file = n_file+1

        else: 
            raise NameError('not yet coded!')

        avg_Pk = np.array([sum_Pk[i]/np.float(n_file) for i in range(len(sum_Pk))])    # average P(k)
        
        # P(k) Comparison
        if correction['name'].lower() == 'true':
            # Should be first and the demoninator otherwise code crashes 
            avg_Pk_true = avg_Pk
    
            # since Tiling Mock does not have Delta P, we scale it from QPM 
            if catalog_name.lower() in ('tilingmock', 'qpm'): 
                qpm_cat_corr = {'catalog': {'name': 'qpm'}, 'correction': {'name': 'true'}, 'spec': {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360}} 
                qpm_Pk = np.array([fc.avgP_interp(avg_k[i], qpm_n_mock, **qpm_cat_corr) for i in range(len(avg_k))])
                qpm_delP = np.array([fc.deltaP_interp(avg_k[i], qpm_n_mock, **qpm_cat_corr) for i in range(len(avg_k))])
                delP = (qpm_delP/qpm_Pk) * avg_Pk
            elif catalog_name.lower() == 'lasdamasgeo':
                ldg_cat_corr = {'catalog': {'name': 'lasdamasgeo'}, 'correction': {'name': 'true'}, 'spec': {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360}} 
                ldg_delP = np.array([fc.deltaP_interp(avg_k[i], ldg_n_mock, **ldg_cat_corr) for i in range(len(avg_k))])
                #print ldg_delP
                delP = ldg_delP
            else: 
                raise NameError("blabadklfjalsdfkj") 
        else: 
            chi2_total = (avg_Pk - avg_Pk_true)**2/delP**2

            chi2 = np.zeros(len(avg_k))
            for i_k, kval in enumerate(avg_k): 
                chi2[i_k] = np.sum(chi2_total[0:i_k])/np.float(len(chi2_total[0:i_k]))
    
            # Plot label 
            if correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'allpeakshot', 'vlospeakshot'):
                # set residual label
                resid_label = ''.join([r"$\mathtt{", correction['fit'], ", ", correction['name'], ",", 
                    ''.join([str(correction['sigma']), str(correction['fpeak'])]), "}$"])
            else:
                resid_label = ''.join([r"$\mathtt{", correction['name'], "}$"])

            sub.scatter(avg_k, chi2, color=pretty_colors[i_corr], label=resid_label)
            print resid_label 
            print chi2[-5:-1]

        try: 
            corr_str 
        except NameError: 
            corr_str = ''.join(['_', correction['name']]) 
        else: 
            if correction['name'] in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'vlospeakshot'): 
                corr_str = ''.join([corr_str, '_', correction['fit'], '_', correction['name'], 
                    '_sigma', str(correction['sigma']), '_fpeak', str(correction['fpeak'])]) 
            else: 
                corr_str = ''.join([corr_str, '_', correction['name']]) 

        del avg_k
        del avg_Pk

    # set axes
    #if catalog_name.lower() == 'tilingmock': ylimit = [0.0,1.0] 
    #elif catalog_name.lower() == 'qpm': ylimit = [0.0, 0.3] 
    #elif catalog_name.lower() == 'lasdamasgeo': ylimit = [0.0, 0.2] 
    ylimit=[10**-4, 10**2]
    ytext = 1.15
    ylabel = r"$\mathtt{\frac{1}{N[k<k_{max}]}\sum_{k< k_{max}} (\overline{P(k)}-\overline{P(k)_{\rm{True}}})^2/(\Delta P)^2}$"

    fig_name = ''.join(['pk_chisquared_', catalog_name.lower(), '_fibcoll_', corr_str, '_comparison.png'])     

    #sub.text(1.5*10**-3, np.mean(ylimit), ''.join([str(n_file), ' ', catalog_name.upper()]))              # number of mocks + Catalog name 
    sub.set_xscale('log')
    sub.set_yscale('log')
    sub.set_xlim([10**-3,10**0])
    sub.set_ylim(ylimit)
    sub.set_xlabel(r'$\mathtt{k_{max}}$', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)
    sub.legend(loc='upper left', scatterpoints=1, prop={'size':14})
    
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]), bbox_inches="tight")
    fig.clear()

# galaxy environment 
def plot_nearest_neighbor(n=3, **cat_corr): 
    ''' Plot the distribution of nth nearest neighbor distances
    '''

    fig = plt.figure(1, figsize=(15, 15)) 
    sub = fig.add_subplot(111) 
    
    # get nth nearest neighbor distance 
    neigh_dist = genv.n_nearest(n = n, **cat_corr) 

    env_hist, env_binedges = np.histogram(neigh_dist, range=[0.0, 50.0], bins=25) 
    env_bin_low = env_binedges[:-1]
    env_bin_high = env_binedges[1:]
    env_bin_mid = [ 0.5*(env_bin_low[i] + env_bin_high[i])
            for i in range(len(env_bin_low)) ]

    sub.plot(env_bin_mid, env_hist, lw=4) 
    fig.savefig('nearest_neighbor_dist.png', bbox_inches='tight')

def plot_catalog_nearest_neighbor(n=3, cat='lasdamasgeo'):
    ''' Plot the combined distribution of nth nearest neighbor distances for a catalog
    '''
    fig = plt.figure(1, figsize=(15, 15)) 
    sub = fig.add_subplot(111) 

    catalog = {} 
    correction = {'name': 'upweight'} 
    if cat == 'lasdamasgeo': 
        catalog['name'] = cat

        for i_mock in range(1, 11): 
            for letter in ['a', 'b', 'c', 'd']: 
                catalog['n_mock'] = i_mock 
                catalog['letter'] = letter 

                i_cat_corr = {
                        'catalog': catalog, 
                        'correction': correction} 

                print catalog['name'], catalog['n_mock'], catalog['letter']
                
                try: 
                    neigh_dist
                except NameError: 
                    neigh_dist = genv.n_nearest(n = n, **i_cat_corr)
                else: 
                    neigh_dist = np.concatenate((neigh_dist, genv.n_nearest(n = n, **i_cat_corr))) 
    else: 
        raise NameError('asdf')  

    env_hist, env_binedges = np.histogram(neigh_dist, range=[0.0, 50.0], bins=25) 
    env_bin_low = env_binedges[:-1]
    env_bin_high = env_binedges[1:]
    env_bin_mid = [ 0.5*(env_bin_low[i] + env_bin_high[i])
            for i in range(len(env_bin_low)) ]

    sub.plot(env_bin_mid, env_hist, lw=4) 
    fig.savefig('combined_nearest_neighbor_dist.png', bbox_inches='tight')

#--------------------------------------------------------------------------------
# Functions used in plotting
def residual(arr1, arr2): 
    if len(arr1) != len(arr2): 
        raise TypeError("Input array lengths do not match.")
    else: 
        resid = np.array([arr1[i]/arr2[i] for i in range(len(arr1))])
    return resid

def classify_triangles(k1, k2, k3, triangle='equilateral'): 
    '''
    Given k1, k2, k3, returns indices for (k1,k2,k3) that satify the specified triangle type  
    '''
    maxk = np.array([np.max([k1[i], k2[i], k3[i]]) for i in range(len(k1))])
    mink = np.array([np.min([k1[i], k2[i], k3[i]]) for i in range(len(k1))])
    if triangle == 'equilateral':       # only keep equilateral triangles
        triangle_index = (k1 == k2) & (k2 == k3) 
    elif triangle == 'acute':           # acute triangle
        triangle_index = (k1**2 + k2**2 > k3**2) & (k2**2 + k3**2 > k1**2) & (k3**2 + k1**2 > k2**2)
    elif triangle == 'obtuse':          # obtuse triangle
        triangle_index = (k1**2 + k2**2 < k3**2) | (k2**2 + k3**2 < k1**2) | (k3**2 + k1**2 < k2**2)
    elif triangle == 'extended':        # extended triangle   
        triangle_index = maxk/mink > 3.0
    return triangle_index

#--------------------------------------------------------------------------------
# Figures for fc_paper

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
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
        elif catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 
        else: 
            raise NameError('not coded yet') 

        # Gaussian fit
        if catalog['name'].lower() == 'lasdamasgeo': 
            corr = {'name': 'peakshot', 'sigma': 6.5, 'fpeak': 0.76, 'fit':'gauss'}
        elif catalog['name'].lower() == 'qpm': 
            corr = {'name': 'peakshot', 'sigma': 4.4, 'fpeak':0.65, 'fit':'gauss'}
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
                    #print Pk_var/avg_Pk_true
                    sub.plot(avg_k, Pk_var/avg_Pk_true+1.0, lw=2, ls='--', color=pretty_colors[0], label='Sample Variance')
                else: 
                    sub.plot(avg_k, Pk_var/avg_Pk_true, lw=4, color=pretty_colors[0], label='Sample Variance')
            else: 
                corr_label = 'Hahn et al.'
                corr_color = 7  

                # set residual label
                #resid_label = ''.join([r"$\mathtt{", "\overline{P(k)_{", corr_label,"}}/\overline{P(k)_{True}}}$"])
                
                # plot residual 
                sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=pretty_colors[corr_color], label=corr_label)
                #print resid_label 
                print catalog['name'].lower() 
                print np.min(residual(avg_Pk, avg_Pk_true)), np.max(residual(avg_Pk, avg_Pk_true))
    

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
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]), bbox_inches="tight")
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

def chi_squared(): 
    # Chi-squared calculation 
    plot_pk_chi2_comp('lasdamasgeo', 40, 
            [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma':6.5, 'fpeak':0.76, 'fit':'gauss'}]) 

    plot_pk_chi2_comp('qpm', 100, 
            [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}]) 

    plot_pk_chi2_comp('tilingmock', 1, 
            [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma':4.8, 'fpeak':0.62, 'fit':'gauss'}])

if __name__=="__main__":
    cat_corr = {'catalog': {'name': 'lasdamasgeo', 'n_mock': 1, 'letter': 'a'}, 
            'correction': {'name': 'upweight'}} 

    plot_catalog_nearest_neighbor(n=3, cat='lasdamasgeo') 
    '''
    # --------------------------------------------------------------------------------------------------------------------------------------------
    #qpm_corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name':'shotnoise'}, {'name': 'peaknbar', 'sigma': 4.4, 'fpeak': 1.0, 'fit': 'gauss'}, {'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}]#, {'name': 'vlospeakshot', 'sigma':580, 'fpeak':0.62, 'fit':'gauss'}]#, {'name': 'vlospeakshot', 'sigma':580, 'fpeak':0.65, 'fit':'gauss'}] 
    qpm_corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name': 'floriansn'}, {'name': 'hectorsn'}, {'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}, {'name': 'peakshot', 'fit': 'true', 'fpeak': 0.65}]
    #, {'name': 'peaknbar', 'sigma': 4.4, 'fpeak': 1.0, 'fit': 'gauss'}, 
    #, {'name': 'vlospeakshot', 'sigma':580, 'fpeak':0.62, 'fit':'gauss'}]#, {'name': 'vlospeakshot', 'sigma':580, 'fpeak':0.65, 'fit':'gauss'}] 
    #plot_pk_fibcol_comp('qpm', 100, qpm_corr_methods, resid='True') 
    #qpm_corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}, {'name': 'floriansn'}, {'name': 'hectorsn'}]
    qpm_corr_methods = [{'name': 'true'}, {'name':'upweight'}]#, {'name':'tailupw'}]#{'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'true'}]#, {'name': 'peakshot', 'fit': 'true', 'fpeak': 0.63}]
    #plot_pk_fibcol_comp('qpm', 49, qpm_corr_methods, resid='True') 
    #plot_p2k_fibcol_comp('qpm', 10, qpm_corr_methods, resid='True') 
    qpm_corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name': 'floriansn'}]
    #plot_pk_fibcol_comp('qpm', 49, qpm_corr_methods, resid='True') 

    #plot_p2k_fibcol_comp('qpm', 100, corr_methods, resid='False') 
    ldg_corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name': 'peakshot', 'sigma': 6.5, 'fpeak': 0.76, 'fit': 'gauss'}, {'name': 'floriansn'}, {'name': 'hectorsn'}]
    #ldg_corr_methods = [{'name': 'true'},  {'name': 'shotnoise'}, {'name': 'peaknbar', 'sigma': 6.5, 'fpeak': 1.0, 'fit': 'gauss'}] #, {'name': 'floriansn'}, {'name': 'hectorsn'}]
    #plot_pk_fibcol_comp('lasdamasgeo', 39, ldg_corr_methods, resid='True') 
    #plot_p2k_fibcol_comp('lasdamasgeo', 1, ldg_corr_methods, resid='True') 
    #plot_pk_fibcol_comp('qpm', 100, corr_methods, resid='True', quad=True) 
    #plot_pk_fibcol_comp('lasdamasgeo', 1, corr_methods, resid='True') 
    #
    tm_corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name': 'peakshot', 'sigma':4.8, 'fpeak':0.63, 'fit':'gauss'}, {'name': 'peakshot', 'sigma':4.8, 'fpeak':1.0, 'fit':'gauss'}, {'name': 'peakshot', 'sigma':4.8, 'fpeak':0.0, 'fit':'gauss'}]#, {'name': 'floriansn'}, {'name': 'hectorsn'}]
    ##{'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma':5.4, 'fpeak':0.57, 'fit':'expon'}, {'name': 'vlospeakshot', 'sigma':620, 'fpeak':0.58, 'fit':'gauss'}] 
    #plot_pk_fibcol_comp('tilingmock',1, corr_methods, resid='True') 
    #plot_p2k_fibcol_comp('tilingmock',1, tm_corr_methods, resid='True') 
    #plot_p2k_fibcol_comp('tilingmock',1, tm_corr_methods, resid='False') 
    
    #plot_z_dist_fcpaper(cat_corrs)
    #plot_pk_fcpaper(['lasdamasgeo', 'qpm', 'tilingmock', 'cmass']) 
    #plot_pk_shotnoiseonly_fcpaper(['lasdamasgeo', 'qpm', 'tilingmock'])
    #corr_methods = [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma':5.4, 'fpeak':0.57, 'fit':'expon'}, {'name': 'vlospeakshot', 'sigma':700, 'fpeak':0.65, 'fit':'expon'}]
    #plot_pk_fibcol_comp('tilingmock', 1, corr_methods, resid='True')
    '''
