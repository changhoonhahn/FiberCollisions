'''

Plotting codes for FiberCollisions project


Author(s): ChangHoon Hahn 


'''


import numpy as np 
import scipy as sp
import os.path
import subprocess
import cosmolopy as cosmos
from matplotlib.collections import LineCollection

# --- Local --- 
import fibcol_data as fc_data
import fibcol_nbar as fc_nbar
import fibcol_spec as fc_spec
import fibcol_utility as fc_util
import fibercollisions as fc
import galaxy_environment as genv

# Plot P(k) monopole or quadrupole ---------------------
def plot_pk_fibcol_comp(cat_corrs, n_mock, quad=False, type='ratio', **kwargs): 
    ''' Plots the comparison of avg(P(k)) (monopole or quadrupole) for multiple 
    fibercollision correciton methods
    
    Paramters
    ---------
    cat_corrs : list of catalog correction dictionary 
    n_mock : number of mocks 
    quad : If True, then plot quadrupole. If False plot monopole
    type : 'regular' compares the actual P2(k) values, 'ratio' compares the ratio 
    with the true P2(k), 'residual' compares the difference with the true P2(k) 

    Notes
    -----
    * Make sure 'true' is first in the correction method list!

    '''
    prettyplot()                         # set up plot 
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(7, 8)) # set up figure 
    sub = fig.add_subplot(111)
    if 'Ngrid' in kwargs.keys():
        Ngrid = kwargs['Ngrid']
    else: 
        Ngrid = 360         # default Ngrid


    for i_corr, cat_corr in enumerate(cat_corrs):     # loop through correction methods
    
        catalog = cat_corr['catalog']
        correction = cat_corr['correction']

        # power/bispectrum properties 
        if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz', 'qpm', 'nseries', 'bigmd', 
                'patchy'): 

            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 
                    'grid': Ngrid, 'quad': quad} 

        elif catalog['name'].lower() == 'tilingmock': 

            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 
                    'grid': Ngrid, 'quad': quad} 
        else: 
            raise NotImplementedError('not coded yet') 

        try: 
            if len(n_mock) == len(cat_corrs): 
                n_mock_i = n_mock[i_corr] 
            else: 
                raise NameError('they need to match') 

        except TypeError: 
            n_mock_i = n_mock 

        if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'):           # LasDamasGeo 

            # compute average[P(k)] for each correction method
            n_file = 0
            for i_mock in range(1, n_mock_i+1):
                for letter in ['a', 'b', 'c', 'd']: 
                    # set catalog correction dictionary for specific file 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_catalog['letter'] = letter
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                    
                    # import power spectrum 
                    power_i = fc_spec.Spec('power', **i_cat_corr)
                    print power_i.file_name
                    power_i.readfile()

                    if quad == True: 
                        PK = power_i.P2k
                    else: 
                        PK = power_i.Pk
                    
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_Pk = PK 
                    else: 
                        sum_Pk += PK 

                    n_file = n_file+1

        elif catalog['name'].lower() in ('tilingmock', 'bigmd'):       # Tiling Mocks

            i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}

            power = fc_spec.Spec('power', **i_cat_corr)         # read tiling mock P(k)
            power.readfile()

            avg_k = power.k

            if quad == True: 
                sum_Pk = power.P2k
            else: 
                sum_Pk = power.Pk

            n_file = 1          # only one file 

        elif catalog['name'].lower() in ('qpm', 'nseries', 'patchy'): 

            n_file = 0  
            for i_mock in range(1, n_mock_i+1): 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock
                i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                    
                power_i = fc_spec.Spec('power', **i_cat_corr)
                power_i.readfile()

                if quad: 
                    PK = power_i.P2k
                else: 
                    PK = power_i.Pk
                    
                try: 
                    avg_k 
                except NameError: 
                    avg_k = power_i.k
                    sum_Pk = PK
                else: 
                    sum_Pk += PK 

                n_file += 1

        else: 
            raise NameError('not yet coded!')

        # calculate the average P(k)
        avg_Pk = [ sum_Pk[i]/np.float(n_file) for i in range(len(sum_Pk)) ] 
        
        if type == 'regular':                       # P(k) Comparison
            
            if i_corr > 0: 
                lstyle = '--' 
            else: 
                lstyle = '-' 

            if correction['name'].lower() == 'true':    # P_true(k) 
            
                lbl = ' '.join([catalog['name'].upper(), correction['name'].upper()])
                sub.plot(avg_k, avg_Pk, 
                        color=pretty_colors[i_corr+1], ls=lstyle,  
                        label=lbl, lw=4)

            elif correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 
                    'peakshot', 'allpeakshot', 'vlospeakshot'):

                if correction['name'] == 'peakshot': 
                    corr_name = 'Hahn' 
                else: 
                    corr_name = correction['name']  
                
                if correction['fit'].lower() in ('expon', 'gauss'): 
                    resid_label = ''.join([
                        corr_name, ': ', 
                        str(correction['sigma']), ',', str(correction['fpeak'])
                        ]) 

                elif correction['fit'].lower() in ('true'): 
                    resid_label = ''.join([
                        corr_name, ': ', str(correction['fpeak'])
                        ]) 
            
                sub.scatter(avg_k, avg_Pk, 
                        color=pretty_colors[i_corr+1], label=resid_label)

            else: 
                lbl = ' '.join([catalog['name'].upper(), correction['name'].upper()])
                sub.scatter(avg_k, avg_Pk, 
                        color=pretty_colors[i_corr+1], label=lbl)

        elif type == 'ratio':                       # P_corr(k)/P_true comparison 
            if i_corr == 0 :        
                avg_Pk_denom = avg_Pk        # P_true(k) 

            else: 
                if correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 
                        'peakshot', 'allpeakshot', 'vlospeakshot'):

                    if correction['name'] == 'peakshot': 
                        corr_name = 'Hahn'
                    else: 
                        corr_name = correction['name']

                    if correction['fit'].lower() in ('expon', 'gauss'): # labels
                        # exponential and gaussian fits 
                        resid_label = ''.join([
                            corr_name, ': ', correction['fit'], ', ' 
                            ','.join([str(correction['sigma']), str(correction['fpeak'])])
                            ])
                        #resid_label = ''.join([r"$\mathtt{", 
                        #    "\overline{P_2(k)_{", correction['fit'], ", ", 
                        #    correction['name'], ",", ''.join([str(correction['sigma']), 
                        #        str(correction['fpeak'])]), "}}/\overline{P_2(k)_{True}}}$"])

                    elif correction['fit'].lower() in ('true'): 
                        # true distribution fits
                        resid_label = ''.join([
                            corr_name, ': ', correction['fit'], ', ', str(correction['fpeak'])
                            ])
                        #resid_label = ''.join([r"$\mathtt{", 
                        #    "\overline{P_2(k)_{", correction['fit'], ", ", 
                        #    correction['name'], ",", str(correction['fpeak']), 
                        #    "}}/\overline{P_2(k)_{True}}}$"])
                    else:
                        raise NameError('asdflkj')
                    
                    # plot ratio
                    sub.scatter(avg_k, residual(avg_Pk, avg_Pk_denom), \
                            color=pretty_colors[i_corr+1], label=resid_label)
                    
                    '''
                    print avg_Pk[-5:-1]
                    print avg_Pk_true[-5:-1]
                    print resid_label 
                    print residual(avg_Pk, avg_Pk_true)[(avg_k > 0.15) & (avg_k < 0.2)]
                    print avg_k[(avg_k > 0.15) & (avg_k < 0.2)]
                    '''
                else:
                    resid_label = ''.join([catalog['name'], ' ', correction['name']])

                    sub.scatter(avg_k, residual(avg_Pk, avg_Pk_denom), \
                            color=pretty_colors[i_corr+1], label=resid_label)
                    '''
                    print resid_label 
                    print residual(avg_P2k, avg_P2k_true)[(avg_k > 0.15) & (avg_k < 0.2)]
                    print avg_k[(avg_k > 0.15) & (avg_k < 0.2)]
                    '''
        
        elif type == 'residual':        # |P_corr(k) - P_true(k)|/Delta P(k) 

            if correction['name'].lower() == 'true': # P_true(k) 
                true_cat_corr = {'catalog': catalog, 
                        'correction': {'name': 'true'}, 'spec': spec} 

                avg_Pk_true = avg_Pk
                delta_P = fc.deltaP(n_mock, **true_cat_corr) 
                delP = delta_P[1]
            else: 
                if correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 
                        'peakshot', 'allpeakshot', 'vlospeakshot'):

                    if correction['name'] == 'peakshot': 
                        corr_name = 'Hahn'
                    else: 
                        corr_name = correction['name']  

                    if correction['fit'].lower() in ('expon', 'gauss'): # residual labels
                        # exponential and gaussian fits 
                        resid_label = ''.join([
                            corr_name, ': ', correction['fit'], ', ' 
                            ','.join([str(correction['sigma']), str(correction['fpeak'])])
                            ])

                    elif correction['fit'].lower() in ('true'): 
                        # true distribution fits
                        resid_label = ''.join([
                            corr_name, ': ', correction['fit'], ',', str(correction['fpeak'])
                            ]) 
                    else:
                        raise NameError('asdflkj')
                    
                    resids = [ np.abs(avg_Pk[i] - avg_Pk_true[i])/delP[i]
                            for i in range(len(avg_Pk)) ] 

                    # plot residual 
                    sub.scatter(avg_k, resids, \
                            color=pretty_colors[i_corr+1], label=resid_label)

                else:
                    resid_label = correction['name'] 

                    resids = [ np.abs(avg_Pk[i] - avg_Pk_true[i])/delP[i]
                            for i in range(len(avg_Pk)) ] 
                    
                    sub.scatter(avg_k, resids, \
                            color=pretty_colors[i_corr+1], label=resid_label)

        else: 
            raise NotImplementedError('asdfasdfasdfadf') 

        try: 
            if correction['name'].lower() in ('peak', 'peaknbar', 'peakshot'): 
                corr_str += ''.join(['_', catalog['name'].lower(), '_', 
                    correction['fit'], '_', correction['name'], 
                    '_sigma', str(correction['sigma']), '_fpeak', str(correction['fpeak'])]) 
            else: 
                corr_str += ''.join(['_', catalog['name'].lower(), '_', correction['name']]) 

        except NameError: 
            corr_str = ''.join(['_', catalog['name'].lower(), '_', correction['name']]) 

        del avg_k
        del avg_Pk

    # set axes
    if type == 'regular':
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = 10**.5*min(ylimit) 
        else: 
            ylimit = [10**2,10**5.5]
            yytext = 10**2.5
        
        if 'ylabel' in kwargs.keys(): 
            ylabel = kwargs['ylabel']
        else: 
            if quad == True: 
                ylabel = r'$\mathtt{P_2(k)}$'
            else: 
                ylabel = r'$\mathtt{P_0(k)}$'

        if 'xscale' in kwargs.keys(): 
            sub.set_xscale(kwargs['xscale']) 
        else: 
            sub.set_xscale('log')

        if 'yscale' in kwargs.keys(): 
            sub.set_yscale(kwargs['yscale'])
        else: 
            sub.set_yscale('log')

        resid_str = ''

    elif type == 'ratio': 
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = 0.05 + min(ylimit) 
        else: 
            ylimit = [0.5, 2.0] 
            yytext = 0.55

        if quad == True: 
            ylabel = r"$\mathtt{\overline{P_2(k)}/\overline{P_2(k)_{\rm{True}}}}$"
        else: 
            ylabel = r"$\mathtt{\overline{P_0(k)}/\overline{P_0(k)_{\rm{True}}}}$"
        
        if 'xscale' in kwargs.keys(): 
            sub.set_xscale(kwargs['xscale']) 
        else: 
            sub.set_xscale('log')

        if 'yscale' in kwargs.keys(): 
            sub.set_yscale(kwargs['yscale'])
        else: 
            pass

        resid_str = '_ratio'

    elif type == 'residual':
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = np.mean(ylimit) 
        else: 
            ylimit = [0.0, 5.0]
            yytext = 2.5

        if quad == True: 
            ylabel = r"$\mathtt{|\overline{P_2(k)} - \overline{P_2(k)_{\rm{True}}}|/\Delta P_2}$"
        else: 
            ylabel = r"$\mathtt{|\overline{P_0(k)} - \overline{P_0(k)_{\rm{True}}}|/\Delta P_0}$"

        if 'xscale' in kwargs.keys(): 
            sub.set_xscale(kwargs['xscale']) 
        else: 
            sub.set_xscale('log')

        if 'yscale' in kwargs.keys(): 
            sub.set_yscale(kwargs['yscale'])
            if (kwargs['yscale'] == 'log') and \
                    ('yrange' not in kwargs.keys()):
                ylimit = [10**-3, 10**1] 
        else: 
            pass

        resid_str = '_residual'
    else: 
        raise NotImplementedError('asdfasdfasdf') 

    if 'xrange' in kwargs.keys():   # specify x-range 
        sub.set_xlim(kwargs['xrange']) 
        yxtext = 1.5*min(kwargs['xrange'])
    else: 
        sub.set_xlim([10**-3,10**0])
        yxtext = 1.5*10**-3

    try: 
        sub.text(yxtext, yytext, 
                '\n'.join([
                    ' '.join([str(n_mock[ii]), ((cat_corrs[ii])['catalog'])['name'].upper()]) 
                    for ii in range(len(n_mock))]))
    except TypeError: 
        sub.text(yxtext, yytext, 
                ''.join([str(n_file), ' ', catalog['name'].upper()]))  # number of mocks + Catalog name 

    sub.set_ylim(ylimit)
    sub.set_xlabel('k', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)

    #if resid == 'False': 
    #    sub.legend(loc='lower left', scatterpoints=1, prop={'size':14})
    #elif resid == 'True': 
    sub.legend(loc='upper left', scatterpoints=1, prop={'size':14})
    
    try: 
        n_mock_str = '_'.join([str(nm) for nm in n_mock]) 
    except TypeError:
        n_mock_str = str(n_mock) 

    # fig file name 
    if quad == True: 
        fig_name = ''.join(['p2k_', n_mock_str, 'mock',
            '_fibcoll_', corr_str, resid_str, '_comparison_Ngrid', str(Ngrid) , '.png'])     
    else: 
        fig_name = ''.join(['p0k_', n_mock_str, 'mock',
            '_fibcoll_', corr_str, resid_str, '_comparison_Ngrid', str(Ngrid), '.png'])     
    fig.savefig(''.join(['figure/', fig_name]), bbox_inches="tight")

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

# Average P(k) (monopole and quadrupole) ------------------------------
def plot_avg_pk_fibcol(catalog_name, n_mock, corr_method, quad=False, **kwargs):
    ''' Plot average P(k) monopole or quadrupole of n_mock realizations

    Paramters
    ---------
    catalog_name : Name of mock catalog 
    n_mock : number of mocks 
    corr_method : correction method

    '''
    prettyplot()                         # set up plot 
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(7, 8)) # set up figure 
    sub = fig.add_subplot(111)

    catalog = {'name': catalog_name}    # catalog
    correction = corr_method            # correction 

    # get spectrum info for cat_corr dictionary 
    if catalog_name.lower() in ('lasdamasgeo', 'ldgdownnz', 'qpm'): 
        spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 
                'grid':360, 'quad': quad} 
    elif catalog_name.lower() == 'tilingmock': 
        spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360, 
                'quad': quad} 
    else: 
        raise NotImplementedError('not coded yet') 

    if catalog_name.lower() in ('lasdamasgeo', 'ldgdownnz'):           # LasDamasGeo 

        # plot P(k) for each mock realization
        for i_mock in range(1, n_mock+1):                       
            for letter in ['a', 'b', 'c', 'd']: 

                # set catalog correction dictionary for specific file 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock
                i_catalog['letter'] = letter
                i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                
                # import power spectrum 
                power_i = fc_spec.Spec('power', **i_cat_corr)
                power_i.readfile()
                
                if quad == True: 
                    power_Pk = power_i.P2k
                else: 
                    power_Pk = power_i.Pk 

                sub.plot(power_i.k, power_Pk, 
                        color=pretty_colors[i_mock % 20], lw=2) 

    elif catalog_name.lower() == 'tilingmock':          # Tiling Mocks

        i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}

        power = fc_spec.Spec('power', **i_cat_corr)         # read tiling mock P(k)
        power.readfile()

        if quad == True: 
            power_Pk = power_i.P2k
        else: 
            power_Pk = power_i.Pk 

        sub.plot(power_i.k, power_Pk, 
                color=pretty_colors[0], lw=2) 

    elif catalog_name.lower() == 'qpm': 
        n_file = 0  
        for i_mock in range(1, n_mock+1): 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock
            i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                
            power_i = fc_spec.Spec('power', **i_cat_corr)
            power_i.readfile()
                
            if quad == True: 
                power_Pk = power_i.P2k
            else: 
                power_Pk = power_i.Pk 

            sub.plot(power_i.k, power_Pk, 
                    color=pretty_colors[i_mock % 20], lw=2) 

            n_file += 1

    else: 
        raise NameError('not yet coded!')

    # read in average P(k)  
    avg_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
    avg_k, avg_Pk = fc.avgPk(n_mock, clobber=True, **avg_cat_corr)
    
    avg_k, delPk = fc.deltaP(n_mock, clobber=True, **avg_cat_corr) 

    sub.plot(avg_k, avg_Pk, color='black', lw=6, label=r'$\mathtt{Average\;P(k)}$')
    #sub.fill_between(avg_k, avg_Pk - delPk, avg_Pk + delPk) 
    #sub.errorbar(avg_k, avg_Pk, yerr=delPk, fmt=None, capsize=4, color='black', 
    #        elinewidth=4)#, label=r'$\mathtt{Average\;P(k)}$')


    # set axes
    ylimit = [10**2,10**5.5]
    yytext = 10**5.2
    if quad == True: 
        quad_str = '2'
        ylabel = r'$\mathtt{P_2(k)}$'
    else: 
        quad_str = '0'
        ylabel = r'$\mathtt{P_0(k)}$'

    sub.set_yscale('log')
    resid_str = ''

    if 'xrange' in kwargs.keys():   # specify x-range 
        sub.set_xlim(kwargs['xrange']) 
        yxtext = 1.5*min(kwargs['xrange'])
    else: 
        sub.set_xlim([10**-3,10**0])
        yxtext = 1.5*10**-3

    sub.text(yxtext, yytext, 
            ''.join([str(n_mock), ' ', catalog_name.upper()]))  # number of mocks + Catalog name 
    sub.set_xscale('log')

    sub.set_ylim(ylimit)
    sub.set_xlabel('k', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)

    sub.legend(loc='upper left', scatterpoints=1, prop={'size':14})
    
    fig_name = ''.join(['avg_', 'p', quad_str, 'k_', catalog_name.lower(), 
        '_fibcoll_', correction['name'], '.png'])     
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', 
                fig_name]), 
            bbox_inches="tight")
    fig.clear()

# others 
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

def plot_comdis2z_test(): 
    ''' Plot to test interpolation of comdis2z

    Notes
    -----
    * linear, quadratic, cubic: cubic fits best 
    '''
    prettyplot() 
    pretty_colors = prettycolors()  
    
    # comsology 
    omega_m = 0.31
    cosmo = {} 
    cosmo['omega_M_0'] = omega_m 
    cosmo['omega_lambda_0'] = 1.0 - omega_m 
    cosmo['h'] = 0.676
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 
    
    z_arr = np.array(np.arange(0.0, 1.05, 0.05))
    dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo)*cosmo['h']
    
    z_fine_arr = np.array(np.arange(0.43, 0.7, 0.01))
    dm_fine_arr = cosmos.distance.comoving_distance(z_fine_arr, **cosmo)*cosmo['h']

    fig = plt.figure()
    sub = fig.add_subplot(111) 

    for i_interp, interp in enumerate(['linear', 'quadratic', 'cubic']): 

        dmz_interp = sp.interpolate.interp1d(dm_arr, z_arr, kind=interp) 

        z_dm_interp = np.array([dmz_interp(dm_fine) for dm_fine in dm_fine_arr])

        sub.plot(dm_fine_arr, z_fine_arr - z_dm_interp, 
                lw=4, c=pretty_colors[i_interp], label=interp)
    
    sub.set_xlabel(r'$D_c$') 
    sub.set_ylabel(r'$\Delta z$') 
    sub.legend()

    fig.savefig('figure/comdis2z_test.png')
    fig.clear()
    plt.close(fig)

# nbar(z) ------------------------------------------------------------
def plot_nbar_comparison(cat_corr_list, type='ratio', **kwargs):
    ''' Plot n(z) values for a list of correction methods 
    
    
    Paramters
    ---------
    cat_corr_list : list of catalog, correction dictionary  
    type : Comparison type ('regular', 'ratio', etc) 

    Notes
    -----
    * If type == 'ratio', put denominator first (the true) 

    '''
    # make pretty 
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(1) 
    sub = fig.add_subplot(111)

    nbar_list, nbar_label = [], [] 
    cat_corr_str = ''
    for i_cc, cat_corr in enumerate(cat_corr_list): 
        
        catalog = cat_corr['catalog']
        correction = cat_corr['correction']
        cat_corr_label = ''.join([catalog['name'].lower(), ': ', correction['name']]) 
        cat_corr_str += ''.join(['_', catalog['name'].lower(), '.', correction['name']]) 
        
        nb = fc_nbar.nbar(**cat_corr)   # import nbar class 

        if type == 'regular': 
            # plot nbar(z) 
            sub.plot(nb.zmid, nb.nbar, c=pretty_colors[i_cc], lw=4, label=cat_corr_label) 
    
        elif type == 'ratio': 
            if i_cc == 0: 
                denom_nbar = nb.nbar
            else: 
                nbar_list.append(nb.nbar) 
                nbar_label.append(cat_corr_label) 
                nbar_z = nb.zmid
    
    # plot nbar(z) ratio 
    if type == 'ratio': 
        for i_nbz, nbz in enumerate(nbar_list): 
            sub.plot(nbar_z, nbz/denom_nbar, lw=4, c=pretty_colors[i_nbz], label=nbar_label[i_nbz])

    sub.legend() 

    if 'xrange' in kwargs.keys(): 
        sub.set_xlim(kwargs['xrange'])
    else: 
        sub.set_xlim([0.43, 0.7])

    if 'yrange' in kwargs.keys(): 
        sub.set_ylim(kwargs['yrange'])
    else: 
        if type =='ratio': 
            sub.set_ylim([0.8, 1.2])
    
    sub.set_xlabel(r'$z$ (Redshift)') 
    if type == 'ratio': 
        sub.set_ylabel(r'$\mathtt{\bar{n}(z)/\bar{n}_{true}(z)}$')
        
    fig_file = ''.join(['figure/', 
        'nbar', cat_corr_str, '_', type, '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 

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

def chi_squared(): 
    # Chi-squared calculation 
    plot_pk_chi2_comp('lasdamasgeo', 40, 
            [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma':6.5, 'fpeak':0.76, 'fit':'gauss'}]) 

    plot_pk_chi2_comp('qpm', 100, 
            [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}]) 

    plot_pk_chi2_comp('tilingmock', 1, 
            [{'name': 'true'}, {'name':'upweight'}, {'name':'floriansn'}, {'name':'hectorsn'}, {'name': 'peakshot', 'sigma':4.8, 'fpeak':0.62, 'fit':'gauss'}])

if __name__=="__main__":

    #cat_corr_list = [] 
    #catalog = {'name': 'nseries'}
    #for corr in ['true', 'upweight', 'scratch_peakknown']: 
    #    cat_corr_list.append({'catalog': catalog, 'correction': {'name': corr}}) 

    #plot_nbar_comparison(cat_corr_list, type='ratio', xrange=[0.16, 0.4], yrange=[0.95, 1.05])

    #cat_corr = {'catalog': {'name': 'lasdamasgeo', 'n_mock': 1, 'letter': 'a'}, 
    #        'correction': {'name': 'upweight'}} 

    #plot_catalog_nearest_neighbor(n=3, cat='lasdamasgeo') 
    #plot_avg_pk_fibcol('lasdamasgeo', 40, {'name': 'true'}, quad=True)   
    #plot_avg_pk_fibcol('qpm', 100, {'name': 'true'}, quad=True)   

    #plot_avg_pk_fibcol('lasdamasgeo', 40, {'name': 'true'}, quad=False)   
    #plot_avg_pk_fibcol('qpm', 100, {'name': 'true'}, quad=False)   
    #plot_comdis2z_test()
    catcorr_methods = [
            {'catalog': {'name': 'bigmd'}, 'correction': {'name': 'upweight'}},
            {'catalog': {'name': 'patchy'}, 'correction': {'name': 'upweight'}}, 
            {'catalog': {'name': 'qpm'}, 'correction': {'name': 'upweight'}} 
            ]

    plot_pk_fibcol_comp(catcorr_methods, [10,10,1],
            quad=False, Ngrid=360, type='regular', 
            xrange=[0.001, 1.0], yrange=[10**3, 3*10**5])
    plot_pk_fibcol_comp(catcorr_methods, [10,10,1],
            quad=False, Ngrid=360, type='ratio', 
            xrange=[0.001, 1.0], yrange=[0.8, 1.2])

    ''' 
    catcorr_methods = [
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'true'}}, 
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}}, 
            {'catalog': {'name': 'bigmd'}, 'correction': {'name': 'true'}}, 
            {'catalog': {'name': 'bigmd'}, 'correction': {'name': 'upweight'}}
            ]

    plot_pk_fibcol_comp(catcorr_methods, [10, 10, 1, 1],
            quad=False, Ngrid=360, type='regular', 
            xrange=[0.001, 1.0], yrange=[10**3, 3*10**5])
    plot_pk_fibcol_comp(catcorr_methods, [10, 10, 1, 1],
            quad=False, Ngrid=360, type='ratio', 
            xrange=[0.001, 1.0], yrange=[0.8, 1.2])
    
    catcorr_methods = [
            {'catalog': {'name': 'bigmd'}, 'correction': {'name': 'true'}}, 
            {'catalog': {'name': 'bigmd'}, 'correction': {'name': 'upweight'}}
            ]

    plot_pk_fibcol_comp(catcorr_methods, 1,
            quad=False, Ngrid=360, type='regular', 
            xrange=[0.001, 1.0], yrange=[10**3, 3*10**5])
    plot_pk_fibcol_comp(catcorr_methods, 1,
            quad=False, Ngrid=360, type='ratio', 
            xrange=[0.001, 1.0], yrange=[0.8, 1.2])
    '''
    """
    nseries_corr_methods = [{'name': 'true'}, {'name':'upweight'}] #,
            #{'name': 'peakshot', 'sigma':3.8, 'fpeak': 0.7, 'fit': 'gauss'}, 
            #{'name': 'scratch_peakknown'},  
            #{'name': 'scratch_peakknown_gauss'},  
            #{'name': 'scratch_peakknown_gauss_divide'} 
            #]
            #{'name': 'peakshot', 'fpeak': 0.7, 'fit': 'true'} 
            #]

    plot_pk_fibcol_comp('nseries', 84, nseries_corr_methods, 
            quad=True, Ngrid=360, type='regular', 
            xrange=[0.001, 1.0], yrange=[10**3, 3*10**5])
    plot_pk_fibcol_comp('nseries', 84, nseries_corr_methods, 
            quad=True, Ngrid=360, type='ratio', 
            xrange=[0.001, 1.0], yrange=[0.2, 2.0]) 
        
    plot_pk_fibcol_comp('nseries', 84, nseries_corr_methods, 
            quad=True, Ngrid=360, type='residual', yscale='log', 
            xrange=[0.001, 1.0], yrange=[0.02, 100.0]) 
    
    ldg_corr_methods = [{'name': 'true'}, {'name':'upweight'},
            {'name': 'peakshot', 'sigma':6.5, 'fpeak': 0.76, 'fit': 'gauss'}, 
            {'name': 'scratch_peakknown'}]
            
            #, {'name': 'scratch_peakknown_ang'}, {'name': 'scratch_peakknown_gauss'}]
            #{'name': 'peakshot', 'fpeak': 0.7, 'fit': 'true'} 
            #]

    plot_pk_fibcol_comp('lasdamasgeo', 10, ldg_corr_methods, 
            quad=True, Ngrid=360, type='regular', 
            xrange=[0.001, 1.0], yrange=[10**3, 3*10**5])
    plot_pk_fibcol_comp('lasdamasgeo', 10, ldg_corr_methods, 
            quad=True, Ngrid=360, type='ratio', 
            xrange=[0.001, 1.0], yrange=[0.2, 2.0]) 
    """
