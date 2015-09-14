"""

Plot power/bispectrum

Author(s): ChangHoon Hahn

"""
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt
import os.path
from matplotlib.collections import LineCollection

# --- Local --- 
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from spec.spec import Spec

def plot_pk_comp(cat_corrs, n_mock, quad=False, type='ratio', **kwargs): 
    ''' Plot comparison of average power spectrum monopole or quadrupole (avg(P(k))) 
    for multiple a list of catalog and correction specifications. Main use is to 
    compare the effects of fiber collisions correction method. However, it can be
    used to compare any power spectra as long as cat_corr dictionary is specified.

    --------------------------------------------------------------------------
    Paramters
    --------------------------------------------------------------------------
    cat_corrs : list of catalog correction dictionary 
    n_mock : number of mocks 
    quad : If True, then plot quadrupole. If False plot monopole
    type : 'regular' compares the actual P2(k) values, 'ratio' compares the ratio 
    with the true P2(k), 'residual' compares the difference with the true P2(k) 
    
    --------------------------------------------------------------------------
    Notes
    --------------------------------------------------------------------------
    * Long ass code with a lot of idiosyncracies.
    * Make sure k values agree with each other. 
    
    --------------------------------------------------------------------------
    Example
    --------------------------------------------------------------------------
    cat_corrs = [ 
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'true'}},
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}},
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 3.9, 'fpeak': 0.68}} 
            ] 

    plot_pk_comp(cat_corrs, 84, quad=False, type='Pk')
    plot_pk_comp(cat_corrs, 84, quad=False, type='ratio')

    '''

    if 'Ngrid' in kwargs.keys():
        Ngrid = kwargs['Ngrid']
    else: 
        Ngrid = 360 

    if isinstance(n_mock, int): 
        n_mock_list = [ n_mock for i in xrange(len(cat_corrs)) ] 
    else: 
        if len(n_mock) != len(cat_corrs): 
            raise ValueError()
        else: 
            n_mock_list = n_mock

    corr_str = ''

    prettyplot()                         # set up plot 
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(7, 8)) # set up figure 
    sub = fig.add_subplot(111)

    for i_corr, cat_corr in enumerate(cat_corrs):

        catdict = cat_corr['catalog']
        corrdict = cat_corr['correction']
        specdict = {
                'P0': 20000,
                'Lbox': 3600, 
                'Ngrid': Ngrid, 
                'quad': quad 
                }

        n_mock_i = n_mock_list[i_corr]
        
        catcorr_mocks = [{
                    'catalog': {
                        'name': catdict['name'], 
                        'n_mock': i_mock
                        }, 
                    'correction': corrdict, 
                    'spec': specdict
                    }
                for i_mock in xrange(1, n_mock_i + 1)]
                        
        # calculate average P(k) (both monopole or quadrupole) 
        for catcorr_mock in catcorr_mocks: 

            specclass = Spec('pk', catcorr_mock) 
            specclass.read()
            
            k_arr = specclass.k

            if not quad: 
                pk_arr = specclass.p0k
            
            try:
                sum_pk += pk_arr
            except UnboundLocalError: 
                sum_pk = pk_arr

        avg_pk = sum_pk/np.float(n_mock_i) 
        del sum_pk
        
        # Compare P(k) to each other 
        if type == 'Pk':

            sub.plot( 
                    k_arr, avg_pk, 
                    color = pretty_colors[i_corr + 1], 
                    label = plot_label(cat_corr),
                    lw = 4
                    ) 
        
        # Compare k^1.5 * P(k) with each other. Enhances the 
        # BAO signature? (Upon Paco's request). 
        elif type == 'kPk': 

            kPk = k_arr**1.5 * avg_pk
                
            sub.scatter(
                    k_arr, kPk, 
                    color = pretty_colors[i_corr+1], 
                    label = plot_label(cat_corr)
                    )

        # Compare the ratio of the power spectra (P/P_denom)
        elif type == 'ratio':                       

            if i_corr == 0 :        
                avg_pk_denom = avg_pk
                denom_cat = catdict['name']

            else: 
                sub.scatter(
                        k_arr, 
                        avg_pk/avg_pk_denom, 
                        color = pretty_colors[i_corr+1], 
                        label = plot_label(cat_corr)
                        )
                #print avg_Pk[-5:-1]
                #print avg_Pk_true[-5:-1]
                #print resid_label 
        
        del avg_pk
        
        # Specify corrections for figure file name  
        if 'dlospeak' in corrdict['name']: 
            corr_str += ''.join([
                catdict['name'], '_', 
                corrdict['name'], '_', 
                corrdict['fit'], '_',
                '_sigma', str(corrdict['sigma']), 
                'fpeak', str(corrdict['fpeak'])
                ]) 
        else: 
            corr_str += ''.join([ 
                catdict['name'], '_', 
                corrdict['name']
                ]) 
    
    # Dictate the x-range and y-range of the plotting
    # based on type of comparison 
    if type == 'Pk':
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = 10**.5*min(ylimit) 
        else: 
            ylimit = [10**2,10**5.5]
            yytext = 10**2.5
        
        if 'ylabel' in kwargs.keys(): 
            ylabel = kwargs['ylabel']
        else: 
            if quad: 
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
    
    elif type == 'kPk':
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = 10**.5*min(ylimit) 
        else: 
            ylimit = [10**0,10**2.0]
            yytext = 10**0.1

        if 'ylabel' in kwargs.keys(): 
            ylabel = kwargs['ylabel']
        else: 
            ylabel = r'$\mathtt{k^{1.5} P_0(k)}$'

        if 'xscale' in kwargs.keys(): 
            sub.set_xscale(kwargs['xscale']) 
        else: 
            sub.set_xscale('log')

        if 'yscale' in kwargs.keys(): 
            sub.set_yscale(kwargs['yscale'])
        else: 
            sub.set_yscale('log')

        resid_str = '_kPk'

    elif type == 'ratio': 
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = 0.05 + min(ylimit) 
        else: 
            ylimit = [0.5, 1.5] 
            yytext = 0.55

        if quad: 
            ylabel = ''.join([
                r"$\mathtt{\overline{P_2(k)}/\overline{P_2(k)_{\rm{", denom_cat, "}}}}$"
                ])
        else: 
            ylabel = ''.join([
                r"$\mathtt{\overline{P_0(k)}/\overline{P_0(k)_{\rm{", denom_cat, "}}}}$"
                ])
        
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

    sub.set_ylim(ylimit)
    sub.set_xlabel('k (h/Mpc)', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)
    
    # Display the number of mocks for given catalog so that
    # I know how many mocks the P(k) is averaged over.
    n_mock_text = '\n'.join([
                ' '.join([
                    str(n_mock_list[ii]), 
                    ((cat_corrs[ii])['catalog'])['name'].upper()
                    ]) 
                for ii in xrange(len(n_mock_list))
                ])
    sub.text(yxtext, yytext, n_mock_text)

    sub.legend(loc='upper left', scatterpoints=1, prop={'size':14})
    
    try: 
        n_mock_str = '_'.join([str(nm) for nm in n_mock]) 
    except TypeError:
        n_mock_str = str(n_mock) 
    
    if quad: 
        spec_str = 'p2k_'
    else: 
        spec_str = 'p0k_'

    fig_name = ''.join([
        spec_str, 
        n_mock_str, 
        'mock_fibcoll_', 
        corr_str, 
        resid_str, 
        '_comparison_Ngrid', 
        str(Ngrid), 
        '.png'
        ])     

    fig_dir = '/home/users/hahn/powercode/FiberCollisions/figure/'
    fig.savefig(
            ''.join([fig_dir, fig_name]), 
            bbox_inches="tight"
            )
    plt.close()

    return None 

def plot_label(cat_corr): 
    """ Labels for plotting that specify the catalog and 
    correction names. Written for plot_pk_comp but can be used
    elsewhere.
    """
    catdict = cat_corr['catalog'] 
    corrdict = cat_corr['correction']

    if 'dlospeak' in corrdict['name']: 

        label = ' '.join([
            catdict['name'].upper(), 
            ''.join([
                corrdict['name'].upper(), ':', 
                str(corrdict['sigma']), ',', str(corrdict['fpeak'])
                ])
            ])

    elif corrdict['name'] in ('true', 'upweight'): 

        label = ' '.join([ 
            catdict['name'].upper(), 
            corrdict['name'].upper()
            ])
    else: 
        raise NotImplementedError()

    return label 

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

if __name__=='__main__': 
    
    cat_corrs = [ 
            {
                'catalog': {'name': 'nseries'}, 
                'correction': {'name': 'true'}
                },
            {
                'catalog': {'name': 'nseries'}, 
                'correction': {'name': 'upweight'}
                },
            {
                'catalog': {'name': 'nseries'}, 
                'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 3.9, 'fpeak': 0.68}
                } 
            ] 

    plot_pk_comp(cat_corrs, 84, quad=False, type='Pk')
    plot_pk_comp(cat_corrs, 84, quad=False, type='ratio')


"""
        if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'):             # LasDamasGeo 
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

        elif catalog['name'].lower() in ('tilingmock', 'cmass'):       # Tiling Mocks
            i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}

            power = fc_spec.Spec('power', **i_cat_corr)         # read tiling mock P(k)
            power.readfile()
            print power.file_name

            avg_k = power.k

            if quad == True: 
                sum_Pk = power.P2k
            else: 
                sum_Pk = power.Pk

            n_file = 1          # only one file 
        
        elif 'bigmd' in catalog['name'].lower():                        # BigMD
            i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}

            power = fc_spec.Spec('power', **i_cat_corr)         # read tiling mock P(k)
            power.readfile()
            print power.file_name

            avg_k = power.k

            if quad == True: 
                sum_Pk = power.P2k
            else: 
                sum_Pk = power.Pk

            n_file = 1          # only one file 

        elif catalog['name'].lower() in ('qpm', 'nseries', 'patchy'):           # QPM/PATCHY
            n_file = 0  
            for i_mock in range(1, n_mock_i+1): 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock
                i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                    
                power_i = fc_spec.Spec('power', **i_cat_corr)
                #print power_i.file_name
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
"""
    
