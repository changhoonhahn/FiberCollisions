'''

Plotting codes for FiberCollisions project


Author(s): ChangHoon Hahn 


'''
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt
import os.path
from matplotlib.collections import LineCollection

# --- Local --- 
from utility.plotting import prettyplot
from utility.plotting import prettycolors 

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
    cat_corrs = [{'catalog': {'name': 'nseries', 'n_mock': 1}, 
        'correction': {'name': 'photozenvpeakshot', 'fit': 'gauss', 'sigma': 4.0, 'fpeak': 0.69, 'n_NN':5}}]
    plot_peakcorrection_dlos_check(cat_corrs)
   
    '''
    catcorr_methods = [
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'true'}}, 
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}}, 
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'peakshot', 
                'sigma': '4.0', 'fpeak': '0.69', 'fit':'gauss'}}, 
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'scratch_peakknown'}}, 
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'photozpeakshot', 
                'sigma': '4.0', 'fpeak': '0.69', 'fit':'gauss'}} 
            ]

    catcorr_methods = [
            {'catalog': {'name': 'patchy'}, 'correction': {'name': 'true'}}, 
            {'catalog': {'name': 'patchy'}, 'correction': {'name': 'upweight'}}
            ]
    n_mock_list = 100
    plot_pk_fibcol_comp(catcorr_methods, n_mock_list, \
            quad=False, Ngrid=960, type='regular', 
            xrange=[0.01, 1.0], yrange=[10**2, 3*10**5])
    plot_pk_fibcol_comp(catcorr_methods, n_mock_list, 
            quad=False, Ngrid=960, type='ratio', 
            xrange=[0.01, 1.0], yrange=[0.0, 2.0])
    #plot_pk_fibcol_comp(catcorr_methods, n_mock_list, 
    #        quad=True, Ngrid=960, type='kPk', 
    #        xrange=[0.01, 1.0], yrange=[10**0, 3*10**3])
    ''' 
    #catcorr_methods = [
    #        {'catalog': {'name': 'cmass', 'cosmology': 'fiducial'}, 
    #            'correction': {'name': 'peakshot', 'sigma':6.9, 'fpeak':0.7, 'fit':'gauss'}},
    #        {'catalog': {'name': 'bigmd3'}, 'correction': {'name': 'true'}}
    #        ]

    #plot_pk_fibcol_comp(catcorr_methods, n_mock_list, 
    #        quad=False, Ngrid=960, type='kPk', 
    #        xrange=[0.001, 1.0], yrange=[10**0, 3*10**3])
