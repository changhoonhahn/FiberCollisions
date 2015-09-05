'''

Make dLOS distribution peak corrected CMASS and LOWZ 
combined galaxy catalogs 

'''
from spec.data import Data
from dlos.dlos import Dlos
from dlos import fitting as dlos_fit
from plotting.dlos_plot import Plotdlos
from dlos.fitting import peak_gauss 
from plotting.test_dlos_plots import dlospeak_dlos_test
    
def build_combinedsample(): 
    """ Build CMASS LOWZ combined sample for three separate sector regions and two 
    redshfit bins
    """
    for zbin_str in ['_low', '_high']:
        for specifier_str in ['', 'e2', 'e3']: 
            cat_corr = {
                    'catalog': {'name': 'cmasslowz'+specifier_str+zbin_str}, 
                    'correction': {'name': 'upweight'} 
                    }
            corrdata = Data('data', cat_corr)
            print 'Constructing ', corrdata.file_name
            corrdata.build()

    return None 

def build_combinedsample_random(): 
    """ Build CMASS LOWZ combined sample for three separate sector regions and two 
    redshfit bins
    """
    for zbin_str in ['_low', '_high']:
        for specifier_str in ['', 'e2', 'e3']: 
            cat_corr = {
                    'catalog': {'name': 'cmasslowz'+specifier_str+zbin_str}, 
                    'correction': {'name': 'upweight'} 
                    }
            corrdata = Data('random', cat_corr)
            print 'Constructing ', corrdata.file_name
            corrdata.build()

    return None 

def build_dlospeak_corr_combinedsample(): 
    """ Build "dlospeak" corrected CMASS LOWZ combined sample for three separate sector regions and two 
    redshfit bins
    """
    for zbin_str in ['_low', '_high']:
        for specifier_str in ['', 'e2', 'e3']: 

            if 'low' in zbin_str: 
                sigma_z = 6.9
                fpeak_z = 0.72
            elif 'high' in zbin_str: 
                sigma_z = 6.3
                fpeak_z = 0.7

            cat_corr = {
                    'catalog': {
                        'name': 'cmasslowz'+specifier_str+zbin_str
                        }, 
                    'correction': {
                        'name': 'dlospeak', 
                        'fit': 'gauss', 
                        'sigma': sigma_z, 
                        'fpeak': fpeak_z
                        } 
                    }
            corrdata = Data('data', cat_corr)
            print 'Constructing ', corrdata.file_name
            corrdata.build()

def build_dlos(): 
    """ Calculate line-of-sight displacement for each of the CMASS LOWZ 
    combined sample 
    """

    for zbin_str in ['_low', '_high']:
        for specifier_str in ['', 'e2', 'e3']: 
            cat_corr = {
                    'catalog': {'name': 'cmasslowz'+specifier_str+zbin_str}, 
                    'correction': {'name': 'upweight'} 
                    }
            dlosclass = Dlos(cat_corr, clobber=True)
            print dlosclass.file_name 
            dlosclass.read()

def dlos_peak_fit_combined_dlos(zbin, **kwargs): 
    """ Fit peak of the combined dLOS distribution of all three 
    CMASS LOWZ combined sample sectors for the two redshift bins
    separately
    """

    zbin_str = '_'+zbin 
    
    if 'combined_dlos' in kwargs.keys(): 
        combined_dlos = kwargs['combined_dlos'] 

    for specifier_str in ['', 'e2', 'e3']: 
        cat_corr = {
                'catalog': {'name': 'cmasslowz'+specifier_str+zbin_str}, 
                'correction': {'name': 'upweight'} 
                }
        dlosclass = Dlos(cat_corr)
        dlosclass.read()

        dlos_i = dlosclass.dlos 

        try: 
            combined_dlos += list(dlos_i)
        except NameError: 
            combined_dlos = list(dlos_i)
    
    bestfit_fpeak, bestfit_sigma, bestfit_amp = dlos_fit.dlospeak_fit( 
            np.array(combined_dlos), 
            fit = 'gauss', 
            peak_range=[-15.0, 15.0]
            ) 

    return [bestfit_fpeak, bestfit_sigma, bestfit_amp]

def cmasslowz_dlospeak_fit_test(zbin): 
    """ Test dLOS distribution peaking fitting for CMASS LOWZ 
    combined sample
    """
    combined_dlos = [] 

    bestfit_fpeak, bestfit_sigma, bestfit_amp = dlos_peak_fit_combined_dlos(
            zbin,
            combined_dlos=combined_dlos
            )

    dlos_fig = Plotdlos()
    
    dlos_fig.plotfrom_dlos(
            combined_dlos, 
            binsize = 'fd_binsize', 
            lw = 4, 
            color = prettycolors()[1], 
            label = r'$\mathtt{d_{LOS}}$ Distribution' 
            )

    fit_label = ''.join([
        "Gauss ", 
        r"$(f_{peak} = ", str(round(bestfit_fpeak,2)), 
        ", \sigma = ", str(round(bestfit_sigma,2)), 
        ", A=", str(round(bestfit_amp)),
        ")$"
        ]) 

    dlos_fig.sub.plot(
            np.arange(-100.0, 100.0, 0.1), 
            peak_gauss(np.arange(-100.0, 100.0, 0.1), [bestfit_amp, bestfit_sigma]), 
            ls = '--', 
            lw = 4, 
            color = 'k',
            label = fit_label
            )

    dlos_fig.set_range(xrange=[-20.0, 20.0])
    dlos_fig.set_legend()

    fig_name = ''.join([
        'figure/', 
        'CMASSLOWZ_combinedsample_',
        zbin, 
        '_dlos_peak_fit.png'
        ])
    dlos_fig.save_fig(fig_name)

def dlospeak_sampled_dlos_test(): 
    """ Build "dlospeak" corrected CMASS LOWZ combined sample for three separate sector regions and two 
    redshfit bins
    """
    for zbin_str in ['_low', '_high']:
        dlos_fig = Plotdlos()

        if 'low' in zbin_str: 
            sigma_z = 6.9
            fpeak_z = 0.72
        elif 'high' in zbin_str: 
            sigma_z = 6.3
            fpeak_z = 0.7

        for specifier_str in ['', 'e2', 'e3']: 

            cat_corr = {
                    'catalog': {
                        'name': 'cmasslowz'+specifier_str+zbin_str
                        }, 
                    'correction': {
                        'name': 'dlospeak', 
                        'fit': 'gauss', 
                        'sigma': sigma_z, 
                        'fpeak': fpeak_z
                        } 
                    }
            dlosclass = Dlos(cat_corr)
            dlosclass.read()

            dlos_i = dlosclass.dlos 
            
            # read in sampled peak dLOS values 
            dataclass = Data('data', cat_corr)
            datafile = dataclass.file_name

            sampled_dlos = np.loadtxt(
                    datafile + '.dlos', 
                    unpack = True, 
                    usecols = [0] 
                    )

            try: 
                combined_dlos += list(dlos_i)
                combined_sampled_dlos += list(sampled_dlos)
            except NameError: 
                combined_dlos = list(dlos_i)
                combined_sampled_dlos = list(sampled_dlos)

        # histogram is rescaled due to the fact that CMASS dLOS 
        # file only contain *resolved* fiber collided pairs
        n_combined_dlos = np.float(len(combined_dlos))
        n_combined_sampled_dlos = np.float(len(combined_sampled_dlos))
        
        dlos_fig.plotfrom_dlos(
                combined_dlos, 
                binsize = 0.5, 
                lw = 4, 
                color = prettycolors()[1],
                label = 'Measured dLOS',
                rescale = n_combined_sampled_dlos/(fpeak_z * n_combined_dlos)
                )

        dlos_fig.plotfrom_dlos(
                combined_sampled_dlos, 
                binsize = 0.5, 
                lw = 2, 
                ls = '--', 
                color = 'k', 
                label = 'sampled dLOS'
                )
        dlos_fig.set_range(xrange=[-20.0, 20.0])
        dlos_fig.set_legend()

        fig_name = ''.join([
            'figure/', 
            'qaplot_sampled_dlos_', 
            'cmasslowz_combined_',
            zbin_str, 
            '.png'
            ])
        dlos_fig.save_fig(fig_name)

if __name__=="__main__": 
    build_combinedsample_random()
    #build_dlospeak_corr_combinedsample()
    #dlospeak_sampled_dlos_test()
