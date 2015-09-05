'''

Make dLOS distribution peak corrected CMASS and LOWZ 
combined galaxy catalogs 

'''
from spec.data import Data
from dlos.dlos import Dlos
from dlos import fitting as dlos_fit
from plotting.dlos_plot import Plotdlos
from dlos.fitting import peak_gauss 
    
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

if __name__=="__main__": 
    cmasslowz_dlospeak_fit_test('low')
    cmasslowz_dlospeak_fit_test('high')
