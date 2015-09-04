'''

Plots to test code/analysis of dLOS 

Author(s): ChangHoon Hahn

'''
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt
import os.path
from matplotlib.collections import LineCollection

# --- Local --- 
from spec.data import Data
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from dlos_plot import Plotdlos
from dlos.fitting import catalog_dlospeak_fit
from dlos.fitting import peak_gauss 


def dlospeak_dlos_test(cat_corr):
    ''' Compare peak sampled dLOS generated for the dLOS Peak Correction 
    (outputed during the correction) with dLOS distribution of the mock 
    
    Parameters
    ----------
    * cat_corr : catalog correction dictionary  

    '''

    dlos_fig = Plotdlos()
    
    # plot dLOS distribution from dLOS file for 
    # data corresponding to catalog correction  
    # dictionary
    dlos_fig.plotfrom_catcorr(
            cat_corr, 
            binsize = 0.5, 
            lw = 4, 
            color = prettycolors()[1]
            )

    # read in sampled peak dLOS values 
    dataclass = Data('data', cat_corr)
    datafile = dataclass.file_name

    sampled_dlos = np.loadtxt(
            datafile + '.dlos', 
            unpack = True, 
            usecols = [0] 
            )

    dlos_fig.plotfrom_dlos(
            sampled_dlos, 
            binsize = 0.5, 
            lw = 2, 
            ls = '--', 
            color = 'k', 
            label = 'sampled dLOS' 
            )
    dlos_fig.set_range(xrange=[-20.0, 20.0])
    dlos_fig.set_legend()

    dlos_fig.show_fig()

def catalog_dlospeak_fit_test(catalog_name, fit='gauss', **kwargs): 
    """ Test dLOS distribution peaking fitting 
    """
    combined_dlos = [] 
    dlos_hist = []
    bestfit_fpeak, bestfit_sigma, bestfit_amp = catalog_dlospeak_fit(
            catalog_name, 
            fit='gauss', 
            combined_dlos=combined_dlos, 
            dlos_hist=dlos_hist
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
    dlos_fig.show_fig()

"""
        fig_dir = 'figure/'
        fig_file = ''.join([fig_dir, 
            catalog['name'].lower(), '_', str(n_mocks), 
            'mocks_combined_dlos_peakfit_', fit.lower(), bigfc_flag, '.png'])
        fig.savefig(fig_file, bbox_inches="tight")
        fig.clear() 

    return [sigma, fpeak]
"""

if __name__=="__main__": 
    cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock':1}, 
            'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 4.0, 'fpeak':0.69}
            }
    catalog_dlospeak_fit_test('nseries', fit='gauss')
    #dlospeak_dlos_test(cat_corr)
