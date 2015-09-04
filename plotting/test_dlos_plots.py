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

def dlospeak_dlos_check(cat_corr):
    ''' Compare peak sampled dLOS generated for the dLOS Peak Correction 
    (outputed during the correction) with dLOS distribution of the mock 
    
    Parameters
    ----------
    * cat_corr : catalog correction dictionary  

    '''

    dlos_fig = Plotdlos()
    
    print cat_corr
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
    print cat_corr
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
    dlos_fig.set_range()
    dlos_fig.set_legend()

    dlos_fig.show_fig()

if __name__=="__main__": 
    cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock':1}, 
            'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 4.0, 'fpeak':0.69}
            }
    dlospeak_dlos_check(cat_corr)
