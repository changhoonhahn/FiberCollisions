'''

dLOS distribution plots 


'''
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt
import os.path
from matplotlib.collections import LineCollection

# --- Local --- 
from utility.plotting import prettyplot
from utility.plotting import prettycolors 
from dlos.dlos import Dlos

class Plotdlos(object): 

    def __init__(self, **kwargs): 
        """ Class that describes dLOS distribution plots 
        """    
        self.kwargs = kwargs

        # pretty figure up
        prettyplot() 
        pretty_colors = prettycolors()  
    
        fig = plt.figure(1) 
        sub = fig.add_subplot(111) 
   
    def from_catcorr(self, cat_corr, **pltkwargs): 
        """ dLOS distribution from catalog correction dictionary
        """

        catdict = cat_corr['catalog']
        corrdict = cat_corr['correction']

        dlosclass = Dlos(cat_corr)
        dlos_file = dlosclass.file_name 

        xmid, dlos_hist = dlosclass.dlos_dist( binsize = 0.5 )
    
        if 'label' not in pltkwargs.keys(): 
            pltkwargs['label'] = #INSERT SOMETHING!

        sub.plot(xmid, dlos_hist, label=catcorr_label, **pltkwargs) 

    def from_dlos(self, dlosdata, **pltkwargs): 
        """ dLOS distribution from dLOS data  
        """
        
        dlosclass = Dlos(cat_corr)
        dlosclass.dlos = dlosdata
        
        xmid, dlos_hist = dlosclass.dlos_dist( binsize = 0.5 ) 

        sub.plot(xmid, dlos_hist, label=dlos_label, **pltkwargs) 

    for i_cat, cat in enumerate(cat_list): 
        cat_corr = {'catalog': cat, 'correction': None}
        # read in dLOS values for cat_corr 
        los_disp_i = fc_dlos.dlos(**cat_corr)
        # calculate dLOS histogram 
        dlos_hist, mpc_binedges = np.histogram(los_disp_i.dlos, bins=n_bins, range=[x_min, x_max])
        if np.sum(peakdlos_hist)/np.sum(dlos_hist) > 3.: 
            # plot dLOS distribution of mock 
            sub.plot(xmid, (1.0/fpeak_list[i_cat])*(np.sum(peakdlos_hist)/np.sum(dlos_hist))*dlos_hist, 
                    lw=2, ls='--', color='k', label='dLOS of Mock Catalog') 
        else: 
            # plot dLOS distribution of mock 
            sub.plot(xmid, dlos_hist, 
                    lw=2, ls='--', color='k', label='dLOS of Mock Catalog') 
    sub.set_xlim([-20., 20.])
    sub.set_xlabel(r"$d_{LOS}$ Mpc") 
    sub.legend() 
 


