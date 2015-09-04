'''

dLOS distribution plots 


'''
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt
import os.path
from matplotlib.collections import LineCollection

# --- Local --- 
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from dlos.dlos import Dlos

class Plotdlos(object): 

    def __init__(self, **kwargs): 
        """ Class that describes dLOS distribution plots 
        """    
        self.kwargs = kwargs

        # pretty figure up
        prettyplot() 
        pretty_colors = prettycolors()  
    
        self.fig = plt.figure(1) 
        self.sub = self.fig.add_subplot(111) 

        self.hist_max = 0.0 
   
    def plotfrom_catcorr(self, cat_corr, **pltkwargs): 
        """ dLOS distribution from catalog correction dictionary
        """

        catdict = cat_corr['catalog']
        corrdict = cat_corr['correction']

        dlosclass = Dlos(cat_corr)
        dlos_file = dlosclass.file_name 

        if 'binsize' in pltkwargs.keys(): 
            binsize = pltkwargs['binsize']
            pltkwargs.pop('binsize', None) # remove from dictionary
        else: 
            binsize = 0.5   # (default)

        xmid, dlos_hist = dlosclass.dlos_dist( binsize = binsize )
    
        if 'label' not in pltkwargs.keys(): 
            pltkwargs['label'] = ''.join([
                catdict['name'], ':', 
                corrdict['name']
                ])
        elif pltkwargs['label'] == False: 
            pass

        self.sub.plot(xmid, dlos_hist, **pltkwargs) 
    
        self.hist_max = max([ dlos_hist.max(), self.hist_max ]) 

        return None 

    def plotfrom_dlos(self, dlosdata, **pltkwargs): 
        """ dLOS distribution from dLOS data  
        """

        # arbitrary catalog correction dictionary 
        cat_corr = {
                'catalog': {'name': 'nseries', 'n_mock': 1}, 
                'correction': {'name': 'upweight'} 
                }

        dlosclass = Dlos(cat_corr)
        dlosclass.dlos = dlosdata
        
        if 'binsize' in pltkwargs.keys(): 
            binsize = pltkwargs['binsize']
            pltkwargs.pop('binsize', None) # remove from dictionary
        else: 
            binsize = 0.5   # (default)
        
        xmid, dlos_hist = dlosclass.dlos_dist( binsize = binsize ) 
        
        self.sub.plot(xmid, dlos_hist, **pltkwargs) 
        
        self.hist_max = max([ dlos_hist.max(), self.hist_max ]) 
        
        return None 

    def set_range(self, **rangekwargs): 
        """ Set x and y range for the plot 
        """

        if 'xrange' in rangekwargs.keys(): 
            xrange = rangekwargs['xrange']
        else: 
            xrange = [-50.0, 50.0] # (default)

        if 'yrange' in rangekwargs.keys(): 
            yrange = rangekwargs['yrange']
        else: 
            yrange = [0.0, 1.25 * self.hist_max]

        self.sub.set_xlim(xrange) 
        self.sub.set_ylim(yrange) 

        self.sub.set_xlabel(r"$\mathtt{d_{LOS}}$ (Mpc/h)", fontsize=20)

        return None 

    def set_legend(self, **lgdkwargs): 
        """ Set the legend of the figure 
        """

        if 'loc' not in lgdkwargs.keys(): 
            lgdkwargs['loc'] = 'upper right'
        
        if 'scatterpoints' not in lgdkwargs.keys(): 
            lgdkwargs['scatterpoints'] = 1 

        self.sub.legend(**lgdkwargs) 
 
        return None 
    
    def save_fig(self, fig_name): 
        """ Save figure
        """
        self.fig.savefig(fig_name, bbox_inches="tight")
        self.fig.clear() 

        return None 

    def show_fig(self): 
        """ Show figure
        """
        plt.show()


