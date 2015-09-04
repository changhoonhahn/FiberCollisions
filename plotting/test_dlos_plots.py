'''

Plots to test code/analysis 

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

def dlospeak_dlos_check(cat_corrs):
    ''' Compare peak sampled dLOS generated for the Peak Correction Methods (outputed during 
    the correction) with dLOS distribution of the mock 
    
    Parameters
    ----------
    * cat_corrs : list of Catalog + Correction dictionaries

    Notes
    -----
    * Function entirely edited Aug 18, 2015 

    '''
    # pretty figure up
    prettyplot() 
    pretty_colors = prettycolors()  
    
    # set up figure 
    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 
    
    cat_list, fpeak_list = [], [] 
    for i_cat_corr, cat_corr in enumerate(cat_corrs): # loop through cat_corr dictionaries
        # catalog and correction 
        catalog = cat_corr['catalog']
        correction = cat_corr['correction'] 
        if catalog not in cat_list: 
            # keep track of all the catalogs used
            cat_list.append(catalog)
            fpeak_list.append(correction['fpeak'])
        
        corr_label = 'fpeak'+str(correction['fpeak'])+'.sigma'+str(correction['sigma'])    
        print corr_label

        # read in sampled peak dLOS values 
        data_file = fc_data.get_galaxy_data_file('data', **cat_corr)
        peakdlos_file = ''.join([data_file, '.dlosvalues']) 
        peakdlos = np.loadtxt(peakdlos_file, unpack=True, usecols=[0]) 
            
        # calculate dLOS histogram 
        x_min, x_max = -1000.0, 1000.0
        n_bins = int((x_max-x_min)/0.5) 
        peakdlos_hist, mpc_binedges = np.histogram(peakdlos, bins=n_bins, range=[x_min, x_max]) 
        # x values for plotting 
        xlow = mpc_binedges[:-1]
        xhigh = mpc_binedges[1:] 
        xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
        
        # plot sampled peak dLOS distribution  
        sub.plot(xmid, peakdlos_hist, lw=4, color=pretty_colors[i_cat_corr], label=corr_label) 

        # read in dLOS values for cat_corr 
        los_disp_i = fc_dlos.dlos(**cat_corr)
        # calculate dLOS histogram 
        dlos_hist, mpc_binedges = np.histogram(los_disp_i.dlos, bins=n_bins, range=[x_min, x_max])
        # plot dLOS distribution of mock 

        inpeak = np.where((xmid > -10.0) & (xmid < 10.0)) 
        #print dlos_hist[inpeak] - peakdlos_hist[inpeak]
        #print np.sum(dlos_hist[inpeak] - peakdlos_hist[inpeak])
    
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
    plt.show() 

