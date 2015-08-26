'''

Plotting for fiber collisions with photoz and environment data

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
import fibcol_dlos as fc_dlos
import fibcol_utility as fc_util
import galaxy_environment as genv
import dlos_env
import mpfit as mpfit

def plot_dLOS_env(cat_corr, n=3, **kwargs):
    ''' Plot combined dLOS distribution of in bins of galaxy environment 
    
    Notes
    -----
    * Only implemented for Nth Nearest Neighbor galaxy environment tracer 

    '''
    # combined dLOS for catalog and correction 
    comb_dlos = dlos_env.dlos_env(cat_corr, n=n, **kwargs) 

    combined_dlos = comb_dlos.dlos
    combined_dNN = comb_dlos.env

    # set up figure 
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=(14,5))
    sub = fig.add_subplot(1,1,1)
    
    print 'd'+str(n)+'NN, Minimum ', min(combined_dNN), ' Maximum ', max(combined_dNN)

    # bin dNN values
    dnn_step = 2    # currently hardcoded
    istart = np.int(np.floor(min(combined_dNN)))
    iend = np.int(np.floor(max(combined_dNN)/np.float(dnn_step)))
    dNN_bins = [ (dnn_step*i, dnn_step*(i+1)) for i in np.arange(istart, iend) ] 

    for i_bin, dNN_bin in enumerate(dNN_bins): 
        # dNN bin 
        try: 
            bin_index = np.where(( combined_dNN >= dNN_bin[0] ) & ( combined_dNN < dNN_bin[1] )) 
            bin_label = ''.join([r'$', str(dNN_bin[0]), ' < d_{', str(n), 'NN} < ', str(dNN_bin[1])]) 
        except TypeError: 
            bin_index = np.where((combined_dNN >= dNN_bin)) 
            bin_label = ''.join([r'$', str(dNN_bin), ' < d_{', str(n), 'NN}']) 

        bin_dlos = combined_dlos[bin_index] 
        if len(bin_dlos) < 50: 
            # if the bin contains too few dLOSs then skip
            continue
        
        # fraction of total dLOS values 
        bin_perc = np.float( len(bin_dlos) )/np.float( len(combined_dlos) ) * 100.0
    
        # calculate the histogram and fit peak of the distribution 
        dlos_hist, mpc_mid, peak_param  = \
                fc_dlos.dlos_hist_peak_fit(bin_dlos, fit='gauss', peak_range=[-15.0, 15.0])
        
        # plot dLOS distribution 
        sub.plot(mpc_mid, dlos_hist, 
                lw=4, color=pretty_colors[((i_bin+1) % 20)], 
                label = bin_label+',\;( '+('%.1f' % bin_perc)+'\%) $') 
        
        # plot best-fit to dLOS 
        fit_label = r'$\sigma ='+('%.2f' % peak_param['sigma'])+\
                ', f_{peak} = '+('%.2f' % peak_param['fpeak'])+'$'
        sub.plot(mpc_mid, fc_dlos.peak_gauss(mpc_mid, [peak_param['amp'], peak_param['sigma']]), 
                lw=4, ls='--', color=pretty_colors[((i_bin+1) % 20)], label=fit_label)

    sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
    sub.set_xlim([-50.0, 50.0])
    if len(dNN_bins) < 5: 
        sub.legend(loc='upper left') 
    sub.set_ylim([0.0, 0.1]) 
    fig_file = ''.join(['figure/', 
        'dlos_env_dependence_d', str(n), 'NN_', catalog['name'], '.png']) 

    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear()
    
    return None 

def plot_dLOS_fpeak_env(cat_corr, n=3, **kwargs):
    ''' Plot best-fit fpeak of combined dLOS distribution of in bins of galaxy environment 
    
    Notes
    -----
    * Only implemented for Nth Nearest Neighbor galaxy environment tracer 

    '''
    catalog = cat_corr['catalog']
    # combined dLOS for catalog and correction 
    comb_dlos = dlos_env.dlos_env(cat_corr, n=n, **kwargs) 

    combined_dlos = comb_dlos.dlos
    combined_dNN = comb_dlos.env

    # set up figure 
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=(8,8))
    sub = fig.add_subplot(1,1,1)
    
    print 'd'+str(n)+'NN, Minimum ', min(combined_dNN), ' Maximum ', max(combined_dNN)

    # bin dNN values
    dnn_step = 2    # currently hardcoded
    istart = np.int(np.floor(min(combined_dNN)))
    iend = np.int(np.floor(max(combined_dNN)/np.float(dnn_step)))
    dNN_bins = [ (dnn_step*i, dnn_step*(i+1)) for i in np.arange(istart, iend) ] 

    for i_bin, dNN_bin in enumerate(dNN_bins): 
        # dNN bin 
        try: 
            bin_index = np.where(( combined_dNN >= dNN_bin[0] ) & ( combined_dNN < dNN_bin[1] )) 
        except TypeError: 
            bin_index = np.where((combined_dNN >= dNN_bin)) 

        bin_dlos = combined_dlos[bin_index] 
        if len(bin_dlos) < 50: 
            # if the bin contains too few dLOSs then skip
            continue
        
        # fraction of total dLOS values 
        bin_perc = np.float( len(bin_dlos) )/np.float( len(combined_dlos) ) * 100.0
    
        # calculate the histogram and fit peak of the distribution 
        dlos_hist, mpc_mid, peak_param  = \
                fc_dlos.dlos_hist_peak_fit(bin_dlos, fit='gauss', peak_range=[-15.0, 15.0])
        
        try:        # save avg_dNN and fpeak values 
            dNN_avg.append( 0.5 * np.float(dNN_bin[0] + dNN_bin[1]) ) 
            fpeaks.append( peak_param['fpeak'] ) 
        except NameError: 
            dNN_avg = [0.5 * np.float(dNN_bin[0] + dNN_bin[1])]
            fpeaks = [peak_param['fpeak']]
        except TypeError: 
            dNN_avg.append( np.float(dNN_bin) ) 
            fpeaks.append( peak_param['fpeak'] ) 

    sub.scatter( dNN_avg, fpeaks, s=6 )       # fpeaks vs dNN-avg

    # MPfit ----
    p0 = [-0.01, 0.8]            # initial guesses for p0[0] * x + p0[1]
    fa = {'x': np.array(dNN_avg), 'y': np.array(fpeaks)}
    fit_param = mpfit.mpfit(fc_dlos.mpfit_linear, p0, functkw=fa, nprint=0)
    
    fit_slope = fit_param.params[0]
    fit_yint = fit_param.params[1]

    # plot best line fit 
    sub.plot( np.array(dNN_avg), fc_dlos.fit_linear(np.array(dNN_avg), fit_param.params), lw=4, ls='--', c='k')       

    sub.set_xlabel('$\mathtt{d_{'+str(n)+'NN}}$', fontsize=20) 
    sub.set_ylabel('$\mathtt{f_{peak}}$', fontsize=20) 
    sub.set_xlim([0.0, 50.0]) 
    sub.set_ylim([0.0, 1.0])

    fig_file = ''.join(['../figure/', 
        'dlos_fpeak_d', str(n), 'NN_', catalog['name'], '.png'])
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear()

    return None 

if __name__=='__main__':
    cat_corr = {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}}
    plot_dLOS_fpeak_env(cat_corr, n=3)
