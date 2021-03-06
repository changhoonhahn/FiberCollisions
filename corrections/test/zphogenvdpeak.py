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
import plot_fibcol as fc_plot
import galaxy_environment as genv
import dlos_env
import mpfit as mpfit

# LOS Displacement ------
def plot_dLOS_env(cat_corr, n_NN=3, **kwargs):
    ''' Plot combined dLOS distribution of in bins of galaxy environment 
    
    Notes
    -----
    * Only implemented for Nth Nearest Neighbor galaxy environment tracer 

    '''
    # combined dLOS for catalog and correction 
    comb_dlos = dlos_env.DlosEnv(cat_corr, n_NN=n_NN, **kwargs) 
    combined_dlos = comb_dlos.dlos
    combined_dNN = comb_dlos.env

    # set up figure 
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=(14,5))
    sub = fig.add_subplot(1,1,1)
    
    print 'd'+str(n_NN)+'NN, Minimum ', min(combined_dNN), ' Maximum ', max(combined_dNN)

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

def plot_dLOS_fpeak_env(cat_corr, n_NN=3, **kwargs):
    ''' Plot best-fit fpeak of combined dLOS distribution of in bins of galaxy environment 
    
    Notes
    -----
    * Only implemented for Nth Nearest Neighbor galaxy environment tracer 

    '''
    if not isinstance(n_NN, list): 
        n_NN_list = [ n_NN ] 
    else:
        n_NN_list = n_NN
    
    # set up figure 
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=(8,8))
    sub = fig.add_subplot(1,1,1)

    catalog = cat_corr['catalog']
    for i_nNN, nNN in enumerate(n_NN_list):   # for different n_NN values 
        # combined dLOS for catalog and correction 
        comb_dlos = dlos_env.DlosEnv(cat_corr, n_NN=nNN, **kwargs) 
        comb_dlos.Read()

        combined_dlos = comb_dlos.dlos
        combined_dNN = comb_dlos.env

        print 'd'+str(nNN)+'NN, Minimum ', min(combined_dNN), ' Maximum ', max(combined_dNN)
        dNN_label = ''.join(['d', str(nNN), 'NN'])
    
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

        sub.scatter( dNN_avg, fpeaks, s=6, c=pretty_colors[i_nNN])       # fpeaks vs dNN-avg

        # MPfit ----
        p0 = [-0.01, 0.8]            # initial guesses for p0[0] * x + p0[1]
        fa = {'x': np.array(dNN_avg), 'y': np.array(fpeaks)}
        fit_param = mpfit.mpfit(fc_dlos.mpfit_linear, p0, functkw=fa, nprint=0)
        
        fit_slope = fit_param.params[0]
        fit_yint = fit_param.params[1]

        dNN_label += ': '+str(round(fit_slope, 4))+', '+str(round(fit_yint, 4))
        # plot best line fit 
        sub.plot( np.array(dNN_avg), fc_dlos.fit_linear(np.array(dNN_avg), fit_param.params), 
                lw=4, ls='--', c=pretty_colors[i_nNN], label=dNN_label)       
        
        del dNN_avg, fpeaks

    sub.set_xlabel('$\mathtt{d_{NN}}$', fontsize=20) 
    sub.set_ylabel('$\mathtt{f_{peak}}$', fontsize=20) 
    sub.set_xlim([0.0, 50.0]) 
    sub.set_ylim([0.0, 1.0])
    sub.legend(loc='upper right')

    fig_file = ''.join(['../figure/', 
        'dlos_fpeak', ''.join(['_d'+str(ni)+'NN' for ni in n_NN_list]), '_', catalog['name'], '.png'])
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear()

    return None 

def plot_dLOS_envdist(cat_corr, n_NN=3, **kwargs): 
    ''' Plot environment distribution of dLOS 
    '''
    if not isinstance(n_NN, list): 
        n_NN_list = [ n_NN ] 
    else:
        n_NN_list = n_NN
    
    # set up figure 
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=(8,8))
    sub = fig.add_subplot(1,1,1)

    catalog = cat_corr['catalog']
    for i_nNN, nNN in enumerate(n_NN_list):   # for different n_NN values 
        # combined dLOS for catalog and correction 
        comb_dlos = dlos_env.DlosEnv(cat_corr, n_NN=nNN, **kwargs) 
        comb_dlos.Read()

        combined_dlos = comb_dlos.dlos
        combined_dNN = comb_dlos.env

        dNN_label = ''.join(['d', str(nNN), 'NN'])

        dNN_min, dNN_max = min(combined_dNN), max(combined_dNN)
        if 'stepsize' in kwargs.keys():
            stepsize = kwargs['stepsize']
        else: 
            stepsize = 2.0
        n_bins = int((dNN_max - dNN_min)/stepsize)

        dNN_dist, dNN_binedges = np.histogram(combined_dNN, bins=n_bins, range=[dNN_min, dNN_max]) 

        dNN_low = dNN_binedges[:-1]
        dNN_high = dNN_binedges[1:]
        dNN_mid = np.array([0.5 * (dNN_low[i] + dNN_high[i]) for i in range(len(dNN_low))]) 
    
        sub.step(dNN_low, dNN_dist, lw=4, color=pretty_colors[i_nNN], label = dNN_label)

    sub.set_xlabel('$\mathtt{d_{NN}}$', fontsize=20) 
    sub.set_ylabel(r'$\mathtt{N_{gal}}$', fontsize=20) 
    sub.set_xlim([0.0, 50.0]) 
    sub.legend(loc='upper right')

    fig_file = ''.join(['../figure/', 
        'dlos_envdist', ''.join(['_d'+str(ni)+'NN' for ni in n_NN_list]), '_', catalog['name'], '.png'])
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear()

    return None 

def plot_dLOS_env_fpeakdist(cat_corr, n_NN=3, **kwargs): 
    ''' Plot fpeak distribution of dLOS + Env

    Notes
    -----
    * Function designed to investigate which n value for nth nearest neighbor has the greatest variance 
    
    '''
    if not isinstance(n_NN, list): 
        n_NN_list = [ n_NN ] 
    else:
        n_NN_list = n_NN
    
    # set up figure 
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=(8,8))
    sub = fig.add_subplot(1,1,1)

    catalog = cat_corr['catalog']
    for i_nNN, nNN in enumerate(n_NN_list):   # for different n_NN values 
        # combined dLOS for catalog and correction 
        comb_dlos = dlos_env.DlosEnv(cat_corr, n_NN=nNN, **kwargs) 
        comb_dlos.Read()

        combined_dNN = comb_dlos.env
        
        # calculate fpeak as a function dNN 
        fpeaks = dlos_env.fpeak_dNN(combined_dNN, cat_corr, n_NN=nNN)

        fpeak_label = ''.join(["d", str(nNN), 'NN'])

        fpeak_min, fpeak_max = min(fpeaks), max(fpeaks)
        if 'stepsize' in kwargs.keys():
            stepsize = kwargs['stepsize']
        else: 
            stepsize = 0.025 
        n_bins = int((fpeak_max - fpeak_min)/stepsize)

        fpeak_dist, fpeak_binedges = np.histogram(fpeaks, bins=n_bins, range=[fpeak_min, fpeak_max]) 

        fpeak_low = fpeak_binedges[:-1]
        fpeak_high = fpeak_binedges[1:]
        fpeak_mid = np.array([0.5 * (fpeak_low[i] + fpeak_high[i]) for i in range(len(fpeak_low))]) 
    
        fpeak_label += '; sigma = '+str(round(np.std(fpeaks), 4))
        sub.step(fpeak_low, fpeak_dist, lw=4, color=pretty_colors[i_nNN], label = fpeak_label)

    sub.set_xlabel(r"$\mathtt{f_{peak}(dNN)}$", fontsize=20)
    sub.set_ylabel(r'$\mathtt{N_{gal}}$', fontsize=20) 
    sub.set_xlim([0.0, 1.0]) 
    sub.legend(loc='upper right')

    fig_file = ''.join(['../figure/', 
        'dlos_env_fpeakdist', ''.join(['_d'+str(ni)+'NN' for ni in n_NN_list]), '_', catalog['name'], '.png'])
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear()

    return None 

# Powerspectrum -------
def plot_photozenv_pk(mock_name, **kwargs): 
    ''' Plot powerspectrum for PhotoZ + Galaxy Environment + Peak + Shot Noise correction method. 

    Notes
    -----
    * Wrapper that uses plot_fibcol.py from FiberCollisions
    * Optimal correction methods are hardcoded for each mock catalog

    '''
    if 'Ngrid' in kwargs.keys():        # Ngrid
        Ngrid = kwargs['Ngrid']
    else: 
        Ngrid = 360     # default Ngrid

    if 'quad' in kwargs.keys(): 
        quad = kwargs['quad']
    else: 
        quad = False    # default is monopole

    if mock_name == 'nseries':          # Nseries
        peakshot_corr = {'name': 'peakshot', 'sigma': '4.0', 'fpeak': '0.69', 'fit':'gauss'}
        photozpeakshot_corr = {'name': 'photozpeakshot', 'sigma': '4.0', 'fpeak': '0.69', 'fit':'gauss'} 
        photozenvpeakshot_corr = {'name': 'photozenvpeakshot', 'fit': 'gauss', 'sigma': 4.0, 'fpeak': 0.69, 'n_NN': 5}

        # list of corrections
        catcorr_methods = [
                {'catalog': {'name': mock_name}, 'correction': {'name': 'true'}}, 
                {'catalog': {'name': mock_name}, 'correction': {'name': 'upweight'}}, 
                {'catalog': {'name': mock_name}, 'correction': peakshot_corr}, 
                {'catalog': {'name': mock_name}, 'correction': photozpeakshot_corr}, 
                {'catalog': {'name': mock_name}, 'correction': photozenvpeakshot_corr}
                ]
                #{'catalog': {'name': mock_name}, 'correction': {'name': 'scratch_peakknown'}},
        n_mock_list = 20
    else: 
        raise NotImplementedError("Not Implemented")
        
    fc_plot.plot_pk_fibcol_comp(catcorr_methods, n_mock_list, 
            quad=quad, Ngrid=Ngrid, type='regular', 
            xrange=[0.001, 1.0], yrange=[10**3, 3*10**5])
    
    fc_plot.plot_pk_fibcol_comp(catcorr_methods, n_mock_list, 
            quad=quad, Ngrid=Ngrid, type='ratio', 
            xrange=[0.001, 1.0], yrange=[0.0, 2.0])

    return None 

if __name__=='__main__':
    plot_photozenv_pk('nseries', Ngrid=960) 
