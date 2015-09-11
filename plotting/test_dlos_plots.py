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
from dlos.fitting import catalog_dlospeak_env_fit
from dlos.fitting import peak_gauss 
from dlos.fitting import dlos_envbin_peakfit
from dlos.fitting import dlos_peakfit_fpeak_env_fit
from dlos.fitting import dlos_peakfit_sigma_env_fit


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
    return dlos_fig

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
    #dlos_fig.show_fig()
    
    fig_file = ''.join([
        '/home/users/hahn/powercode/FiberCollisions/figure/' , 
        catalog_name, 
        '_combined_dlos_peakfit_', 
        fit, '.png'
        ])

    dlos_fig.save_fig(fig_file)

    return [sigma, fpeak]

def catalog_dlospeak_env_fit_test(catalog_name, n_NN=3, fit='gauss', **kwargs): 
    """ Test dLOS distribution as a function of environment peaking fitting 
    """

    combined_dlos, combined_env = [], [] 
    bestfit_fpeaks, bestfit_sigmas, bestfit_amps, bestfit_envbins = catalog_dlospeak_env_fit(
            catalog_name, 
            n_NN=n_NN,
            fit='gauss', 
            combined_dlos=combined_dlos, 
            combined_env=combined_env
            )

    combined_dlos = np.array(combined_dlos)
    combined_env = np.array(combined_env)
        
    dlos_fig = Plotdlos()
    
    for i_env, envbin in enumerate(bestfit_envbins): 

        in_envbin = np.where(
                (combined_env >= envbin[0]) & 
                (combined_env < envbin[1])
                )

        envbin_label = ''.join([
            r'$\mathtt{d_{LOS}}$ Distribution $', 
            str(envbin[0]), 
            ' < d_{', str(n_NN), 'NN} < ', 
            str(envbin[1]), '$'
            ])
        
        dlos_fig.plotfrom_dlos(
                combined_dlos[in_envbin], 
                binsize = 'fd_binsize', 
                lw = 4, 
                color = prettycolors()[i_env % 15], 
                label = envbin_label
                )

        fit_label = ''.join([
            "Gauss ", 
            r"$(f_{peak} = ", str(round(bestfit_fpeaks[i_env], 2)), 
            ", \sigma = ", str(round(bestfit_sigmas[i_env], 2)), 
            ", A =", str(round(bestfit_amps[i_env])),
            ")$"
            ]) 

        dlos_fig.sub.plot(
                np.arange(-100.0, 100.0, 0.1), 
                peak_gauss(
                    np.arange(-100.0, 100.0, 0.1), 
                    [bestfit_amps[i_env], bestfit_sigmas[i_env]]
                    ), 
                ls = '--', 
                lw = 4, 
                color = prettycolors()[i_env % 15], 
                label = fit_label
                )

    dlos_fig.set_range(xrange=[-20.0, 20.0])
    #dlos_fig.set_legend()
        
    fig_file = ''.join([
        '/home/users/hahn/powercode/FiberCollisions/figure/' , 
        'qaplot_', 
        catalog_name, 
        '_combined_dlos_d', str(n_NN), 'NN', 
        '_peakfit_', fit, 
        '.png'
        ])

    dlos_fig.save_fig(fig_file)
    plt.close()

    return None

def dlos_envbin_peakfit_test(cat_corr, n_NN=3, **kwargs): 
    """
    """
    catdict = cat_corr['catalog']
    
    prettyplot()

    fpeak_fig = plt.figure(1)
    fpeak_sub = fpeak_fig.add_subplot(111)
    sigma_fig = plt.figure(2)
    sigma_sub = sigma_fig.add_subplot(111)

    if not isinstance(n_NN, list): 
        n_NN_list = [n_NN]
    else: 
        n_NN_list = n_NN
        
    for i_nNN, nNN in enumerate(n_NN_list): 

        env_low, env_high, fpeaks, sigmas, amps, nbins =  dlos_envbin_peakfit(
                cat_corr, 
                n_NN = nNN,
                **kwargs
                )

        env_mid = 0.5 * ( env_low + env_high )

        fpeak_errs = np.sqrt( fpeaks / nbins )
        sigma_errs = np.sqrt( sigmas / nbins ) 

        #fpeak_sub.plot(env_mid, fpeaks, lw = 4, c = prettycolors()[i_nNN], label = r"$\mathtt{d_{"+str(nNN)+"NN}}$")
        fpeak_sub.errorbar(
                env_mid, fpeaks, 
                yerr = fpeak_errs, 
                c = prettycolors()[i_nNN],
                label = r"$\mathtt{d_{"+str(nNN)+"NN}}$"
                )
        # best fit line to fpeak(env)
        fpeak_slope, fpeak_yint = dlos_peakfit_fpeak_env_fit(cat_corr, n_NN=nNN, fit='gauss', **kwargs)
        print fpeak_slope, fpeak_yint
        fpeak_sub.plot(
                env_mid, fpeak_slope * env_mid + fpeak_yint,
                lw = 4, 
                ls = '--', 
                c = prettycolors()[i_nNN]
                )

        #sigma_sub.plot(env_mid, sigmas, lw=4, c=prettycolors()[i_nNN], label=r"$\mathtt{d_{"+str(nNN)+"NN}}$")
        sigma_sub.errorbar(
                env_mid, sigmas, 
                yerr = sigma_errs, 
                c = prettycolors()[i_nNN], 
                label = r"$\mathtt{d_{"+str(nNN)+"NN}}$"
                )
        # best fit line to sigma(env)
        sigma_slope, sigma_yint = dlos_peakfit_sigma_env_fit(cat_corr, n_NN=nNN, fit='gauss', **kwargs)
        print sigma_slope, sigma_yint
        sigma_sub.plot(
                env_mid, sigma_slope * env_mid + sigma_yint,
                lw = 4, 
                ls = '--', 
                c = prettycolors()[i_nNN]
                )

    fpeak_sub.set_xlabel(''.join([r"$\mathtt{d_{", str(n_NN), "NN}}$"]), fontsize = 20)
    fpeak_sub.set_ylabel(''.join([r"$\mathtt{f_{peak}}$"]), fontsize = 20)
    fpeak_sub.set_xlim([0.0, env_high[-1]])
    fpeak_sub.set_ylim([0.0, 1.0]) 
    fpeak_sub.legend()

    sigma_sub.set_xlabel(''.join([r"$\mathtt{d_{", str(n_NN), "NN}}$"]), fontsize = 20)
    sigma_sub.set_ylabel(''.join([r"$\mathtt{\sigma}$"]), fontsize = 20)
    sigma_sub.set_xlim([0.0, env_high[-1]])
    sigma_sub.set_ylim([0.0, 5.0])
    sigma_sub.legend()
    
    fpeak_fig_file = ''.join([
        '/home/users/hahn/powercode/FiberCollisions/figure/' , 
        'qaplot_', 
        catdict['name'], 
        '_fpeak_d', ','.join([str(nNN) for nNN in n_NN_list]), 'NN',
        '.png'
        ])
    fpeak_fig.savefig(fpeak_fig_file, bbox_inches="tight")

    sigma_fig_file = ''.join([
        '/home/users/hahn/powercode/FiberCollisions/figure/' , 
        'qaplot_', 
        catdict['name'], 
        '_sigma_d', ','.join([str(nNN) for nNN in n_NN_list]), 'NN',
        '.png'
        ])
    sigma_fig.savefig(sigma_fig_file, bbox_inches="tight")

    plt.close()

    return None


if __name__=="__main__": 
    cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock':1}, 
            'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 4.0, 'fpeak':0.69}
            }

    #for nNN in [1,3,5,7,10]: 
    #    catalog_dlospeak_env_fit_test('nseries', n_NN=nNN, fit='gauss') 
    dlos_envbin_peakfit_test(cat_corr, n_NN=[3,5,7,10])
    #catalog_dlospeak_fit_test('nseries', fit='gauss')
    #dlospeak_dlos_test(cat_corr)
