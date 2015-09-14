import numpy as np
import math 
import os 
from scipy.stats import norm 

import mpfit
from dlos import Dlos
from util.direc import direc
from dlos_env import DlosEnv
from fitting import mpfit_linear
from fitting import mpfit_peak_gauss
from fitting import mpfit_peak_expon
from fitting import dlos_envbin_peakfit

def dlospeak_fit_test(dlos, fit = 'gauss', peak_range = [-15.0, 15.0], **kwargs): 
    """ TEST: Fit the peak of the dLOS distribution to specified 
    fit type. Function uses MPfit. 
    """
    
    # dummy catalog correction dictionary
    dum_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'upweight'}
            }

    dlosclass = Dlos(dum_cat_corr)
    dlosclass.dlos = dlos

    # Freedman Diaconis binsize
    fd_binsize = dlosclass.fd_binsize(
            peak_range = peak_range
            )
    
    dlos_mid, dlos_dist = dlosclass.dlos_dist(
            binsize = fd_binsize, 
            **kwargs
            )
    
    # only fit the peak of the dLOS distribution  
    inpeak = np.where(
            (dlos_mid < peak_range[1]) & 
            (dlos_mid > peak_range[0])
            ) 

    amp_guess = np.mean( dlos_dist[ np.where( (dlos_mid < 1.5) & (dlos_mid > -1.5) ) ] )

    p0 = [amp_guess, 5.0]    # initial guess

    fa = {'x': dlos_mid[inpeak], 'y': dlos_dist[inpeak]}

    if fit == 'expon': 
        mpfit_func = mpfit_peak_expon
    elif fit == 'gauss': 
        mpfit_func = mpfit_peak_gauss
    else: 
        raise NotImplementedError("Fit not yet coded") 

    fit_param = mpfit.mpfit(mpfit_func, p0, functkw=fa, quiet=1)
        
    best_amp = fit_param.params[0]
    best_sigma = fit_param.params[1]

    if fit == 'gauss': 
        # analytic integral of gaussian with bestfit amplitude and sigma 
        bestfit_func_integral = np.sqrt(2.0 * np.pi) * best_amp * best_sigma 
    else: 
        raise NotImplementedError()

    dlos_dist_integral = np.sum(dlos_dist) * fd_binsize
    
    # peak fraction calculated from integral of bestfit 
    # function (should be lower than the peak fraction
    # from the distribution itself)
    fpeak_bestfit_func = bestfit_func_integral / dlos_dist_integral
        
    # peak fraction calculated based on the area of the 
    # dLOS distribution in the peak dLOS range (should be
    # higher than previous measurement)
    dlos_inpeak = np.where(
            (dlos > -3.0 * best_sigma) & 
            (dlos < 3.0 * best_sigma)
            )
    fpeak_dlos_dist = np.float(len(dlos[dlos_inpeak]))/np.float(len(dlos))  # computed from the distribution
    
    print 'fpeak from best-fit = ', fpeak_bestfit_func
    print 'fpeak from dlosdist = ', fpeak_dlos_dist  

    return [fpeak_bestfit_func, best_sigma, best_amp]

def compare_dlospeak_fit_test_to_normfit(): 
    """ Kind of stupid test to compare the dlospeak_fit function to 
    scipy.stat.norm.fit. Obviously the values for sigma are not the
    same because sigma =/= standard deviation of dLOS peak values 
    because dLOS peak is not a Gaussian. It has heavier tails. 

    
    Also tests that normlizating the dLOS distribution does NOT 
    change the sigma or fpeak fit values! MPfit does a good job. 

    """

    cat_corr = { 
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'upweight'}
            }

    dlosclass = Dlos(cat_corr)
    dlosclass.read() 

    print dlospeak_fit_test(dlosclass.dlos, fit = 'gauss', peak_range = [-15.0, 15.0])
    print dlospeak_fit_test(dlosclass.dlos, fit = 'gauss', peak_range = [-15.0, 15.0], normed=True)
    
    inpeak = np.where(
            (dlosclass.dlos < 15.0) & 
            (dlosclass.dlos > -15.0)
            ) 
    print norm.fit(dlosclass.dlos[inpeak])
    print np.std(dlosclass.dlos[inpeak])

    return None 

def dlosenv_peakfit_fpeak_env_fit_test(cat_corr, n_NN=3, fit='gauss', **kwargs): 
    """ Linear fit of the best-fit fpeak as a function of environment. 
    Bestfit fpeak values are computed for the dLOS distribution in bins of 
    galaxy environment. Linear fit is done using MPFit.

    """

    env_low, env_high, fpeaks, fpeak_errs, sigmas, sigma_errs, amps, nbins =  \
            dlos_envbin_peakfit(
                    cat_corr, 
                    n_NN = n_NN,
                    fit = fit, 
                    **kwargs
                    )
    
    env_mid = np.array( 0.5 * (env_low + env_high) )
    
    p0 = [ -0.01, 0.8 ] # guess
    fa = {'x': env_mid, 'y': fpeaks, 'err': fpeak_errs}
    
    fit_param = mpfit.mpfit(mpfit_linear, p0, functkw=fa, quiet=1)
        
    best_slope = fit_param.params[0]
    best_yint = fit_param.params[1]
    
    print 'Entire range'
    print best_slope, best_yint
    
    for range_test in [20.0, 30.0, 40.0, 50.0]: 

        test_range = np.where(
                env_mid < range_test 
                )
        fa = {'x': env_mid[test_range], 'y': fpeaks[test_range], 'err': fpeak_errs[test_range]}
        
        fit_param = mpfit.mpfit(mpfit_linear, p0, functkw=fa, quiet=1)
            
        best_slope = fit_param.params[0]
        best_yint = fit_param.params[1]
        print 'dNN < ', str(range_test)
        print best_slope, best_yint
    
    return best_slope, best_yint

def dlosenv_peakfit_sigma_env_fit_test(cat_corr, n_NN=3, fit='gauss', **kwargs): 
    """ Linear fit of the best-fit sigma as a function of environment. 
    Bestfit sigma values are computed for the dLOS distribution in bins of 
    galaxy environment. Linear fit is done using MPFit.

    """

    env_low, env_high, fpeaks, fpeak_errs, sigmas, sigma_errs, amps, nbins =  \
            dlos_envbin_peakfit(
                    cat_corr, 
                    n_NN = n_NN,
                    fit = fit, 
                    **kwargs
                    )
        
    env_mid = np.array( 0.5 * (env_low + env_high) )
    p0 = [ -0.03, 4.0 ] # guess
    fa = {'x': env_mid, 'y': sigmas, 'err': sigma_errs}
    
    fit_param = mpfit.mpfit(mpfit_linear, p0, functkw=fa, quiet=1)
        
    best_slope = fit_param.params[0]
    best_yint = fit_param.params[1]

    print 'Entire range'
    print best_slope, best_yint
    
    for range_test in [20.0, 30.0, 40.0, 50.0]: 

        test_range = np.where(
                env_mid < range_test 
                )
        fa = {'x': env_mid[test_range], 'y': sigmas[test_range], 'err': sigma_errs[test_range]}
        
        fit_param = mpfit.mpfit(mpfit_linear, p0, functkw=fa, quiet=1)
            
        best_slope = fit_param.params[0]
        best_yint = fit_param.params[1]
        
        print 'dNN < ', str(range_test)
        print best_slope, best_yint

    return best_slope, best_yint

if __name__=="__main__": 
    cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock':1}, 
            'correction': {'name': 'upweight'}
            }
    
    for nNN in [5]: 
        fpeak_slope, fpeak_yint = dlosenv_peakfit_fpeak_env_fit_test(
                cat_corr, 
                n_NN=nNN, 
                fit='gauss'
                )
        fpeak_slope, fpeak_yint = dlosenv_peakfit_sigma_env_fit_test(
                cat_corr, 
                n_NN=nNN, 
                fit='gauss'
                )
