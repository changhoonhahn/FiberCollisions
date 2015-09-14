"""

Fitting the peak of the line-of-sight displacement distribution


"""
import numpy as np
import math 
import os 

import mpfit
from dlos import Dlos
from util.direc import direc
from dlos_env import DlosEnv

def dlospeak_fit(dlos, fit = 'gauss', peak_range = [-15.0, 15.0], **kwargs): 
    """ Fit the peak of the dLOS distribution to specified 
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
            binsize = fd_binsize
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

    fit_param = mpfit.mpfit(mpfit_func, p0, functkw=fa)
        
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

def catalog_dlospeak_fit(catalog_name, fit='gauss', **kwargs): 
    """ Fit peak of the dLOS distribution for catalog

    Parameters
    ----------

    Returns 
    ------- 
    sigma : sigma value of gaussian/exponential function  
    fpeak : peak fraction of the dLOS distribution 

    """
    
    # depends on catalog 
    if catalog_name in ('nseries'): 
        n_mocks = 84 
        catdict_list = [ 
                {'name': catalog_name, 'n_mock': i_mock} for i_mock in range(1, n_mocks+1)
                ]
        corrdict = {'name': 'upweight'}

    if 'combined_dlos' in kwargs.keys(): 
        combined_dlos = kwargs['combined_dlos'] 

    for catdict in catdict_list: 

        i_cat_corr = { 
                'catalog': catdict, 
                'correction': corrdict
                }

        dlosclass = Dlos(i_cat_corr)
        dlosclass.read() 

        dlos_i = dlosclass.dlos

        # Combine dLOS values 
        try: 
            combined_dlos += list(dlos_i)
        except NameError: 
            combined_dlos = list(dlos_i)

    fpeak, sigma, amp = dlospeak_fit(np.array(combined_dlos), fit = 'gauss', peak_range = [-15.0, 15.0])   

    return fpeak, sigma, amp

def catalog_dlospeak_env_fit(catalog_name, n_NN=3, fit='gauss', writeout=True, **kwargs): 
    """ Fit peak of the dLOS distribution as a function of environment 
    for various catalogs. dLOS distribution is calculated for bins of 
    galaxy environment. For each galaxy environment bin, dLOS distributions 
    are individually fit using MPFIT 

    Parameters
    ----------
    catalog_name : 'nseries'
    n_NN : n in n-th nearest neighbor 
    fit : 'gauss' best fit function type 

    Returns 
    ------- 
    sigma : sigma value of gaussian/exponential function  
    fpeak : peak fraction of the dLOS distribution 

    """
    
    if 'dNN_stepsize' in kwargs.keys(): 
        dNN_stepsize = kwargs['dNN_stepsize']
    else:
        dNN_stepsize = 2.0
    
    # this is to pass combined_dlos and combined_env as a 
    # global variable when this function is called
    if 'combined_dlos' in kwargs.keys(): 
        combined_dlos = kwargs['combined_dlos'] 
    
    if 'combined_env' in kwargs.keys(): 
        combined_env = kwargs['combined_env'] 
    
    # for each catalog set list of catalog correction 
    # dictionaries which are used to read in dLOS files
    if catalog_name in ('nseries'): 
        n_mocks = 84 
        catdict_list = [ 
                {'name': catalog_name, 'n_mock': i_mock} 
                for i_mock in range(1, n_mocks+1)
                ]
        corrdict = {'name': 'upweight'}
    
    for catdict in catdict_list: 

        i_cat_corr = { 
                'catalog': catdict, 
                'correction': corrdict
                }

        dlosclass = DlosEnv(i_cat_corr, n_NN=n_NN)
        dlosclass.read() 

        dlos_i = dlosclass.dlos
        env_i = dlosclass.env

        # Combine dLOS and env values 
        try: 
            combined_dlos += list(dlos_i)
            combined_env += list(env_i)
        except NameError: 
            combined_dlos = list(dlos_i)
            combined_env = list(env_i)

    comb_dlos = np.array(combined_dlos)
    comb_env = np.array(combined_env)

    print 'd'+str(n_NN)+'NN, Minimum ', comb_env.min(), ' Maximum ', comb_env.max()
    
    # Calculate fpeak and sigma for dLOS distribution in bins of 
    # galaxy environment. 
    istart = np.int( np.floor( comb_env.min()/np.float(dNN_stepsize) ) )
    iend = np.int( np.ceil( comb_env.max()/np.float(dNN_stepsize) ) )

    env_bins = [ 
            (dNN_stepsize * i, dNN_stepsize * (i+1))
            for i in np.arange(istart, iend) 
            ] 

    fpeaks, sigmas, amps, envbins, nbins = [], [], [], [], [] 

    for i_bin, env_bin in enumerate(env_bins): 
        
        in_envbin = np.where(
                (comb_env >= env_bin[0]) & 
                (comb_env < env_bin[1])
                )

        if len(in_envbin[0]) < 10: 
            continue 
    
        try: 
            fpeak, sigma, amp = dlospeak_fit(
                    comb_dlos[in_envbin], 
                    fit = 'gauss', 
                    peak_range = [-15.0, 15.0]
                    )   
        except ZeroDivisionError: 
            continue

        envbins.append( env_bin ) 
        fpeaks.append(fpeak)
        sigmas.append(sigma)
        amps.append(amp)
        nbins.append(np.float(len(in_envbin[0])))

    # calculate jack knife error estimates for fpeak and sigma by 
    # dividing the combined dLOS list into n_mock parts and calucating
    # fpeak(dNN) and env(dNN) for parts of combined dLOS not excluded
    
    jumbled_indices = np.arange(len(comb_dlos))
    np.random.shuffle( jumbled_indices ) 
    
    jackknife_fpeaks  = np.zeros(len(envbins))
    jackknife_sigmas  = np.zeros(len(envbins))
    n_jackknife_fpeak = np.zeros(len(envbins))
    n_jackknife_sigma = np.zeros(len(envbins))
    
    omit_binsize = int(
            np.floor( np.float(len(comb_dlos)) / np.float(n_mocks) )
            )

    for i_jack in range(n_mocks): 

        keep_indices = list( 
                set( range(len(comb_dlos)) ) -  \
                        set( range(i_jack * omit_binsize, (i_jack + 1) * omit_binsize) )
                        ) 
        keep_jumbled_indices = jumbled_indices[keep_indices]
        
        keep_comb_dlos = comb_dlos[keep_jumbled_indices]
        keep_comb_env  = comb_env[keep_jumbled_indices]
        
        for i_bin, env_bin in enumerate(envbins): 
        
            in_envbin = np.where(
                    (keep_comb_env >= env_bin[0]) & 
                    (keep_comb_env < env_bin[1])
                    )
            
            try: 
                fpeak, sigma, amp = dlospeak_fit(
                        keep_comb_dlos[in_envbin], 
                        fit = 'gauss', 
                        peak_range = [-15.0, 15.0]
                        )   

                if not math.isnan((fpeak - fpeaks[i_bin])**2): 
                    jackknife_fpeaks[i_bin] += (fpeak - fpeaks[i_bin])**2
                    n_jackknife_fpeak[i_bin] += 1.

                if not math.isnan((sigma - sigmas[i_bin])**2): 
                    jackknife_sigmas[i_bin] += (sigma - sigmas[i_bin])**2
                    n_jackknife_sigma[i_bin] += 1.

            except ZeroDivisionError: 
                pass 

    jackknife_fpeaks = np.sqrt( (n_jackknife_fpeak - 1.0)/n_jackknife_fpeak * jackknife_fpeaks ) 
    jackknife_sigmas = np.sqrt( (n_jackknife_sigma - 1.0)/n_jackknife_sigma * jackknife_sigmas ) 
    
    # write out the bestfit parameters for each of the environment bins
    if writeout: 

        envbins_low = np.array([ebin[0] for ebin in envbins])
        envbins_high = np.array([ebin[1] for ebin in envbins])
        
        data_list = [
                envbins_low, 
                envbins_high, 
                np.array(fpeaks), 
                jackknife_fpeaks, 
                np.array(sigmas), 
                jackknife_sigmas, 
                np.array(amps), 
                np.array(nbins)
                ]
        data_hdr = "Columns : dNN_low, dNN_high, fpeak, fpeak_err, sigma, sigma_err, amplitude, nbin"
        data_fmt = ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f']

        output_file = ''.join([
            (dlosclass.dlos_file).rsplit('/', 1)[0], '/', 
            'DLOS_env_d', str(n_NN), 'NN_bin_dist_bestfit.dat'
            ])

        np.savetxt(
                output_file, 
                (np.vstack(np.array(data_list))).T, 
                fmt = data_fmt, 
                delimiter = '\t', 
                header = data_hdr
                ) 
    
    return data_list 

def dlos_envbin_peakfit(cat_corr, n_NN=3, fit='gauss', **kwargs): 
    """ Best-fit parameters for the peak of the dLOS distribution
    classified into bins of galaxy environment. This function is 
    mostly a wrapper for the function catalog_dlospeak_env_fit.
    """
    catdict = cat_corr['catalog']

    dlosclass = Dlos(cat_corr)
    dlosclass.file_name
       
    # file that contains the bestfit parameters of the 
    # dLOS distribution for dLOS in bins of galaxy 
    # environment
    dlos_envbin_fit_file = ''.join([
        (dlosclass.file_name).rsplit('/', 1)[0], '/', 
        'DLOS_env_d', str(n_NN), 'NN_bin_dist_bestfit.dat'
        ])
    
    if not os.path.isfile(dlos_envbin_fit_file): 

        bestfit_list = catalog_dlospeak_env_fit(
                catdict['name'], 
                n_NN=n_NN,
                fit=fit, 
                writeout=True
                )

    env_low, env_high, fpeaks, fpeak_err, sigmas, sigma_err, amps, nbins = np.loadtxt(
            dlos_envbin_fit_file, 
            skiprows = 1, 
            unpack = True, 
            usecols = [0,1,2,3,4,5,6,7]
            )

    return [env_low, env_high, fpeaks, fpeak_err, sigmas, sigma_err, amps, nbins] 

def dlosenv_peakfit_fpeak_env_fit(cat_corr, n_NN=3, fit='gauss', **kwargs): 
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

    p0 = [ -0.01, 0.8 ] # guess
    fa = {'x': np.array( 0.5 * (env_low + env_high) ), 'y': fpeaks, 'err': fpeak_errs}
    
    fit_param = mpfit.mpfit(mpfit_linear, p0, functkw=fa)
        
    best_slope = fit_param.params[0]
    best_yint = fit_param.params[1]
    
    try: 
        if kwargs['writeout']: 
            dlos_fpeak_envbin_fit_file = ''.join([
                direc('data', cat_corr), 
                'DLOS_fpeak_env_d', str(n_NN), 'NN_bin_bestfit.dat'
                ])
        
            data_hdr = "Columns : bestfit_slope, bestfit_yint"
            data_fmt = ['%10.5f', '%10.5f']

            np.savetxt(
                    dlos_fpeak_envbin_fit_file, 
                    np.c_[best_slope, best_yint],
                    fmt = data_fmt, 
                    delimiter = '\t', 
                    header = data_hdr
                    ) 

    except KeyError: 
        pass

    return best_slope, best_yint

def dlosenv_peakfit_sigma_env_fit(cat_corr, n_NN=3, fit='gauss', **kwargs): 
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
        
    p0 = [ -0.03, 4.0 ] # guess
    fa = {'x': np.array( 0.5 * (env_low + env_high) ), 'y': sigmas, 'err': sigma_errs}
    
    fit_param = mpfit.mpfit(mpfit_linear, p0, functkw=fa)
        
    best_slope = fit_param.params[0]
    best_yint = fit_param.params[1]
    
    try: 
        if kwargs['writeout']: 
            dlos_fpeak_envbin_fit_file = ''.join([
                direc('data', cat_corr), 
                'DLOS_sigma_env_d', str(n_NN), 'NN_bin_bestfit.dat'
                ])
        
            data_hdr = "Columns : bestfit_slope, bestfit_yint"
            data_fmt = ['%10.5f', '%10.5f']

            np.savetxt(
                    dlos_fpeak_envbin_fit_file, 
                    np.c_[best_slope, best_yint],
                    fmt = data_fmt, 
                    delimiter = '\t', 
                    header = data_hdr
                    ) 

    except KeyError: 
        pass

    return best_slope, best_yint

#---- fit functions -----
def fit_linear(x, p): 
    ''' Linear function y = a * x + b 
    p[0] = a , p[1] = b
    '''
    return (p[0] * x) + p[1]

def peak_expon(x, p): 
    ''' Exponential function for the peak 
    '''
    return p[0]*np.exp(-np.abs(x)/p[1])

def peak_gauss(x, p): 
    """ Exponential function for the peak 
    """
    return p[0]*np.exp(-0.5*x**2/(p[1])**2)

# --- MPfit wrappers ---
def mpfit_linear(p, fjac=None, x=None, y=None, err=None): 
    model = fit_linear(x, p) 
    status = 0 
    if err == None: 
        err = np.array([1.0 for i in range(len(x))])
    return([status, (y-model)/err]) 

def mpfit_peak_expon(p, fjac=None, x=None, y=None): 
    model = peak_expon(x, p) 
    status = 0 
    return([status, (y-model)]) 

def mpfit_peak_gauss(p, fjac=None, x=None, y=None): 
    model = peak_gauss(x, p) 
    status = 0 
    return([status, (y-model)]) 

if __name__=="__main__":
    cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock':1}, 
            'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 4.0, 'fpeak':0.69}
            }
    
    for nNN in [1,3,5,7,10]: 
        fpeak_slope, fpeak_yint = dlos_peakfit_fpeak_env_fit(
                cat_corr, 
                n_NN=nNN, 
                fit='gauss', 
                writeout=True
                )
        sigma_slope, sigma_yint = dlos_peakfit_sigma_env_fit(
                cat_corr, 
                n_NN=nNN, 
                fit='gauss', 
                writeout=True
                )
