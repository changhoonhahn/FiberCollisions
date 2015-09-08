"""

Fitting the peak of the line-of-sight displacement distribution


"""
import numpy as np

import mpfit
from dlos import Dlos

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
