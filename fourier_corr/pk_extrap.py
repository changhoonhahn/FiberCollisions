'''

Extrapolation of P_l(k) used in calculation delP_fc^corr

'''
import pickle
import numpy as np
# --- Local --- 
import mpfit

def pk_powerlaw_bestfit(k_in, pk_in, k_fit=0.25, k_fixed=0.6, quiet=True): 
    '''
    Find best fit power law for P_l(k) extrapolatin
    '''

    fit_range = np.where(k_in > k_fit)
    if len(fit_range[0]) == 0: 
        raise ValueError('k_fit range does not make sense')

    fa = {'x': k_in[fit_range], 'y': pk_in[fit_range], 'k_fixed': k_fixed}

    param_guess = [pk_in[-1], -2.]

    bestfit = mpfit.mpfit(pk_powerlaw_mpfit, param_guess, functkw=fa, quiet=quiet)

    alpha = bestfit.params[0]
    n = bestfit.params[1]

    return [alpha, n]

def pk_powerlaw(k, param, k_fixed=0.6): 
    '''
    P(k) = alpha * (k/k_fixed)^n

    param = [alpha, n]
    '''
    return param[0] * (k/k_fixed)**param[1]

def pk_powerlaw_mpfit(param, fjac=None, x=None, y=None, k_fixed=None): 
    '''
    mpfit wrapper for pk_powerlaw
    '''
    model = pk_powerlaw(x, param, k_fixed=k_fixed)
    status = 0 

    return([status, y-model])
