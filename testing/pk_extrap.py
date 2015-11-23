'''

Extrapolate P_l(k)

'''
import pickle
import numpy as np

import mpfit

def average_Pk(l, n_mocks, Ngrid=360, quiet=True): 
    '''
    Return averaged Nseries P_l(k) from data
    '''
    if l == 0: 
        data_cols = [0, 1]
    elif l == 2: 
        data_cols = [0, 2]
    elif l == 4: 
        data_cols = [0, 4]
    else: 
        raise ValueError
    
    n_files = 0.
    for i_mock in xrange(1, n_mocks+1): 

        pk_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/', 
            'POWER_Q_CutskyN', str(i_mock), '.fidcosmo.dat.grid', str(Ngrid), 
            '.P020000.box3600'
            ])
        try: 
            k_i, Plk_i = np.loadtxt(
                    pk_file, 
                    unpack = True, 
                    usecols = data_cols
                    )
            if quiet: 
                print pk_file

            n_files += 1.

        except IOError: 

            continue
    
        try: 
            Plk_sum += Plk_i
        except UnboundLocalError: 
            k = k_i 
            Plk_sum = Plk_i

    print int(n_files), ' P(k) files have beeen read'
    Plk = Plk_sum/n_files

    return [k, Plk]

def average_Pk_upw(l, n_mocks, Ngrid=360, quiet=True): 
    '''
    Return averaged Nseries P_l(k) from data
    '''
    if l == 0: 
        data_cols = [0, 1]
    elif l == 2: 
        data_cols = [0, 2]
    elif l == 4: 
        data_cols = [0, 4]
    else: 
        raise ValueError
    
    n_files = 0.
    for i_mock in xrange(1, n_mocks+1): 

        pk_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/', 
            'POWER_Q_CutskyN', str(i_mock), '.fidcosmo.fibcoll.dat.grid', str(Ngrid), 
            '.P020000.box3600'
            ])
        try: 
            k_i, Plk_i = np.loadtxt(
                    pk_file, 
                    unpack = True, 
                    usecols = data_cols
                    )
            if quiet: 
                print pk_file

            n_files += 1.

        except IOError: 

            continue
        
        try: 
            Plk_sum += Plk_i
        except UnboundLocalError: 
            k = k_i 
            Plk_sum = Plk_i

    print int(n_files), ' P(k) files have beeen read'
    Plk = Plk_sum/n_files

    return [k, Plk]

def pk_powerlaw(k, param, k_fixed=0.6): 
    return param[0] * (k/k_fixed)**param[1]

def pk_powerlaw_mpfit(param, fjac=None, x=None, y=None, k_fixed=None): 
    '''
    P(k) = alpha * (k/k_fixed)^n

    param = [alpha, n]
    '''
    model = pk_powerlaw(x, param, k_fixed=k_fixed)
    status = 0 

    return([status, y-model])

def pk_bestfit(k_in, pk_in, k_max=0.25, k_fixed=0.6, quiet=True): 
    '''
    Find best fit power law for P_l(k) extrapolatin
    '''

    fit_range = np.where(k_in > k_max)

    fa = {'x': k_in[fit_range], 'y': pk_in[fit_range], 'k_fixed': k_fixed}

    param_guess = [pk_in[-1], -2.]

    bestfit = mpfit.mpfit(pk_powerlaw_mpfit, param_guess, functkw=fa, quiet=quiet)

    alpha = bestfit.params[0]
    n = bestfit.params[1]

    return [alpha, n]

if __name__=="__main__": 
    for l_i in [0, 2, 4]: 
        #for k_max_i in np.arange(0.4, 0.65, 0.05): 
        k_max_i = 0.6
        print 'l = ', l_i, 'k_max = ', k_max_i
        k_i, Plk_i = average_Pk(l_i, 10, Ngrid=720)
        print pk_bestfit(k_i, Plk_i, k_max=k_max_i, k_fixed=0.6)
