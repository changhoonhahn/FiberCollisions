'''

Fourier transform of 2PCF top hat Fiber collision Correction 

'''
import time
import pickle
import numpy as np

# --- Polynomials ---
from scipy.special import legendre
from scipy.special import j0 as J0
from scipy.special import j1 as J1
from scipy.special import jv
from scipy.integrate import quad as quad_int 
from scipy.interpolate import interp1d

import pk_extrap 

def delP_uncorr(k, l, fs=1.0, rc=0.4):
    '''
    Fiber collision correction term to Powerspectrum multipole l. The uncorrelated term.

    delP_uncorr(k)_l = - (2 pi f_s) (pi rc^2) (2 l + 1)/2 Leg_l(0) W_2d(k*rc)/k

    '''
    leg_l = legendre(l)

    alpha = -1.0 * np.pi**2 * rc**2 * fs 

    delP = alpha * (2.*np.float(l)+1) * leg_l(0) * W_2d(k*rc)/k

    return delP

def delP_corr(k, Pk, l, fs=1.0, rc=0.4, extrap_params=[[3345.0, -1.6], [400.0, -4.]], k_fixed=None): 
    '''
    Fiber collision correction term to Powerspectrum multipole l. The correlated term.
    
    delP_corr(k)_l = - (2 pi f_s) (2 pi^2 rc^2) Sum_l' [int^inf_0 q Pl'(q) f_ll'(q*rc,k*rc) dq]
    
    Leg_l(0) W_2d(k*rc)/k

    extrap_params : [[alpha_0, n_0], [alpha_2, n_2], [alpha_4, n_4]]

    '''

    alpha = -0.5 * fs * rc**2

    if not isinstance(Pk, list): 
        raise ValueError('Pk input has to be a list containing the multipoles. e.g. [P0k, P2k, P4k]')

    if len(Pk) != len(extrap_params): 
        raise ValueError('Extrapolation parameters need to be specified for each multipole') 
    
    lvalues = [0,2,4,6,8]

    sum_delP = np.zeros(len(k))

    for i_lp, lp in enumerate(lvalues[:len(Pk)]):
        
        # Cubic Spline interpolated function of P(k)
        Pk_interp = interp1d(k, Pk[i_lp], kind='cubic')

        # Extrapolated power law for P(k) k > k_max
        Pk_kplus_extrap = lambda q_or_k: \
                pk_extrap.pk_powerlaw(
                        q_or_k, 
                        extrap_params[i_lp], 
                        k_fixed=k_fixed
                        )

        delP = np.zeros(len(k))

        # loop through k values and evaluate delta P for all ks 
        for i_k, k_i in enumerate(k): 

            dPq_integrand = lambda q_var: delPq_integrand(
                    q_var, 
                    k_i, 
                    l, 
                    lp, 
                    Pk_interp, 
                    Pk_kplus_extrap, 
                    rc=rc, 
                    k_min=k[0], 
                    k_max=k[-1]
                    )
            delP_int, delP_int_err = quad_int(dPq_integrand, 0., np.inf, epsrel=0.01)
            delP[i_k] = alpha * delP_int

        #print "delP(k) calculation takes ", time.time() - start_time 
        #print 'l=', l, " l'=", int(lp)
        #print delP 
        sum_delP += delP
    
    return sum_delP 

def delPq_integrand(q, k, l, lp, f_interp, f_extrap, rc=0.4, k_min=0.002, k_max=0.3, noextrap=False):
    '''
    Integrand for del P_l^corr correction 

    q * dq * P_lp(q) * f_llp(q*rc, k_i*rc)
    '''
    if (q > k_min) and (q < k_max): 
        integrand = q * f_interp(q) * f_l_lp_est(q*rc, k*rc, l, lp)
    elif q <= k_min: 
        integrand = 0.0
    elif q >= k_max: 
        if not noextrap: 
            integrand = q * f_extrap(q) * f_l_lp_est(q*rc, k*rc, l, lp)
        else: 
            integrand = 0.0
        #print 'q = ', q, 'integrand = ', integrand
    else: 
        raise ValueError

    return integrand

def f_l_lp_est(x, y, l, lp): 
    '''
    f_l,l' term in the FC correction equation. 

    f_l,l'(x,y) = (2*l+1)/2 * int_-1^min(x/y, 1) du L_l(u) L_l'(y/x*u) * (W_2d(theta) + y^2(1-u^2) (J2(theta)/theta^2 - 1/4 W_2d(theta))

    where theta = sqrt(x^2 - y^2*u^2). In this convention, x = q * rc and y = q * rc. 

    Below are estimates of the integral, without actually computing the integral 
    '''
    
    if l == lp:     # l = l' (diagonal terms)

        if x < y: 
            return W_2d(x) * (x / y)**(l+1)
        else: 
            return W_2d(x) * (y / x)**l

    elif l > lp: 
        if x < y: 
            return (2.0 * l + 1.)/2.0 * (x / y) * fllp_poly(l, lp, x / y) * W_2d(x)
        else: 
            return 0.0

    elif l < lp: 
        if x > y: 
            return (2.0 * l + 1.)/2.0 * fllp_poly(lp, l, y / x) * W_2d(x)
        else: 
            return 0.0 

def fllp_poly(l1_in, l2_in, x):  
    '''
    Polygons for specific situations called by f_l_lp_est function
    '''

    if (l1_in == 2) and (l2_in == 0): 
        return x**2 - 1.
    elif (l1_in == 4) and (l2_in == 0): 
        return (7./4.)*x**4 - (5./2.)*x**2 + 3./4.
    elif (l1_in == 4) and (l2_in == 2): 
        return x**4 - x**2
    else: 
        raise ValueError

def W_secorder(x): 
    '''
    J2(x)/x^2 - J1(x)/2x
    '''
    return J2(x)/(x**2) - J1(x)/(2*x)

def W_2d(x): 
    '''
    W_2d in FC correction equation 

    W_2d = (2 J_1(x)/ x)
    '''
    return 2.*J1(x)/x

def W_1d(x,y): 
    '''
    W_1d(x,y) = int^2pi_0 W_2D(sqrt(x^2 + y^2 - 2xy cos(phi))) dphi/2pi
    '''
    integrand = lambda phi: W_2d( np.sqrt(x**2 + y**2 - 2*x*y*np.cos(phi)) )

    output = 1./(2.*np.pi) * quad_int(integrand, 0., 2.0 * np.pi, epsrel=0.01)[0]
    
    return output

def J2(x): 
    '''
    Spherical Bessel Function of the first kind of order 2 
    '''
    return jv(2, x)

# Specific cases of delP_corr  ----
def delP_corr_noextrap(k, Pk, l, fs=1.0, rc=0.4): 
    '''
    The correlated term of the fiber collision correction to Powerspectrum multipole l but with NO P_l(k) extrapolation.
    In other words P_l'(k) = 0 for k > k_max. This is most likely not accurate, but used for testing.
    '''

    alpha = -0.5 * fs * rc**2

    if not isinstance(Pk, list): 
        raise ValueError(
                'Pk input has to be a list containing the multipoles. e.g. [P0k, P2k, P4k]'
                )
    
    lvalues = [0,2,4,6,8]

    sum_delP = np.zeros(len(k))

    lps, delPlp = [], [] 

    for i_lp, lp in enumerate(lvalues[:len(Pk)]):
        
        # Cubic Spline interpolated function of P(k)
        Pk_interp = interp1d(k, Pk[i_lp], kind='cubic')
        # Extrapolated power law for P(k) k > k_max
        Pk_kplus_extrap = lambda q_or_k: 0.0 * q_or_k

        delP = np.zeros(len(k))

        start_time = time.time()
        
        # loop through k values and evaluate delta P for all ks 
        for i_k, k_i in enumerate(k): 

            dPq_integrand = lambda q_var: delPq_integrand(
                    q_var, 
                    k_i, 
                    l, 
                    lp, 
                    Pk_interp, 
                    Pk_kplus_extrap, 
                    rc=rc, 
                    k_min=k[0], 
                    k_max=k[-1], 
                    noextrap=True
                    )
            delP_int, delP_int_err = quad_int(dPq_integrand, 0., np.inf, epsrel=0.01)
            delP[i_k] = alpha * delP_int

        sum_delP += delP
    return sum_delP 

def delP_corr_lp_comp(k, Pk, l, fs=1.0, rc=0.4, extrap_params=[[3345.0, -1.6], [400.0, -4.]], k_fixed=None): 
    '''
    Each of the l' components of the correlated term of the fiber collision correction to Powerspectrum multipole l.
    '''

    alpha = -0.5 * fs * rc**2

    if not isinstance(Pk, list): 
        raise ValueError(
                'Pk input has to be a list containing the multipoles. e.g. [P0k, P2k, P4k]'
                )

    if len(Pk) != len(extrap_params): 
        raise ValueError('Extrapolation parameters need to be specified for each multipole') 
    
    lvalues = [0,2,4,6,8]

    sum_delP = np.zeros(len(k))

    lps, delPlp = [], [] 

    for i_lp, lp in enumerate(lvalues[:len(Pk)]):
        
        # Cubic Spline interpolated function of P(k)
        Pk_interp = interp1d(k, Pk[i_lp], kind='cubic')

        # Extrapolated power law for P(k) k > k_max
        Pk_kplus_extrap = lambda q_or_k: \
                pk_extrap.pk_powerlaw(
                        q_or_k, 
                        extrap_params[i_lp], 
                        k_fixed=k_fixed
                        )

        delP = np.zeros(len(k))

        start_time = time.time()
        
        # loop through k values and evaluate delta P for all ks 
        for i_k, k_i in enumerate(k): 

            dPq_integrand = lambda q_var: delPq_integrand(
                    q_var, 
                    k_i, 
                    l, 
                    lp, 
                    Pk_interp, 
                    Pk_kplus_extrap, 
                    rc=rc, 
                    k_min=k[0], 
                    k_max=k[-1]
                    )
            delP_int, delP_int_err = quad_int(dPq_integrand, 0., np.inf, epsrel=0.01)
            delP[i_k] = alpha * delP_int

        #print "delP(k) calculation takes ", time.time() - start_time 
        #print 'l=', l, " l'=", int(lp)
        #print delP 
        sum_delP += delP
    
        lps.append(lp)
        delPlp.append(delP)

    return [lps, delPlp]

# save to pickle files  ----
def delPcorr_estimated(n_mocks=10, k_fixed=0.6, lpcomponent=True, noextrap=False, Ngrid=360):
    '''
    Wrapper function to calculate delta P^corr
    '''
    
    # average P_l(k)
    k0, avg_P0k = pk_extrap.average_Pk(0, n_mocks=n_mocks, Ngrid=Ngrid)
    k2, avg_P2k = pk_extrap.average_Pk(2, n_mocks=n_mocks, Ngrid=Ngrid)
    k4, avg_P4k = pk_extrap.average_Pk(4, n_mocks=n_mocks, Ngrid=Ngrid)
    k = k0

    Pk = [avg_P0k, avg_P2k, avg_P4k]
    
    for l in [0, 2, 4]: 

        delp = delP_corr(
                k, Pk, l, fs=1.0, rc=0.4, 
                extrap_params=[[1009., -2.1], [-524., -0.251], [443., 0.72]],
                k_fixed=k_fixed,
                lpcomponent=lpcomponent, 
                noextrap=noextrap
                )

        if lpcomponent: 
            lpcomponent_str = '_lpcomponent'
        else: 
            lpcomponent_str = ''

        if noextrap: 
            noextrap_str = '_noextrap'
        else: 
            noextrap_str = ''
        
        pickle_file = ''.join([
            'delP', str(l), 'k_corr_estimated', 
            lpcomponent_str, noextrap_str, 
            '_k_fixed', str(round(k_fixed, 1)), 
            '_Ngrid', str(Ngrid), '.p'
            ])
        pickle.dump(delp, open(pickle_file, 'wb'))

def delPcorr_extrapolation(n_mocks=10, k_fixed=0.6, k_max=0.5, Ngrid=360):
    '''
    Wrapper function to calculate delta P^corr
    '''
    # average P_l(k)
    k0, avg_P0k = pk_extrap.average_Pk(0, n_mocks=n_mocks, Ngrid=Ngrid)
    k2, avg_P2k = pk_extrap.average_Pk(2, n_mocks=n_mocks, Ngrid=Ngrid)
    k4, avg_P4k = pk_extrap.average_Pk(4, n_mocks=n_mocks, Ngrid=Ngrid)

    Pk = [avg_P0k, avg_P2k, avg_P4k]
    k = k0

    if isinstance(k_max, list) or isinstance(k_max, np.ndarray):  
        pass
    else: 
        k_max = [k_max]

    for k_max_i in k_max: 

        bestfit_params = [] 
        for i_l, l in enumerate([0, 2, 4]):

            bestfit_params.append(
                    pk_extrap.pk_bestfit(k0, Pk[i_l], k_max=k_max_i, k_fixed=k_fixed)
                    )
        print k_max_i
        print bestfit_params 
    
        for i_l, l in enumerate([0, 2, 4]): 
            
            delp = delP_corr(
                    k, Pk, l, 
                    fs=1.0, 
                    rc=0.4, 
                    extrap_params=bestfit_params,
                    k_fixed=k_fixed
                    )

            pickle_file = ''.join([
                'delP', str(l), 'k_corr', 
                '_k_fixed', str(round(k_fixed, 1)),
                '_kmax', str(round(k_max_i, 2)), 
                '_Ngrid', str(Ngrid), '.p'
                ])
            pickle.dump(delp, open(pickle_file, 'wb'))
