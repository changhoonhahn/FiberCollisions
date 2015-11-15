'''

Fiber collision Correction using the Fourier transformed tophat correction method 

'''
import time
import pickle
import numpy as np
from scipy.special import legendre
from scipy.special import j0 as J0
from scipy.special import j1 as J1
from scipy.special import jv
from scipy.integrate import quad as quad_int 
from scipy.interpolate import interp1d

from util.interruptible_pool import InterruptiblePool as Pewl 

def delP_uncorr(k, l, fs=1.0, rc=0.4):
    '''
    Fiber collision correction term to Powerspectrum multipole l. The uncorrelated term.

    delP_uncorr(k)_l = - (2 pi f_s) (pi rc^2) (2 l + 1)/2 Leg_l(0) W_2d(k*rc)/k

    '''
    leg_l = legendre(l)

    alpha = -1.0 * np.pi**2 * rc**2 * fs 

    delP = alpha * (2.*np.float(l)+1) * leg_l(0) * W_2d(k*rc)/k

    return delP

def delP_corr(k, Pk, l, fs=1.0, rc=0.4, extrap_params=[[3345.0, -1.6], [400.0, -4.]]): 
    '''
    Fiber collision correction term to Powerspectrum multipole l. The correlated term.
    
    delP_corr(k)_l = - (2 pi f_s) (2 pi^2 rc^2) Sum_l' [int^inf_0 q Pl'(q) f_ll'(q*rc,k*rc) dq]
    
    Leg_l(0) W_2d(k*rc)/k

    extrap_params : [[alpha_0, n_0], [alpha_2, n_2], [alpha_4, n_4]]

    '''

    alpha = -4.0 * np.pi**3 * fs * rc**2

    if not isinstance(Pk, list): 
        raise ValueError('Pk input has to be a list containing the multipoles. e.g. [P0k, P2k, P4k]')

    if len(Pk) != len(extrap_params): 
        raise ValueError('Extrapolation parameters need to be specified for each multipole') 
    
    lp_max = len(Pk)

    sum_delP = np.zeros(len(k))
    for i_lp, lp in enumerate(np.arange(0.0, lp_max+0.1, 2)):
        
        Pk_interp = interp1d(k, Pk[i_lp], kind='cubic')
        Pk_kplus_extrap = lambda qork: pk_powerlaw(qork, extrap_params[i_lp])

        delP = np.zeros(len(k))
        for i_k, k_i in enumerate(k): 
            print 'l=', l, " l'=", int(lp), ', k = ', k_i
            start_time = time.time()
            integrand = lambda q_var: delPq_integrand(
                    q_var, k_i, l, int(lp), Pk_interp, Pk_kplus_extrap, rc=rc, k_min=k[0], k_max=k[-1]
                    )
            delP[i_k] = quad_int(integrand, 0., np.inf, epsrel=0.01)[0]
            print "delP(k_i) calculation takes ", time.time() - start_time 

        sum_delP += delP

    sum_delP *= alpha
    
    return sum_delP 

def delPq_integrand(q, k, l, lp, f_interp, f_extrap, rc=0.4, k_min=0.002, k_max=0.3): 
    '''
    Integrand for del P_l^corr correction 

    q * dq * P_lp(q) * f_llp(q*rc, k_i*rc)
    '''
    if (q > k_min) and (q < k_max): 
        integrand = q * f_interp(q) * f_l_lp(q*rc, k*rc, l, lp)
    elif q <= k_min: 
        integrand = 0.0
    elif q >= k_max: 
        integrand = q * f_extrap(q) * f_l_lp(q*rc, k*rc, l, lp)
    else: 
        raise ValueError

    return integrand

def f_l_lp(x, y, l, lp, error=False): 
    '''
    f_l,l' term in the FC correction equation 

    f_l,l'(x,y) = (2*l+1)/2 * int_-1^min(x/y, 1) du L_l(u) L_l'(y/x*u) * (W_2d(theta) + y^2(1-u^2) (J2(theta)/theta^2 - 1/4 W_2d(theta))

    where theta = sqrt(x^2 - y^2*u^2)
    '''
    #start_time = time.time()
    upper_bound = np.min([1.0, x/y])
    lower_bound = np.max([-1.0, -1.0*x/y])
    #print 'upper bound ', upper_bound, ' lower bound ', lower_bound
    
    fllp = lambda mu : f_l_lp_integrand(mu, x, y, l, lp)

    f_llp = 0.5*(2.*np.float(l) + 1.) * quad_int(fllp, lower_bound,  upper_bound, epsrel=0.01)[0]
    f_llp_err = 0.5*(2.*np.float(l) + 1.) * quad_int(fllp, lower_bound, upper_bound, epsrel=0.01)[1]
    
    #print "f_ll' calculation takes ", time.time() - start_time 
    if not error: 
        return f_llp
    else: 
        return [f_llp, f_llp_err]

def f_l_lp_integrand(mu, x, y, l, lp): 
    '''
    Integrand of f_l_lp integration 
    '''
    Leg_l = legendre(l)     # Legendre polynomial of order l 
    Leg_lp = legendre(lp)   # Legendre polynomial of order l'
    
    if np.abs(mu) > np.abs(x/y): 
        raise ValueError

    theta = np.sqrt(x**2 - y**2 * mu**2)

    return Leg_l(mu) * Leg_lp(y * mu/x) * (W_2d(theta) + y**2 * (1. - mu**2) * W_secorder(theta))

def pk_powerlaw(k, param): 
    return param[0] * (k/0.3)**param[1]

def W_secorder(x): 
    '''
    J2(x)/x^2 - J1(x)/2x
    '''
    return J2(x)/x**2 - J1(x)/(2*x)

def W_2d(x): 
    '''
    W_2d in FC correction equation 

    W_2d = (2 J_1(x)/ x)
    '''
    return 2.*J1(x)/x

def J2(x): 
    '''
    Spherical Bessel Function of the first kind of order 2 
    '''
    return jv(2, x)

if __name__=='__main__': 
    k, P0k, P2k, P4k = np.loadtxt('/mount/riachuelo1/hahn/power/Nseries/POWER_Q_CutskyN1.fidcosmo.dat.grid360.P020000.box3600', 
            unpack=True, usecols=[0,1,2,4])

    Pk = [P0k, P2k]#, P4k]
    #for l1 in [0, 2]:
    #    for l2 in [0,2,4]:
    #        for k_i in [1e-8, 1e-3, 0.002, 0.3, 1.0, 2.0, 10.0, 100.0]: 
    #            print l1, l2, k, f_l_lp(k_i, 0.3, l1, l2)

    #print delP_uncorr(k, 0, fs=1.0, rc=0.4)
    #print delP_uncorr(k, 2, fs=1.0, rc=0.4)
    delp0 = delP_corr(k, Pk, 0, fs=1.0, rc=0.4, extrap_params=[[3345.0, -1.6], [400.0, -4.]])
    pickle.dump( delp0, open('testing/delP0k_corr.p', 'wb'))
    delp2 = delP_corr(k, Pk, 2, fs=1.0, rc=0.4, extrap_params=[[3345.0, -1.6], [400.0, -4.]])
    pickle.dump( delp2, open('testing/delP2k_corr.p', 'wb'))
    #print delP_corr(k, Pk, 2, fs=1.0, rc=0.4)
