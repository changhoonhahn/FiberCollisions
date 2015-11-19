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

def delP_uncorr(k, l, fs=1.0, rc=0.4):
    '''
    Fiber collision correction term to Powerspectrum multipole l. The uncorrelated term.

    delP_uncorr(k)_l = - (2 pi f_s) (pi rc^2) (2 l + 1)/2 Leg_l(0) W_2d(k*rc)/k

    '''
    leg_l = legendre(l)

    alpha = -1.0 * np.pi**2 * rc**2 * fs 

    delP = alpha * (2.*np.float(l)+1) * leg_l(0) * W_2d(k*rc)/k

    return delP

def delP_corr(k, Pk, l, fs=1.0, rc=0.4, 
        extrap_params=[[3345.0, -1.6], [400.0, -4.]], 
        first_order=False, 
        w1d=False, 
        lpcomponent=False, 
        noextrap=False): 
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
    
    lp_max = len(Pk)

    sum_delP = np.zeros(len(k))
    lps, delPlp = [], [] 
    print np.arange(0.0, 2.0 * np.float(lp_max) - 0.1, 2)
    for i_lp, lp in enumerate(np.arange(0.0, 2.0 * np.float(lp_max)-0.1, 2)):
        
        Pk_interp = interp1d(k, Pk[i_lp], kind='cubic')
        #print extrap_params[i_lp]
        Pk_kplus_extrap = lambda q_or_k: pk_powerlaw(q_or_k, extrap_params[i_lp])

        delP = np.zeros(len(k))
        start_time = time.time()
        for i_k, k_i in enumerate(k): 
            dPq_integrand = lambda q_var: delPq_integrand(
                    q_var, 
                    k_i, 
                    l, 
                    int(lp), 
                    Pk_interp, 
                    Pk_kplus_extrap, 
                    rc=rc, 
                    k_min=k[0], 
                    k_max=k[-1], 
                    first_order=first_order,
                    w1d=w1d, 
                    noextrap=noextrap
                    )
            delP_int = quad_int(dPq_integrand, 0., np.inf, epsrel=0.01)
            delP[i_k] = alpha * delP_int[0]

        print "delP(k) calculation takes ", time.time() - start_time 
        print 'l=', l, " l'=", int(lp)
        print delP 
        sum_delP += delP
    
        if lpcomponent: 
            lps.append(lp)
            delPlp.append(delP)

    if lpcomponent:
        return [lps, delPlp]
    else: 
        return sum_delP 

def delPq_integrand(q, k, l, lp, f_interp, f_extrap, rc=0.4, k_min=0.002, k_max=0.3, first_order=False, w1d=False, noextrap=False):
    '''
    Integrand for del P_l^corr correction 

    q * dq * P_lp(q) * f_llp(q*rc, k_i*rc)
    '''
    if (q > k_min) and (q < k_max): 
        #integrand = q * f_interp(q) * f_l_lp(q*rc, k*rc, l, lp, first_order=first_order, w1d=w1d)
        integrand = q * f_interp(q) * f_l_lp_est(q*rc, k*rc, l, lp)
    elif q <= k_min: 
        integrand = 0.0
    elif q >= k_max: 
        #integrand = q * f_extrap(q) * f_l_lp(q*rc, k*rc, l, lp, first_order=first_order, w1d=w1d)
        if not noextrap: 
            integrand = q * f_extrap(q) * f_l_lp_est(q*rc, k*rc, l, lp)
        else: 
            integrand = 0.0
        print 'q = ', q, 'integrand = ', integrand
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
    
    if (l1_in == 2) and (l2_in == 0): 
        return x**2 - 1.
    elif (l1_in == 4) and (l2_in == 0): 
        return (7./4.)*x**4 - (5./2.)*x**2 + 3./4.
    elif (l1_in == 4) and (l2_in == 2): 
        return x**4 - x**2
    else: 
        raise ValueError

def f_l_lp(x, y, l, lp, error=False, first_order=False, w1d=False): 
    '''
    f_l,l' term in the FC correction equation 

    f_l,l'(x,y) = (2*l+1)/2 * int_-1^min(x/y, 1) du L_l(u) L_l'(y/x*u) * (W_2d(theta) + y^2(1-u^2) (J2(theta)/theta^2 - 1/4 W_2d(theta))

    where theta = sqrt(x^2 - y^2*u^2)
    '''
    upper_bound = np.min([1.0, x/y])
    lower_bound = np.max([-1.0, -1.0*x/y])
    #print 'upper bound ', upper_bound, ' lower bound ', lower_bound
    
    fllp = lambda mu : f_l_lp_integrand(mu, x, y, l, lp, first_order=first_order, w1d=w1d)

    int_fllp = quad_int(fllp, lower_bound,  upper_bound, epsrel=0.01)
    f_llp = (2.*np.float(l) + 1.)/2. * int_fllp[0]
    
    if not error: 
        return f_llp
    else: 
        f_llp_err = 0.5*(2.*np.float(l) + 1.) * int_fllp[1]
        return [f_llp, f_llp_err]

def f_l_lp_integrand(mu, x, y, l, lp, first_order=False, w1d=False): 
    '''
    Integrand of f_l_lp integration 
    '''
    Leg_l = legendre(l)     # Legendre polynomial of order l 
    Leg_lp = legendre(lp)   # Legendre polynomial of order l'
    
    if np.abs(mu) > np.abs(x/y): 
        raise ValueError

    theta = np.sqrt(x**2 - y**2 * mu**2)
    
    if w1d: 
        output = Leg_l(mu) * Leg_lp(y * mu/x) * W_1d(y * np.sqrt(1. - mu**2), np.sqrt(x**2 - y**2 * mu**2))
    else: 
        if first_order: 
            output = Leg_l(mu) * Leg_lp(y * mu/x) * W_2d(theta)
        else: 
            output = Leg_l(mu) * Leg_lp(y * mu/x) * (W_2d(theta) + y**2 * (1. - mu**2) * W_secorder(theta))

    return output 

def pk_powerlaw(k, param): 
    return param[0] * (k/0.3)**param[1]

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

def delPcorr_estimated(lpcomponent=True, noextrap=False):
    n_mocks = 20
    for i in xrange(1,n_mocks+1): 
        k_i, P0k_i, P2k_i, P4k_i = np.loadtxt(
                ''.join(['/mount/riachuelo1/hahn/power/Nseries/', 
                    'POWER_Q_CutskyN1.fidcosmo.dat.grid360.P020000.box3600']), 
                unpack=True, 
                usecols=[0,1,2,4]
                )
        
        if i == 1: 
            k = k_i 
            P0k_sum = P0k_i
            P2k_sum = P2k_i
            P4k_sum = P4k_i
        else: 
            P0k_sum += P0k_i
            P2k_sum += P2k_i
            P4k_sum += P4k_i

    Pk = [P0k_sum/np.float(n_mocks), P2k_sum/np.float(n_mocks), P4k_sum/np.float(n_mocks)]
    
    for l in [0, 2, 4]: 
        delp = delP_corr(
                k, Pk, l, fs=1.0, rc=0.4, 
                extrap_params=[[3345.0, -1.6], [400.0, -4.], [260.0, -1.0]], 
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
            
        pickle.dump( delp, open('delP'+str(l)+'k_corr_estimated'+lpcomponent_str+noextrap_str+'.p', 'wb'))

if __name__=='__main__': 
    delPcorr_estimated(lpcomponent=True, noextrap=True)
    delPcorr_estimated(lpcomponent=False, noextrap=True)
    delPcorr_estimated(lpcomponent=True, noextrap=False)
    delPcorr_estimated(lpcomponent=False, noextrap=False)
