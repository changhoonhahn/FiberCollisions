
import time
import pickle
import numpy as np
from scipy.special import legendre
from scipy.special import j0 as J0
from scipy.special import j1 as J1
from scipy.special import jv
from scipy.integrate import quad as quad_int 
from scipy.interpolate import interp1d

def delP_corr(l, k_data, P_data, extrap_params=[[3345.0, -1.6], [400.0, -4.]], rc = 0.4, fs=1.0): 
    
    alpha = -4.0 * np.pi**3 * fs * rc**2
    
    delP_corr = np.zeros(len(k_data))
    for i_k, k_i in enumerate(k_data):
        if i_k % 18 == 17:  
            start_time = time.time()
            delP_int = outside_integral(k_i, l, k_data, P_data, extrap_params=extrap_params, rc=rc)
            delP_corr[i_k] = delP_int[0]
            print 'P_'+str(l)+'('+str(k_i)+') took ', time.time() - start_time

    delP_corr *=  alpha * (2.0 * np.float(l) + 1.0)/2.0

    return delP_corr

def outside_integral(k, l, k_data, P_data, extrap_params=[[3345.0, -1.6], [400.0, -4.]], rc = 0.4): 
    
    integrand = lambda mu: outside_integrand(mu, k, l, k_data, P_data, extrap_params=extrap_params, rc = 0.4)
    
    outside_int = quad_int(integrand, -1., 1., epsrel=0.01)

    return outside_int

def outside_integrand(mu, k, l, k_data, P_data, extrap_params=[[3345.0, -1.6], [400.0, -4.]], rc = 0.4): 

    Leg_l = legendre(l)

    integrand = Leg_l(mu) * inside_integral(mu, k, k_data, P_data, extrap_params=extrap_params, lower_bound=np.abs(k * mu), rc=rc)

    return integrand

def inside_integral(mu, k, k_data, P_data, extrap_params=[[3345.0, -1.6], [400.0, -4.]], lower_bound = None, rc=0.4): 
    
    if len(P_data) != len(extrap_params): 
        raise ValueError
    
    inside_integral = 0.0
    lps = [0, 2, 4, 6, 8]
    for i_lp, P_lp in enumerate(P_data): 
        #print 'l'=', lps[i_lp]
        Pk_interp = interp1d(k_data, P_lp, kind = 'cubic')
        Pk_extrap = lambda q_or_k: pk_powerlaw(q_or_k, extrap_params[i_lp])
    
        integrand = lambda q: inside_integrand(mu, q, k, lps[i_lp], Pk_interp, Pk_extrap, rc=rc, q_min=k_data[0], q_max=k_data[-1])

        inside_int = quad_int(integrand, lower_bound, np.inf, epsrel=0.01)

        inside_integral += inside_int[0]
    
    return inside_integral


def inside_integrand(mu, q, k, lp, f_interp, f_extrap, rc=0.4, q_min=0.002, q_max=0.3): 
    
    Leg_lp = legendre(lp)   # Legendre polynomial of order l'

    if q < q_min:
        return 0.0
    elif q > q_max: 
        return f_extrap(q) * Leg_lp(k*mu/q) * W_2d(rc * np.sqrt(q**2 - k**2 * mu**2)) * q
    else: 
        return f_interp(q) * Leg_lp(k*mu/q) * W_2d(rc * np.sqrt(q**2 - k**2 * mu**2)) * q

def pk_powerlaw(k, param): 
    return param[0] * (k/0.3)**param[1]

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

if __name__=="__main__": 
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
    
    for l in [0, 2]: 
        delp = delP_corr(
                l, 
                k, Pk, 
                extrap_params=[[3345.0, -1.6], [400.0, -4.], [260.0, -1.0]], 
                fs=1.0, rc=0.4
                )

        pickle.dump( delp, open('delP'+str(l)+'k_corr_different_integral.p', 'wb') )


