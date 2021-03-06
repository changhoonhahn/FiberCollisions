'''

Fourier transform of 2PCF top hat Fiber collision Correction 

'''
import os
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

# --- Local --- 
from p_k_mu import delPcorr_pkmu_file
from p_k_mu import delPcorr_pkmu_save




def DelPuncorr(ell, fs=0.6, rc=0.43, k_arr=None, clobber=False): 
    ''' Uncorrelated Del P term of the tophat convolution correction. Returns list with 
    k values and Del P^uncorr
    
    Parameters
    ----------
    ell : int
        Specify the multipole. l = 0, 2 for monopole and quadrupole respectively. 

    fs : float
        Fraction of the survey that suffers from fiber collisions. 0 < fs < 1. 
        fs = 0.6 for Nseries and CMASS tiling. 

    rc : float
        Comoving physical scale of the fiber collision angular scale. 

    k_arr : ndarray (optional) 
        Array specifying the k values. The code is fast enough that this feature 
        should be utilized

    clobber : bool
        Boolean that specifies whether to remake the DelP^uncorr pickle file.  
    '''
    uncorr_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/',
        'uncorrdelP', str(ell), 'k', 
        '.fs', str(round(fs,2)), 
        '.rc', str(round(rc,2)), 
        '.p'
        ])

    if os.path.isfile(uncorr_file) and not clobber and k_arr is None: 
        k, delp = pickle.load(open(uncorr_file, 'rb'))  
    else: 
        if k_arr is None: 
            # if k values are not specified then take k values from arbitrary 
            # Nseries P(k) quadrupole 
            pk_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/', 
                'POWER_Q_CutskyN1.fidcosmo.dat.grid960.P020000.box3600'])
            k = np.loadtxt(true_pk_file, unpack = True, usecols =[0])
        else: 
            k = k_arr 
        delp = DelPuncorr_k(
                k,                  # k 
                ell,                # ell 
                fs=fs,              # fs
                rc=rc               # rc
                )

        if clobber: 
            pickle.dump([k, delp], open(uncorr_file, 'wb'))

    return [k, delp]

def DelPuncorr_k(k, l, fs=1.0, rc=0.4):
    ''' The uncorrelated term of the tophat convolved Del P as a function of k. The term 
    is computed as below: 

    Del P_uncorr(k)_l = - (2 pi f_s) (pi rc^2) (2 l + 1)/2 Leg_l(0) W_2d(k*rc)/k

    Parameters
    ----------
    k : float 
        k value to evaluate Del P^uncorr_l(k)
    l : int
        Monopole (l=0) or quadrupole (l=2) 
    fs : float
        The fraction of the survey that suffers from fiber collisions. fs = 0.6 for 
        Nseries and CMASS.
    rc : float
        The comoving physical scale of the fiber collision angule. At redshift z = 0.55
        rc = 0.43. (This is one of the key assumptions of the correction method). 
    '''
    leg_l = legendre(l)
    alpha = -1.0 * np.pi**2 * rc**2 * fs 
    delP = alpha * (2.*np.float(l)+1) * leg_l(0) * W_2d(k*rc)/k

    return delP

def DelPcorr_pkmu(ell, fs=0.6, rc=0.43, fold=10, rebin=20, clobber=False, Nthreads=1, 
        dqdmu=True, dmudq=False, ktrust=None):
    ''' Correlated Del P term of the top hat convolution correction calculated from the 
    Nseries Box P(k, mu).

    Parameters
    ----------
    ell : int
        Specify the multipole. l = 0, 2 for monopole and quadrupole respectively. 

    fs : float
        Fraction of the survey that suffers from fiber collisions. 0 < fs < 1. 
        fs = 0.6 for Nseries and CMASS tiling. 

    rc : float
        Comoving physical scale of the fiber collision angular scale. 

    fold : int 
        Number of folds in the Nseries box for P(k, mu) at extended k_max.

    rebin : int 
        Number of k bins for k > 0.1. The P(k, mu) is rebinned once the folded 
        P(k, mu) is included to reduce the noise. 

    clobber : bool
        Boolean that specifies whether to remake the DelP^uncorr pickle file.  

    Nthreads : int (optional)
        If clobber is True, Nthreads specifies the number of threads to used when 
        recomputing the DelP^corr from P(k, mu). Because it's a 2D integral, 
        multiple threads is recommended. The 2D integral is computed using  
        multiprocessing. 

    dqdmu : bool
        If True, the 2D integral is computed dq first then dmu later. Only dqdmu 
        or only dmudq can be true simultaneously. 

    dmudq : bool
        If True, the 2D integral is computed dmu first then dq later. Only dqdmu 
        or only dmudq can be true simultaneously. 

    ktrust : float (optional)
        
    '''
    corr_file = delPcorr_pkmu_file(ell, fs=fs, rc=rc, fold=fold, rebin=rebin, Ngrid=3600, 
            dqdmu=dqdmu, dmudq=dmudq, ktrust=ktrust) 

    if not os.path.isfile(corr_file) or clobber: 
        # compute Del P^corr from P(k, mu) 
        delPcorr_pkmu_save(ell, fs=fs, rc=rc, fold=fold, rebin=rebin, Ngrid=3600, 
                Nthreads=Nthreads)

    k, delp = pickle.load(open(corr_file, 'rb'))

    return [k, delp]

def DelPcorr_untrust_poly(k, ell=0, ktrust=0.3, order='all', fs=0.6, rc=0.43, fold=10, rebin=20, 
        max_order=10):
    ''' The untrusted portion of Del P^corr(k) calculated from polynomial expansion 
    '''
        
    fold_str = ''
    if fold: 
        fold_str = ''.join([
            '.', str(fold), 'fold', 
            '.', str(rebin), 'rebin'
            ])
    coeff_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'k^n_coeff.corrDelP_ktrust_qmax',
        '.ell', str(ell), 
        '.ktrust', str(round(ktrust,2)), 
        '.rc', str(round(rc,2)), 
        '.fs', str(round(fs,2)),
        fold_str, 
        '.order', str(max_order),
        '.dat'])

    n_k, n_coeff = np.loadtxt(coeff_file, unpack=True, skiprows=1, usecols=[0, 1])

    if order == 'all': 
        n_k_list = n_k 
        coeff_list = n_coeff

        poly = [coeff_list[ii] * k**n_k_list[ii] for ii in range(len(n_k_list))]

        return np.sum(np.vstack(np.array(poly)).T, axis=1)
    elif 'upto' in order: 
        n_order = np.where(n_k <= int(order.rsplit('upto')[-1])) 
        n_k_list = n_k[n_order]
        coeff_list = n_coeff[n_order]

        poly = [coeff_list[ii] * k**n_k_list[ii] for ii in range(len(n_k_list))]
        
        return np.sum(np.vstack(np.array(poly)).T, axis=1)
    else: 
        n_order = np.where(n_k == order)
        if len(n_order[0]) != 1: 
            raise ValueError
        return n_coeff[n_order] * k**order

def delP_corr_qrange(k, Pk, l, q_min=None, q_max=None, fs=1.0, rc=0.4, lp_comp=False): 
    '''
    Fiber collision correction term to Powerspectrum multipole l integrated
    over the q range [q_min, q_max]. The correlated term:
    
    delP_corr(k)_l = alpha * Sum_l' [int^qmax_qmin q P_l'(q) f_l,l'(q*rc,k*rc) dq]
    
    where alpha = - fs * rc^2/2 
    '''
    if q_min is None: 
        raise ValueError
    if q_max is None:
        raise ValueError
    if q_max > k[-1]: 
        raise ValueError
    if q_min < k[0]: 
        raise ValueError

    alpha = -0.5 * fs * rc**2

    if not isinstance(Pk, list): 
        raise ValueError(
                'Pk input has to be a list containing the multipoles. e.g. [P0k, P2k, P4k]')
    
    lvalues = [0, 2, 4, 6, 8, 10, 12, 14]

    sum_delP = np.zeros(len(k))
    
    if lp_comp: 
        lps, delPlp = [], [] 

    for i_lp, lp in enumerate(lvalues[:len(Pk)]):
        print 'delP^corr( l = ', l, ", l'= ", lp, ')'

        # Cubic Spline interpolated function of P(k)
        Pk_interp = interp1d(k, Pk[i_lp], kind='cubic')

        delP = np.zeros(len(k))

        # loop through k values and evaluate delta P for all ks 
        for i_k, k_i in enumerate(k): 
            dPq_integrand = lambda q_var: delPq_integrand(
                    q_var, 
                    k_i, 
                    l, 
                    lp, 
                    Pk_interp, 
                    None, 
                    rc=rc, 
                    k_min=k[0], 
                    k_max=k[-1], 
                    noextrap=True
                    )
            delP_int, delP_int_err = quad_int(dPq_integrand, q_min, q_max, epsrel=0.01)

            delP[i_k] = alpha * delP_int

        sum_delP += delP
         
        if lp_comp: 
            lps.append(lp)
            delPlp.append(delP)

    if lp_comp: 
        return [lps, delPlp]
    else:  
        return sum_delP 

def delP_corr_qmax(k, Pk, l, q_max=None, fs=1.0, rc=0.4, lp_comp=False): 
    '''
    Fiber collision correction term to Powerspectrum multipole l. The correlated term.
    
    delP_corr(k)_l = - (2 pi f_s) (2 pi^2 rc^2) Sum_l' [int^inf_0 q Pl'(q) f_ll'(q*rc,k*rc) dq]
    
    Leg_l(0) W_2d(k*rc)/k

    No extrapolation is assumed. q_max, the upper limit of the dq integral must be specified
    as keyword argument. 
    '''

    if q_max is None:
        raise ValueError
    if q_max > k[-1]: 
        raise ValueError

    alpha = -0.5 * fs * rc**2

    if not isinstance(Pk, list): 
        raise ValueError(
                'Pk input has to be a list containing the multipoles. e.g. [P0k, P2k, P4k]'
                )
    
    lvalues = [0, 2, 4, 6, 8, 10, 12, 14]

    sum_delP = np.zeros(len(k))
    
    if lp_comp: 
        lps, delPlp = [], [] 

    for i_lp, lp in enumerate(lvalues[:len(Pk)]):
        print 'delP^corr( l = ', l, ", l'= ", lp, ')'

        # Cubic Spline interpolated function of P(k)
        Pk_interp = interp1d(k, Pk[i_lp], kind='cubic')

        delP = np.zeros(len(k))

        # loop through k values and evaluate delta P for all ks 
        for i_k, k_i in enumerate(k): 
            dPq_integrand = lambda q_var: delPq_integrand(
                    q_var, 
                    k_i, 
                    l, 
                    lp, 
                    Pk_interp, 
                    None, 
                    rc=rc, 
                    k_min=k[0], 
                    k_max=k[-1], 
                    noextrap=True
                    )
            delP_int, delP_int_err = quad_int(dPq_integrand, 0., q_max, epsrel=0.01)

            delP[i_k] = alpha * delP_int

        sum_delP += delP
         
        if lp_comp: 
            lps.append(lp)
            delPlp.append(delP)

    if lp_comp: 
        return [lps, delPlp]
    else:  
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

    where theta = sqrt(x^2 - y^2*u^2). In this convention, x = q * rc and y = k * rc. 

    Below are estimates of the integral, without actually computing the integral 
    '''
    
    try: 
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
    except ValueError: 
        print 'Integral within integral'
        return f_l_lp(x, y, l, lp)

def fllp_poly(l1_in, l2_in, x):  
    '''
    Polygons for specific situations called by f_l_lp_est function
    '''
    # l2 = 0 
    if (l1_in == 2) and (l2_in == 0): 
        return x**2 - 1.
    elif (l1_in == 4) and (l2_in == 0): 
        return (7./4.)*x**4 - (5./2.)*x**2 + 3./4.
    elif (l1_in == 6) and (l2_in ==0): 
        return (33./8.)*x**6 - (63./8.)*x**4 + (35./8.)*x**2 - 5./8.
    elif (l1_in == 8) and (l2_in == 0): 
        # Poly(8,0,x) = 35/64 - (105 x^2)/16 + (693 x^4)/32 - (429 x^6)/16 + (715 x^8)/64
        return (715./64.)*x**8 - (429./16.)*x**6 + (693./32.)*x**4 - (105./16.)*x**2 + 35./64. 
    elif (l1_in == 10) and (l2_in == 0): 
        # Poly(10,0,x) =  -(63/128) + (1155 x^2)/128 - (3003 x^4)/64 + (6435 x^6)/64 - (12155 x^8)/128 + (4199 x^10)/128
        return (4199./128.)*x**10 - (12155./128.)*x**8 + (6435./64.)*x**6 - (3003./64.)*x**4 + (1155./128.)*x**2 - 63./128.
        # l2 = 2
    elif (l1_in == 4) and (l2_in == 2): 
        return x**4 - x**2
    elif (l1_in == 6) and (l2_in == 2): 
        return (11./4.)*x**6 - (9./2.)*x**4 + (7./4.)*x**2
    elif (l1_in == 8) and (l2_in == 2): 
        # Poly(8,2,x) = -((21 x^2)/8) + (99 x^4)/8 - (143 x^6)/8 + (65 x^8)/8
        return (65./8.)*x**8 - (143./8.)*x**6 + (99./8.)*x**4 - (21./8.)*x**2
    elif (l1_in == 10) and (l2_in == 2): 
        # Poly(10,2,x) =  (231 x^2)/64 - (429 x^4)/16 + (2145 x^6)/32 - (1105 x^8)/16 + (1615 x^10)/64
        return (1615./64.)*x**10 - (1105./16.)*x**8 + (2145./32.)*x**6 - (429./16.)*x**4 + (231./64.)*x**2

    elif (l1_in == 6) and (l2_in == 4): 
        return x**6 - x**4
    elif (l1_in == 8) and (l2_in == 4): 
        # Poly(8,4,x) =  (11 x^4)/4 - (13 x^6)/2 + (15 x^8)/4
        return (15./4.)*x**8 - (13./2.)*x**6 + (11./4.)*x**4
    elif (l1_in == 8) and (l2_in == 6): 
        # Poly(8,6,x) =  -x^6 + x^8
        return x**8 - x**6
    elif (l1_in == 10) and (l2_in == 4): 
        # Poly(10,4,x) =  -((143 x^4)/24) + (195 x^6)/8 - (255 x^8)/8 + (323 x^10)/24
        return (323./24.)*x**10 - (255./8.)*x**8 + (195./8.)*x**6 - (143./24.)*x**4
    elif (l1_in == 10) and (l2_in == 6): 
        # Poly(10,6,x) = (15 x^6)/4 - (17 x^8)/2 + (19 x^10)/4
        return (19./4.)*x**10 - (17./2)*x**8 + (15./4.)*x**6
    elif (l1_in == 10) and (l2_in == 8): 
        # Poly(10,8,x) = -x^8 + x^10
        return x**10 - x**8
    else: 
        raise ValueError

def f_l_lp(x, y, l, lp): #, error=False): 
    '''
    Actual integrals of the f_l,l' term in the FC correction equation 
    '''
    #start_time = time.time()
    upper_bound = np.min([1.0, x/y])
    lower_bound = np.max([-1.0, -1.0*x/y])
    
    fllp = lambda mu : f_l_lp_integrand(mu, x, y, l, lp)

    f_llp_integral = quad_int(fllp, lower_bound,  upper_bound, epsrel=0.01)
    f_llp = 0.5*(2.*np.float(l) + 1.) * f_llp_integral[0]
    
    #print "f_ll' calculation takes ", time.time() - start_time 
    return f_llp

def f_l_lp_integrand(mu, x, y, l, lp): 
    ''' Integrand of f_l_lp integration 
    '''
    Leg_l = legendre(l)     # Legendre polynomial of order l 
    Leg_lp = legendre(lp)   # Legendre polynomial of order l'
    
    if np.abs(mu) > np.abs(x/y): 
        raise ValueError

    theta = np.sqrt(x**2 - y**2 * mu**2)

    return Leg_l(mu) * Leg_lp(y * mu/x) * (W_2d(theta) + y**2 * (1. - mu**2) * W_secorder(theta))

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

def J2(x): 
    '''
    Spherical Bessel Function of the first kind of order 2 
    '''
    return jv(2, x)

"""
    def delP_corr(k, Pk, l, fs=1.0, rc=0.4, extrap_params=[[3345.0, -1.6], [400.0, -4.]], k_fixed=None, lp_comp=False, noextrap=False): 
        '''
        Fiber collision correction term to Powerspectrum multipole l. The correlated term.
        
        delP_corr(k)_l = - (2 pi f_s) (2 pi^2 rc^2) Sum_l' [int^inf_0 q Pl'(q) f_ll'(q*rc,k*rc) dq]
        
        Leg_l(0) W_2d(k*rc)/k

        extrap_params : [[alpha_0, n_0], [alpha_2, n_2], [alpha_4, n_4]]

        '''

        alpha = -0.5 * fs * rc**2       # this is to deal with the normalization 

        if not isinstance(Pk, list): 
            raise ValueError('Pk input has to be a list containing the multipoles. e.g. [P0k, P2k, P4k]')

        if len(Pk) != len(extrap_params): 
            raise ValueError('Extrapolation parameters need to be specified for each multipole') 
        
        lvalues = [0,2,4,6,8, 10, 12, 14]

        sum_delP = np.zeros(len(k))
        
        if lp_comp: 
            lps, delPlp = [], [] 
        if noextrap: 
            print 'no extrapolation'

        for i_lp, lp in enumerate(lvalues[:len(Pk)]):
            print 'delP^corr( l = ', l, ", l'= ", lp, ')'
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
                        k_max=k[-1], 
                        noextrap=noextrap
                        )
                if not noextrap: 
                    delP_int, delP_int_err = quad_int(dPq_integrand, 0., np.inf, epsrel=0.01)
                else: 
                    delP_int, delP_int_err = quad_int(dPq_integrand, 0., k[-1], epsrel=0.01)

                delP[i_k] = alpha * delP_int

            #print "delP(k) calculation takes ", time.time() - start_time 
            #print 'l=', l, " l'=", int(lp)
            #print delP 
            sum_delP += delP
             
            if lp_comp: 
                lps.append(lp)
                delPlp.append(delP)

        if lp_comp: 
            return [lps, delPlp]
        else:  
            return sum_delP 


        def W_1d(x,y): 
            '''
            W_1d(x,y) = int^2pi_0 W_2D(sqrt(x^2 + y^2 - 2xy cos(phi))) dphi/2pi
            '''
            integrand = lambda phi: W_2d( np.sqrt(x**2 + y**2 - 2*x*y*np.cos(phi)) )

            output = 1./(2.*np.pi) * quad_int(integrand, 0., 2.0 * np.pi, epsrel=0.05)[0]
            
            return output

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
"""
