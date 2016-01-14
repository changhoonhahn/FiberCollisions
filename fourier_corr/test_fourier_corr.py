'''

Tests for the fourier_corr.py module 

'''
import time 
import numpy as np

import mpfit
import pickle 
import pk_extrap
import fourier_corr

from scipy.special import j1 as J1
from scipy.interpolate import interp1d
from scipy.integrate import quad as quad_int 

# --- plotting ---
from defutility.plotting import prettyplot 
from defutility.plotting import prettycolors 


def test_f_l_lp_est(l, lp, rc=0.43):
    '''
    Test the polynomial estimates of the f_l_lp integrals. 

    Parameter
    ---------
    - l : 
    - lp : 

    Notes
    -----
    - Consistent for l, lp < 6
    - Possible break-down in the estimates for l or lp = 10? 
    '''
    q_arr = np.logspace(-3, 3, num=100)
    x_arr = q_arr * rc 
    
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(14,8))
    sub = fig.add_subplot(111)
    maxmax, minmin = 0., 0.
    for ik, k in enumerate([0.01, 0.05, 0.1, 0.5]):
        krc = k * rc

        fllp_est = []
        fllp_int = [] 
        for x in x_arr:
            #poly_time = time.time()
            if lp <= 10:
                fllp_est.append(fourier_corr.f_l_lp_est(x, krc, l, lp))
            #print 'polynomial estimate takes ', time.time() - poly_time
            #int_time = time.time()
            fllp_int.append(fourier_corr.f_l_lp(x, krc, l, lp))
            #print 'integration takes ', time.time() - int_time
        
        int_label = 'Integrated $ k = '+str(round(k, 2))+ '$'
        if ik == 0: 
            est_label = 'Polynomial Estimate'
        else: 
            est_label = None
        # plot the results 
        
        sub.plot(q_arr, np.array(fllp_int), c=pretty_colors[ik+1], lw=4, ls='-', label=int_label)
        if lp <= 10:
            sub.plot(q_arr, np.array(fllp_est), c=pretty_colors[0], lw=3, ls='--', label=est_label)
        else: 
            fllp_est = fllp_int

        maxmax = np.max([np.max([np.max(fllp_int), np.max(fllp_est)]), maxmax])
        minmin = np.min([np.min([np.min(fllp_int), np.min(fllp_est)]), minmin])

    sub.set_xscale('log')
    sub.set_xlabel(r"$\mathtt{q}$", fontsize=25)
    sub.set_ylabel(r"$\mathtt{f_{l, l'}(q r_c, k r_c)}$", fontsize=25)
    sub.set_ylim([-1.1, 1.1])

    sub.text(2.*10**-3, 1.01*maxmax, r"l = "+str(l)+", l' = "+str(lp), fontsize=20)
    sub.legend(loc='upper right')
    
    fig.savefig(''.join([
        'figure/', 
        'f_l_lp_estimate', 
        '.l', str(l), '.lp', str(lp), '.png'
        ]), bbox_inches='tight')
    plt.close()

def test_qPqfllp(l, lp, rc=0.43, noextrap=''):
    '''
    Test the polynomial estimates of the f_l_lp integrals. 

    Parameter
    ---------
    - l : 
    - lp : 

    Notes
    -----
    '''
    q_arr = np.logspace(-3, 3, num=100)

    n_mock = 7
    k_fit = 4.
    k_fixed = 4.34
    data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
    for i_mock in xrange(1, 8):
        true_pk_file = ''.join([data_dir, 'power3600z_BoxN', str(i_mock), '.dat'])
        if lp == 0: 
            l_index = -1
        else:
            l_index = 1 + int(lp/2)

        tr_k, tr_pk_i = np.loadtxt(
                    true_pk_file, 
                    unpack = True, 
                    usecols =[0,l_index] 
                    )
        if i_mock == 1: 
            tr_pk = tr_pk_i
        else: 
            tr_pk += tr_pk_i

    tr_pk /= 7.
    tr_specs =  (2.0*np.pi)**3 * tr_pk
    
    # interpolation function 
    Pk_interp = interp1d(tr_k, tr_pk, kind='cubic')
    
    # extrapolation parameter
    tr_extrap_par = pk_extrap.pk_powerlaw_bestfit(tr_k, tr_specs, k_fit=4., k_fixed=4.34)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(14,8))
    sub = fig.add_subplot(111)
    maxmax, minmin = 0., 0.
    for ik, k in enumerate([0.01, 0.05, 0.1, 0.5]):
        krc = k * rc

        qPqfllp = [] 
        for q in q_arr:
            if not noextrap: 
                Pq = pq(q, Pk_interp, tr_extrap_par, k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34)
            else:
                Pq = pq_noextrap(q, Pk_interp, k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34)

            fllp = fourier_corr.f_l_lp(q*rc, krc, l, lp)
            qPqfllp.append(
                    q * Pq * fllp
                    )
        
        int_label = '$ k = '+str(round(k, 2))+ '$'
        
        sub.plot(q_arr, np.array(qPqfllp), c=pretty_colors[ik+1], lw=4, ls='-', label=int_label)

        maxmax = np.max([np.max(qPqfllp), maxmax])
        minmin = np.min([np.min(qPqfllp), minmin])
    
    sub.set_xscale('log')
    sub.set_xlabel(r"$\mathtt{q}$", fontsize=25)
    sub.set_ylabel(r"$\mathtt{q P(q) f_{l, l'}(q r_c, k r_c)}$", fontsize=25)
    if l == 0: 
        sub.vlines(tr_k[-1], -20., 100., color='k', linestyles='--', linewidth=2)
        sub.set_ylim([-20, 20])
    elif l == 2: 
        sub.vlines(tr_k[-1], -25., 5., color='k', linestyles='--', linewidth=2)
        sub.set_ylim([-5, 5])

    sub.text(2.*10**-3, 1.01*maxmax, r"l = "+str(l)+", l' = "+str(lp), fontsize=20)
    sub.legend(loc='upper right')
    
    #plt.show()
    fig.savefig(''.join([
        'figure/', 
        'qPqfllp', noextrap, '.l', str(l), '.lp', str(lp), '.png'
        ]), bbox_inches='tight')
    plt.close()

def test_qPqfllp_k(k, l, rc=0.43, noextrap=''):
    '''
    Test the polynomial estimates of the f_l_lp integrals. 

    Parameter
    ---------
    - l : 
    - lp : 

    Notes
    -----
    '''
    q_arr = np.logspace(-3, 3, num=100)

    n_mock = 7
    krc = k * rc
    k_fit = 4.
    k_fixed = 4.34
    data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
    
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(1, figsize=(14,8))
    sub = fig.add_subplot(111)

    maxmax, minmin = 0., 0.

    for i_lp, ellp in enumerate(range(6)):
        lp = 2 * ellp

        for i_mock in xrange(1, 8):
            true_pk_file = ''.join([data_dir, 'power3600z_BoxN', str(i_mock), '.dat'])
            if lp == 0: 
                l_index = -1
            else:
                l_index = 1 + int(lp/2)

            tr_k, tr_pk_i = np.loadtxt(
                        true_pk_file, 
                        unpack = True, 
                        usecols =[0,l_index] 
                        )
            if i_mock == 1: 
                tr_pk = tr_pk_i
            else: 
                tr_pk += tr_pk_i

        tr_pk /= 7.
        tr_specs =  (2.0*np.pi)**3 * tr_pk
        
        # interpolation function 
        Pk_interp = interp1d(tr_k, tr_specs, kind='cubic')
        
        # extrapolation parameter
        tr_extrap_par = pk_extrap.pk_powerlaw_bestfit(tr_k, tr_specs, k_fit=4., k_fixed=4.34)

        qPqfllp = [] 
        Pqs = [] 
        for q in q_arr:
            if not noextrap: 
                Pq = pq(q, Pk_interp, tr_extrap_par, k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34)
            else:
                Pq = pq_noextrap(q, Pk_interp, k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34)

            fllp = fourier_corr.f_l_lp(q*rc, krc, l, lp)
            qPqfllp.append(
                    q * Pq * fllp
                    )
            Pqs.append(Pq)
        
        int_label = "$ l' = "+str(lp)+ "$"
        
        sub.plot(q_arr, np.array(qPqfllp), c=pretty_colors[i_lp+1], lw=4, ls='-', label=int_label)

        maxmax = np.max([np.max(qPqfllp), maxmax])
        minmin = np.min([np.min(qPqfllp), minmin])
    
    sub.set_xscale('log')
    sub.set_xlabel(r"$\mathtt{q}$", fontsize=25)
    sub.set_ylabel(r"$\mathtt{q P(q) f_{l, l'}(q r_c, k r_c)}$", fontsize=25)

    sub.vlines(tr_k[-1], minmin, maxmax, color='k', linestyles='--', linewidth=2)

    sub.text(2.*10**-3, 1.01*maxmax, r"k = "+str(round(k,2)), fontsize=20)
    sub.legend(loc='upper right')
    
    fig.savefig(''.join([
        'figure/', 
        'qPqfllp', noextrap, '.l', str(l), '.k', str(round(k,2)), '.png'
        ]), bbox_inches='tight')
    plt.close()

def test_fllp_k(k, l, rc=0.43):
    '''
    plot f_l_l' for l = 0, 2 and l' = 0, 2, 4, 6, 8, 10 in order to determine how it behaves

    Parameter
    ---------
    - l : ell value 

    Notes
    -----
    '''
    q_arr = np.logspace(-3, 3, num=100)

    n_mock = 7
    krc = k * rc
    k_fit = 4.
    k_fixed = 4.34
    data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
    
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(1, figsize=(14,8))
    sub = fig.add_subplot(111)

    maxmax, minmin = 0., 0.

    for i_lp, ellp in enumerate(range(20)):
        lp = 2 * ellp
        
        true_pk_file = ''.join([data_dir, 'power3600z_BoxN1.dat'])
        tr_k = np.loadtxt(true_pk_file, unpack=True, usecols =[0])

        fllp = [] 
        for q in q_arr:
            fllp.append(fourier_corr.f_l_lp(q*rc, krc, l, lp))
        
        int_label = "$ l' = "+str(lp)+ "$"
        
        sub.plot(q_arr, np.array(fllp), c=pretty_colors[i_lp+1], lw=4, ls='-', label=int_label)

        maxmax = np.max([np.max(fllp), maxmax])
        minmin = np.min([np.min(fllp), minmin])
    
    sub.set_xscale('log')
    sub.set_xlabel(r"$\mathtt{q}$", fontsize=25)
    sub.set_ylabel(r"$\mathtt{f_{l, l'}(q r_c, k r_c)}$", fontsize=25)

    sub.vlines(tr_k[-1], minmin, maxmax, color='k', linestyles='--', linewidth=2)

    sub.text(2.*10**-3, 1.01*maxmax, r"k = "+str(round(k,2)), fontsize=20)
    sub.legend(loc='upper right')
    
    fig.savefig(''.join([
        'figure/', 
        'fllp.l', str(l), '.k', str(round(k,2)), '.png'
        ]), bbox_inches='tight')
    plt.close()

def test_Pq_extrap():
    '''
    Test the extrapolation of P_l'(q)

    Parameter
    ---------

    Notes
    -----
    '''
    q_arr = np.logspace(-3, 3, num=100)

    n_mock = 7
    k_fit = 4.
    k_fixed = 4.34
    data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
    
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(2, figsize=(14,8))
    sub = fig.add_subplot(111)

    for i_lp, ellp in enumerate(range(6)):
        lp = 2 * ellp

        for i_mock in xrange(1, 8):
            true_pk_file = ''.join([data_dir, 'power3600z_BoxN', str(i_mock), '.dat'])
            if lp == 0: 
                l_index = -1
            else:
                l_index = 1 + int(lp/2)

            tr_k, tr_pk_i = np.loadtxt(
                        true_pk_file, 
                        unpack = True, 
                        usecols =[0,l_index] 
                        )
            if i_mock == 1: 
                tr_pk = tr_pk_i
            else: 
                tr_pk += tr_pk_i

        tr_pk /= 7.
        tr_specs =  (2.0*np.pi)**3 * tr_pk
        
        # interpolation function 
        Pk_interp = interp1d(tr_k, tr_specs, kind='cubic')
        
        # extrapolation parameter
        tr_extrap_par = pk_extrap.pk_powerlaw_bestfit(tr_k, tr_specs, k_fit=4., k_fixed=4.34)

        Pqs = np.array([
            pq(q, Pk_interp, tr_extrap_par, k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34) 
            for q in q_arr])
        int_label = "$ l' = "+str(lp)+ "$"
        sub.plot(q_arr, np.array(Pqs), c=pretty_colors[i_lp+1], lw=4, ls='-', label=int_label)
    
    sub.set_xscale('log')
    #sub.set_yscale('log')
    sub.set_xlabel(r"$\mathtt{q}$", fontsize=25)
    sub.set_ylabel(r"$\mathtt{P(q)}$", fontsize=25)
    sub.legend(loc='upper right')
    
    fig.savefig(''.join([
        'figure/', 
        'Pq_lp_noabs.png'
        ]), bbox_inches='tight')
    plt.close()

def pq(q, f_interp, extrap_param, k_min=0.0003, k_max=4.34, k_fixed=4.34): 
    if (q > k_min) and (q <= k_max): 
        return f_interp(q)
    elif (q <= k_min): 
        return 0.0
    elif (q > k_max): 
        return pk_extrap.pk_powerlaw(q, extrap_param, k_fixed=k_fixed)

def pq_noextrap(q, f_interp, k_min=0.0003, k_max=4.34, k_fixed=4.34): 
    if (q > k_min) and (q <= k_max): 
        return f_interp(q)
    elif (q <= k_min): 
        return 0.0
    elif (q > k_max): 
        return 0.0

def test_B_lp_integrand_nseriesbox(lp, rc=0.43, ktrust=0.5):
    '''
    *** THIS TEST IS INCORRECT; CODE IS LEFT AS REFERENCE ONLY ***
    
    Test the coefficient B_lp of the k^lp term of the polynomial 
    decomposition of delP_corr. Only coded for NseriesBox.

    B_lp integrand = P_lp(q)/q^lp J_1(q_lp * rc)

    Parameter
    ---------
    - lp : 

    Notes
    -----
    '''
    q_arr = np.logspace(-3, 3, num=100)

    n_mock = 7
    k_fit = 4.
    k_fixed = 4.34
    data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
    for i_mock in xrange(1, 8):
        true_pk_file = ''.join([data_dir, 'power3600z_BoxN', str(i_mock), '.dat'])
        if lp == 0: 
            l_index = -1
        else:
            l_index = 1 + int(lp/2)

        tr_k, tr_pk_i = np.loadtxt(
                    true_pk_file, 
                    unpack = True, 
                    usecols =[0,l_index] 
                    )
        if i_mock == 1: 
            tr_pk = tr_pk_i
        else: 
            tr_pk += tr_pk_i

    tr_pk /= 7.
    tr_specs =  (2.0*np.pi)**3 * tr_pk
    
    # interpolation function 
    Pk_interp = interp1d(tr_k, tr_pk, kind='cubic')
    # extrapolation parameter
    tr_extrap_par = pk_extrap.pk_powerlaw_bestfit(tr_k, tr_specs, k_fit=4., k_fixed=4.34)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(14,8))
    sub = fig.add_subplot(111)

    B_lp_integrand = [] 
    for q in q_arr:
        Pq = pq(q, Pk_interp, tr_extrap_par, k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34)

        B_lp_integrand.append(Pq / (q**lp) * J1(q * rc))
    
    ktrust_inf = np.where(q_arr > ktrust)
    sub.plot(
            q_arr[ktrust_inf], 
            np.array(B_lp_integrand)[ktrust_inf], 
            c=pretty_colors[lp+1], 
            lw=4, ls='-', label=r"$l' = "+str(lp)+"$")

    sub.set_xscale('log')
    sub.set_xlabel(r"$\mathtt{q}$", fontsize=25)
    sub.set_ylabel(r"$\mathtt{B_{l'}}$ integrand", fontsize=25)

    sub.vlines(ktrust, -5000, 5000, color='k', 
            linestyles='--', linewidth=2, label=r"$\mathtt{k_{trust}}$")
    sub.set_ylim([-10., 10.])
    sub.legend(loc='upper right')
    
    #plt.show()
    fig.savefig(''.join([
        'figure/', 
        'B_', str(lp), '_integrand_ktrust', str(round(ktrust,2)), '.png'
        ]), bbox_inches='tight')
    plt.close()

def test_k_n_coeff_integrand_nseriesbox(n_k, rc=0.43, ktrust=0.5):
    '''
    Test the coefficient of the k^n term of the polynomial 
    decomposition of delP_corr. Only coded for NseriesBox.

    Parameter
    ---------
    - n_k : 

    Notes
    -----
    '''
    q_arr = np.logspace(-3, 3, num=100)

    n_mock = 7
    k_fit = 4.
    k_fixed = 4.34
    data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'

    for i_mock in xrange(1, 8):
        true_pk_file = ''.join([data_dir, 'power3600z_BoxN', str(i_mock), '.dat'])
        tr_k, tr_pk0_i, tr_pk2_i, tr_pk4_i, tr_pk6_i, tr_pk8_i, tr_pk10_i = np.loadtxt(
                    true_pk_file, 
                    unpack = True, 
                    usecols =[0,-1,2,3,4,5,6] 
                    )
        if i_mock == 1: 
            tr_pk0 = tr_pk0_i
            tr_pk2 = tr_pk2_i
            tr_pk4 = tr_pk4_i
            tr_pk6 = tr_pk6_i
            tr_pk8 = tr_pk8_i
            tr_pk10= tr_pk10_i
        else: 
            tr_pk0 += tr_pk0_i
            tr_pk2 += tr_pk2_i
            tr_pk4 += tr_pk4_i
            tr_pk6 += tr_pk6_i
            tr_pk8 += tr_pk8_i
            tr_pk10+= tr_pk10_i

    tr_pk0 *= (2.0*np.pi)**3/7.
    tr_pk2 *= (2.0*np.pi)**3/7.
    tr_pk4 *= (2.0*np.pi)**3/7.
    tr_pk6 *= (2.0*np.pi)**3/7.
    tr_pk8 *= (2.0*np.pi)**3/7.
    tr_pk10*= (2.0*np.pi)**3/7.
    tr_specs = [tr_pk0 ,tr_pk2 ,tr_pk4 ,tr_pk6 ,tr_pk8 ,tr_pk10]
    
    Pk_interps = [] 
    tr_extrap_pars = [] 
    for i_spec, tr_spec in enumerate(tr_specs):
        # interpolation function 
        Pk_interps.append(interp1d(tr_k, tr_spec, kind='cubic'))
        tr_extrap_pars.append(pk_extrap.pk_powerlaw_bestfit(tr_k, tr_spec, k_fit=4., k_fixed=4.34))

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(14,8))
    sub = fig.add_subplot(111)
    
    alpha = -1.0 * rc**2

    if n_k == 2: 
        Clp = [0., 1., -2.5, 2.5*7./4., 2.5*(-21./8.), 2.5*(231./64.)]
    elif n_k == 4: 
        Clp = [0., 0., 2.5, 2.5*(-9./2.), 2.5*(99./8.), 2.5*(-429./16.)]
    elif n_k == 6: 
        Clp = [0., 0., 0., 2.5*(11./4.), 2.5*(-143./8.), 2.5*(2145./32.)]
    elif n_k == 8: 
        Clp = [0., 0., 0., 0., 2.5*(65./8.), 2.5*(-1105./16.)]
    print Clp

    ktrust_inf = np.where(q_arr > ktrust)
    q_krange = q_arr[ktrust_inf]
    
    # plot the integrands of the coefficients to get an idea of 
    # what they look like 
    for i_spec, tr_spec in enumerate(tr_specs):
        coeff = np.zeros(len(q_krange))
        for i_q, q in enumerate(q_krange):
            Pq = pq(q, Pk_interps[i_spec], tr_extrap_pars[i_spec], 
                    k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34)
            coeff[i_q] = (Clp[i_spec] * Pq) * (J1(q * rc) / (q**n_k))

        sub.plot(
                q_krange, 
                alpha * coeff, 
                c=pretty_colors[i_spec+1], 
                lw=2, ls='--', label=str(2*i_spec))
        
        if i_spec == 0: 
            integrand = coeff
        else: 
            integrand += coeff
    
    sub.plot(
            q_arr[ktrust_inf], 
            alpha * integrand, 
            c=pretty_colors[n_k+1], 
            lw=4, ls='-', label=r"$l = 2$")
    
    # actually integrate to calculate the coefficients 
    coeff = 0.
    for i_spec, tr_spec in enumerate(tr_specs):
        Pq = pq(q, Pk_interps[i_spec], tr_extrap_pars[i_spec], 
                k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34)
    
        if Clp[i_spec] != 0.: 
            coeff_integ = lambda qq: \
                    Clp[i_spec] * \
                    pq(qq, Pk_interps[i_spec], tr_extrap_pars[i_spec], k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34) / \
                    (qq**n_k)

            coeff_int, coeff_int_err = quad_int(coeff_integ, ktrust, np.inf, epsrel=0.01)
        else: 
            coeff_int = 0.

        coeff += coeff_int

    coeff *= alpha

    sub.text(10., 10., str(round(coeff, 2))+r"$\mathtt{k^{"+str(n_k)+"}}$", fontsize=25)
    sub.set_xscale('log')
    sub.set_xlabel(r"$\mathtt{q}$", fontsize=25)
    sub.set_ylabel(
            r"$\mathtt{\Delta P^{corr}}$, $\mathtt{k^{"+str(n_k)+"}}$ Coefficient Integrand", 
            fontsize=25)

    #sub.vlines(ktrust, 1.1*np.min(integrand), 1.1*np.max(integrand), color='k', 
    #        linestyles='--', linewidth=2, label=r"$\mathtt{k_{trust}}$")
    #sub.set_ylim([-10., 10.])
    sub.legend(loc='upper right')
    
    #plt.show()
    fig.savefig(''.join([
        'figure/', 
        'k^', str(n_k), '_coeff_integrand.ktrust', str(round(ktrust,2)), '.png'
        ]), bbox_inches='tight')
    plt.close()

def test_k_n_coeff(ktrust, rc=0.43, ell=2):
    '''
    Calculate the coefficient of the k^n term of the polynomial 
    decomposition of delP_corr. Only coded for NseriesBox.

    Parameter
    ---------

    Notes
    -----
    '''
    n_mock = 7
    k_fit = 4.
    k_fixed = 4.34
    data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
    for i_mock in xrange(1, 8):
        true_pk_file = ''.join([data_dir, 'power3600z_BoxN', str(i_mock), '.dat'])
        tr_k, tr_pk0_i, tr_pk2_i, tr_pk4_i, tr_pk6_i, tr_pk8_i, tr_pk10_i = np.loadtxt(
                    true_pk_file, 
                    unpack = True, 
                    usecols =[0,-1,2,3,4,5,6] 
                    )
        if i_mock == 1: 
            tr_pk0 = tr_pk0_i
            tr_pk2 = tr_pk2_i
            tr_pk4 = tr_pk4_i
            tr_pk6 = tr_pk6_i
            tr_pk8 = tr_pk8_i
            tr_pk10= tr_pk10_i
        else: 
            tr_pk0 += tr_pk0_i
            tr_pk2 += tr_pk2_i
            tr_pk4 += tr_pk4_i
            tr_pk6 += tr_pk6_i
            tr_pk8 += tr_pk8_i
            tr_pk10+= tr_pk10_i

    tr_pk0 *= (2.0*np.pi)**3/7.
    tr_pk2 *= (2.0*np.pi)**3/7.
    tr_pk4 *= (2.0*np.pi)**3/7.
    tr_pk6 *= (2.0*np.pi)**3/7.
    tr_pk8 *= (2.0*np.pi)**3/7.
    tr_pk10*= (2.0*np.pi)**3/7.
    tr_specs = [tr_pk0 ,tr_pk2 ,tr_pk4 ,tr_pk6 ,tr_pk8 ,tr_pk10]
    
    Pk_interps = [] 
    tr_extrap_pars = [] 
    for i_spec, tr_spec in enumerate(tr_specs):
        # interpolation function 
        Pk_interps.append(interp1d(tr_k, tr_spec, kind='cubic'))
        tr_extrap_pars.append(pk_extrap.pk_powerlaw_bestfit(tr_k, tr_spec, k_fit=4., k_fixed=4.34))

    alpha = -1.0 * rc**2
    
    n_ks = [] 
    n_k_coeff = [] 
    for n_k in [2, 4, 6, 8, 10]: 
        if ell == 2: 
            if n_k == 2: 
                Clp = [0., 1., -2.5, 2.5*7./4., 2.5*(-21./8.), 2.5*(231./64.)]
            elif n_k == 4: 
                Clp = [0., 0., 2.5, 2.5*(-9./2.), 2.5*(99./8.), 2.5*(-429./16.)]
            elif n_k == 6: 
                Clp = [0., 0., 0., 2.5*(11./4.), 2.5*(-143./8.), 2.5*(2145./32.)]
            elif n_k == 8: 
                Clp = [0., 0., 0., 0., 2.5*(65./8.), 2.5*(-1105./16.)]
            elif n_k == 10: 
                Clp = [0., 0., 0., 0., 0., 2.5*(1615./64.)]
            else: 
                raise NotImplementedError
        else:
            raise NotImplementedError('Only quadrupole correction implemented') 

        # actually integrate to calculate the coefficients 
        coeff = 0.
        for i_spec, tr_spec in enumerate(tr_specs):

            if Clp[i_spec] != 0.: 
                coeff_integ = lambda qq: \
                        Clp[i_spec] * \
                        pq(qq, Pk_interps[i_spec], tr_extrap_pars[i_spec], k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34) / \
                        (qq**n_k)

                coeff_int, coeff_int_err = quad_int(coeff_integ, ktrust, np.inf, epsrel=0.01)
            else: 
                coeff_int = 0.

            coeff += coeff_int

        coeff *= alpha

        n_ks.append(n_k)
        n_k_coeff.append(coeff)
    
    data_list = [np.array(n_ks), np.array(n_k_coeff)]

    coeff_file = ''.join([
        data_dir, 
        'k^n_coeff',
        '.ell', str(ell), 
        '.ktrust', str(round(ktrust,2)), 
        '.dat'])
    np.savetxt(coeff_file, 
            (np.vstack(np.array(data_list))).T, 
            fmt=['%10.5f', '%10.5f'], 
            delimiter='\t', 
            header='Columns : n_k, coeff'
            )

    return None

def test_k_n_coeff_noextrap(ktrust, rc=0.43, ell=2):
    '''
    Calculate the coefficient of the k^n term of the polynomial 
    decomposition of delP_corr. Only coded for NseriesBox.

    Parameter
    ---------

    Notes
    -----
    '''
    n_mock = 7
    k_fit = 4.
    k_fixed = 4.34
    data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
    for i_mock in xrange(1, 8):
        true_pk_file = ''.join([data_dir, 'power3600z_BoxN', str(i_mock), '.dat'])
        tr_k, tr_pk0_i, tr_pk2_i, tr_pk4_i, tr_pk6_i, tr_pk8_i, tr_pk10_i = np.loadtxt(
                    true_pk_file, 
                    unpack = True, 
                    usecols =[0,-1,2,3,4,5,6] 
                    )
        if i_mock == 1: 
            tr_pk0 = tr_pk0_i
            tr_pk2 = tr_pk2_i
            tr_pk4 = tr_pk4_i
            tr_pk6 = tr_pk6_i
            tr_pk8 = tr_pk8_i
            tr_pk10= tr_pk10_i
        else: 
            tr_pk0 += tr_pk0_i
            tr_pk2 += tr_pk2_i
            tr_pk4 += tr_pk4_i
            tr_pk6 += tr_pk6_i
            tr_pk8 += tr_pk8_i
            tr_pk10+= tr_pk10_i

    tr_pk0 *= (2.0*np.pi)**3/7.
    tr_pk2 *= (2.0*np.pi)**3/7.
    tr_pk4 *= (2.0*np.pi)**3/7.
    tr_pk6 *= (2.0*np.pi)**3/7.
    tr_pk8 *= (2.0*np.pi)**3/7.
    tr_pk10*= (2.0*np.pi)**3/7.
    tr_specs = [tr_pk0 ,tr_pk2 ,tr_pk4 ,tr_pk6 ,tr_pk8 ,tr_pk10]
    
    Pk_interps = [] 
    for i_spec, tr_spec in enumerate(tr_specs):
        # interpolation function 
        Pk_interps.append(interp1d(tr_k, tr_spec, kind='cubic'))

    alpha = -1.0 * rc**2
    
    n_ks = [] 
    n_k_coeff = [] 
    for n_k in [2, 4, 6, 8, 10]: 
        if ell == 2: 
            if n_k == 2: 
                Clp = [0., 1., -2.5, 2.5*7./4., 2.5*(-21./8.), 2.5*(231./64.)]
            elif n_k == 4: 
                Clp = [0., 0., 2.5, 2.5*(-9./2.), 2.5*(99./8.), 2.5*(-429./16.)]
            elif n_k == 6: 
                Clp = [0., 0., 0., 2.5*(11./4.), 2.5*(-143./8.), 2.5*(2145./32.)]
            elif n_k == 8: 
                Clp = [0., 0., 0., 0., 2.5*(65./8.), 2.5*(-1105./16.)]
            elif n_k == 10: 
                Clp = [0., 0., 0., 0., 0., 2.5*(1615./64.)]
            else: 
                raise NotImplementedError
        else:
            raise NotImplementedError('Only quadrupole correction implemented') 

        # actually integrate to calculate the coefficients 
        coeff = 0.
        for i_spec, tr_spec in enumerate(tr_specs):

            if Clp[i_spec] != 0.: 
                coeff_integ = lambda qq: \
                        Clp[i_spec] * \
                        pq_noextrap(qq, Pk_interps[i_spec], k_min=tr_k[0], k_max=tr_k[-1], k_fixed=4.34) / \
                        (qq**n_k)

                coeff_int, coeff_int_err = quad_int(coeff_integ, ktrust, np.inf, epsrel=0.01)
            else: 
                coeff_int = 0.

            coeff += coeff_int

        coeff *= alpha

        n_ks.append(n_k)
        n_k_coeff.append(coeff)
    
    data_list = [np.array(n_ks), np.array(n_k_coeff)]

    coeff_file = ''.join([
        data_dir, 
        'k^n_coeff',
        '.ell', str(ell), 
        '.ktrust', str(round(ktrust,2)), 
        '.noextrap', 
        '.dat'])
    np.savetxt(coeff_file, 
            (np.vstack(np.array(data_list))).T, 
            fmt=['%10.5f', '%10.5f'], 
            delimiter='\t', 
            header='Columns : n_k, coeff'
            )

    return None

def delPcorr_ktrust_inf_bestfit(ell, k_trust=0.5, Ngrid=960):
    '''
    Investigate the contribution of delPcorr integrated over the q range 
    k_trust to infinity. Fit this contribution to a polynomial of the form
    
    C_2 * k^2 + C_4 * k^4 + C_6 * k^6 + ...

    delPcorr is being divided into two parts:
        - 0 to k_trust, which can be reliably calculated 
        - k_trust to infinity, which cannot be reliably calculated because 
        we do not trust the models or the extrapolations beyound this point.
        However they it may be possible to have a polynomial fit. 
    '''
    # total DelP calculated from Nseries mocks (NOT THE BOX MOCKS)
    DelP_tot_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'DelP_', str(ell), 'k_tot.nseries.p'])
    k_upw, DelP_tot = pickle.load(open(DelP_tot_file, 'rb'))
    DelP_tot_k = interp1d(k_upw, DelP_tot)

    # subtract out DelP_uncorr 
    DelP_uncorr_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'uncorrdelP', str(ell), 'k_AVG_power3600z_BoxN.dat.p'
        ])
    k_uncorr, DelP_uncorr = pickle.load(open(DelP_uncorr_file, 'rb'))
    
    krange = np.where(
            (k_uncorr > k_upw[0]) & (k_uncorr < k_upw[-1])
            )
    k_val = k_uncorr[krange]

    DelP_corr = DelP_tot_k(k_val) - DelP_uncorr[krange]

    lps = [0, 2, 4, 6, 8, 10]
    for i_lp, lp in enumerate(lps):
        # delP^corr integrated over q = 0 - k_trust
        pickle_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/Box/', 
            'corrdelP', str(ell), 'k_lp', str(lp), 
            '_qmax', 
            str(round(k_trust,2)),
            '_AVG_power3600z_BoxN.dat.p'
            ])
        k_lp_0_ktrust, corrdelP_lp_0_ktrust = pickle.load(open(pickle_file, 'rb'))
        
        if i_lp == 0: 
            DelP_corr_0_ktrust = corrdelP_lp_0_ktrust
        else: 
            DelP_corr_0_ktrust += corrdelP_lp_0_ktrust
    
    DelP_corr_ktrust_inf = DelP_corr - DelP_corr_0_ktrust[krange]

    fa = {'x': k_val, 'y': DelP_corr_ktrust_inf}

    coeff_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/',
        'k^n_coeff',
        '.ell', str(ell), 
        '.ktrust', str(round(k_trust,2)), 
        '.noextrap.dat'])
    n_ks, n_ks_coeff = np.loadtxt(coeff_file, skiprows=1, unpack=True, usecols=[0, 1])
    param_guess = n_ks_coeff
    bestfit = mpfit.mpfit(delp_poly_mpfit, param_guess, functkw=fa, quiet=True)

    return [n_ks, bestfit.params]

def delp_poly(k, coeff): 
    ''' 
    Generic polynomial (starts with k^2 instead of k^0 because there 
    shouldn't be a k^0 term). 
    '''
    n_coeff = len(coeff)
    poly = [coeff[ii] * k**(2 * (ii+1)) for ii in range(n_coeff)]
    
    return np.sum(np.vstack(np.array(poly)).T, axis=1)

def delp_poly_mpfit(coeff, fjac=None, x=None, y=None): 
    '''
    Mpfit wrapper
    '''
    model = delp_poly(x, coeff)
    status = 0 
    return ([status, y-model])




if __name__=='__main__':
    print delPcorr_ktrust_inf_bestfit(2, k_trust=0.5, Ngrid=960)

    #for ktrust in [0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 2.0, 3.0, 4.0]:
    #    test_k_n_coeff_noextrap(ktrust)
    #test_k_n_coeff_integrand_nseriesbox(h, rc=0.43, ktrust=0.5)
    #test_Pq_extrap()
    #for k_i in [0.3, 0.7]: 
    #    #test_qPqfllp_k(k_i, 0, rc=0.43, noextrap='_noextrap')
    #    #test_qPqfllp_k(k_i, 2, rc=0.43, noextrap='_noextrap')
    #    test_fllp_k(k_i, 0, rc=0.43)
    #    test_fllp_k(k_i, 2, rc=0.43)
    #    #test_qPqfllp_k(k_i, 0, rc=0.43)
    #    #test_qPqfllp_k(k_i, 2, rc=0.43)
    '''
    for ell in [0,2]:
        for ellp in range(6):  
            test_qPqfllp(ell, 2*ellp, rc=0.43, noextrap='_noextrap')
    #        test_f_l_lp_est(ell, 2*ellp, rc=0.43)
    '''
