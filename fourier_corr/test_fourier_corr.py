'''

Tests for the fourier_corr.py module 

'''
import time 
import pk_extrap
import fourier_corr

from scipy.interpolate import interp1d

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

    for i_lp, ellp in enumerate(range(6)):
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
        sub.plot(q_arr, np.abs(np.array(Pqs)), c=pretty_colors[i_lp+1], lw=4, ls='-', label=int_label)
    
    sub.set_xscale('log')
    sub.set_yscale('log')
    sub.set_xlabel(r"$\mathtt{q}$", fontsize=25)
    sub.set_ylabel(r"$\mathtt{P(q)}$", fontsize=25)
    sub.legend(loc='upper right')
    
    fig.savefig(''.join([
        'figure/', 
        'Pq_lp.png'
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

if __name__=='__main__':
    for k_i in [0.1, 0.3, 0.5, 0.7]: 
        #test_qPqfllp_k(k_i, 0, rc=0.43, noextrap='_noextrap')
        #test_qPqfllp_k(k_i, 2, rc=0.43, noextrap='_noextrap')
        test_fllp_k(k_i, 0, rc=0.43)
        test_fllp_k(k_i, 2, rc=0.43)
        #test_qPqfllp_k(k_i, 0, rc=0.43)
        #test_qPqfllp_k(k_i, 2, rc=0.43)

    '''
    for ell in [0,2]:
        for ellp in range(6):  
            test_qPqfllp(ell, 2*ellp, rc=0.43, noextrap='_noextrap')
    #        test_f_l_lp_est(ell, 2*ellp, rc=0.43)
    '''
