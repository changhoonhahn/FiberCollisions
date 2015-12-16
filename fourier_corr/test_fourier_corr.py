'''

Tests for the fourier_corr.py module 

'''
import time 
import fourier_corr

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
        sub.plot(q_arr, np.array(fllp_est), c=pretty_colors[0], lw=3, ls='--', label=est_label)

        maxmax = np.max([np.max([np.max(fllp_int), np.max(fllp_est)]), maxmax])
        minmin = np.min([np.min([np.min(fllp_int), np.min(fllp_est)]), minmin])

    sub.set_xscale('log')
    sub.set_xlabel(r"$\mathtt{q}$", fontsize=25)
    sub.set_ylabel(r"$\mathtt{f_{l, l'}(q r_c, k r_c)}$", fontsize=25)
    sub.set_ylim([1.1*minmin, 1.1*maxmax])

    sub.text(2.*10**-3, 1.01*maxmax, r"l = "+str(l)+", l' = "+str(lp), fontsize=20)
    sub.legend(loc='upper right')
    
    fig.savefig(''.join([
        'figure/', 
        'f_l_lp_estimate', 
        '.l', str(l), '.lp', str(lp), '.png'
        ]), bbox_inches='tight')
    plt.close()


if __name__=='__main__':
    for ell in [0]:
        for ellp in [6]:#, 8, 10]:  
            test_f_l_lp_est(ell, ellp, rc=0.43)
