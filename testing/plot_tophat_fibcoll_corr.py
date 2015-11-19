'''

Plots for xi tohat fiber collision correction methods 

'''
import time
import pickle
import numpy as np 
import scipy as sp

from scipy.interpolate import interp1d
import tophat_fibcoll_corr_test as tophat
from defutility.plotting import prettyplot 
from defutility.plotting import prettycolors 

def plot_delP(l, fs=1.0, rc=0.4, n_mocks=84):
    '''
    Comparison of delP^corr, delP^uncorr, and P^upw_avg - P^true_avg

    '''
    if l == 0: 
        l_cols = [0, 1]
    elif l == 2: 
        l_cols = [0, 2]
    elif l == 4: 
        l_cols = [0, 4]
    else: 
        raise NotImplementedError()
    
    # Calculate average P(k) for upweight and true cases
    for i_mock in xrange(1, n_mocks+1): 

        # true P_l(k)
        k, Pk_i = np.loadtxt(
                ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/', 
                    'POWER_Q_CutskyN'+str(i_mock)+'.fidcosmo.dat.grid360.P020000.box3600'
                    ]),
                unpack = True, 
                usecols = l_cols 
                ) 
        
        # Upweighted P_l(k)
        k, Pk_upw_i = np.loadtxt(
                ''.join(['/mount/riachuelo1/hahn/power/Nseries/',
                    'POWER_Q_CutskyN'+str(i_mock)+'.fidcosmo.fibcoll.dat.grid360.P020000.box3600'
                    ]),
                unpack = True, 
                usecols = l_cols
                ) 

        if i_mock == 1: 
            Pk_sum = Pk_i
            Pk_upw_sum = Pk_upw_i
        else: 
            Pk_sum += Pk_i
            Pk_upw_sum += Pk_upw_i

    Pk = Pk_sum/np.float(n_mocks)
    Pk_upw = Pk_upw_sum/np.float(n_mocks)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    corrdelP = pickle.load(open('delP'+str(l)+'k_corr_estimated.p', 'rb'))
    #corrdelP_diff = pickle.load(open('delP'+str(l)+'k_corr_different_integral.p', 'rb'))
    #print corrdelP[np.where(corrdelP != 0.0)]
    #print corrdelP_diff[np.where(corrdelP_diff != 0.0)]

    corrdelP_lp = pickle.load(open('delP'+str(l)+'k_corr_estimated_lpcomponent.p', 'rb'))[0]
    corrdelP_comp = pickle.load(open('delP'+str(l)+'k_corr_estimated_lpcomponent.p', 'rb'))[1]
    corrdelP_noextrap = pickle.load(open('delP'+str(l)+'k_corr_estimated_noextrap.p', 'rb'))
    corrdelP_noextrap_lp = pickle.load(open('delP'+str(l)+'k_corr_estimated_lpcomponent_noextrap.p', 'rb'))[0]
    corrdelP_noextrap_comp = pickle.load(open('delP'+str(l)+'k_corr_estimated_lpcomponent_noextrap.p', 'rb'))[1]
    #corrdelP_noextrap_firstorder = pickle.load(open('delP'+str(l)+'k_corr_noextrap_firstorder.p', 'rb'))
    uncorrdelP = tophat.delP_uncorr(k, l, fs=fs, rc=rc)
    
    sub.plot(
            k, 
            uncorrdelP, 
            c=pretty_colors[3], 
            lw=4, 
            label="Uncorrelated"
            )

    sub.plot(
            k, 
            corrdelP, 
            c = pretty_colors[1], 
            lw = 4, 
            label = "Correlated"
            )
    
    sub.scatter(
            k, 
            corrdelP_noextrap, 
            c = pretty_colors[2], 
            label = "Correlated (No extrap.)"
            )

    for i_lp in xrange(len(corrdelP_comp)): 
    
        sub.plot(
                k, corrdelP_comp[i_lp],
                c = pretty_colors[7+i_lp],
                lw = 2, 
                ls = '--',
                label = "Correlated: "+r"$l' ="+str(corrdelP_lp[i_lp])+"$"
                )
    
    for i_lp in xrange(len(corrdelP_noextrap_comp)): 
    
        sub.scatter(
                k, corrdelP_noextrap_comp[i_lp],
                c = pretty_colors[7+i_lp],
                label = "No extrapolation: "+r"$l' ="+str(corrdelP_noextrap_lp[i_lp])+"$"
                )
    #nonzero = np.where(corrdelP_noextrap_firstorder != 0.0)
    #sub.scatter(
    #        k[nonzero], 
    #        corrdelP_noextrap_firstorder[nonzero]/(2.0 * np.pi)**0., 
    #        c = pretty_colors[1], 
    #        label = "Correlated (No extrap. First order)"
    #        )

    #sub.plot(
    #        k, 
    #        (uncorrdelP + corrdelP)/(2.0 * np.pi)**0., 
    #        c = pretty_colors[5], 
    #        lw = 4, 
    #        ls = '-.',
    #        label = "Uncorrelated + Correlated Combined "
    #        )
    sub.plot(
            k, 
            Pk_upw - Pk, 
            c = 'k', 
            lw = 4,
            ls = '--', 
            label = r"$\mathtt{P^{upw}(k) - P^{true}(k)}$"
            )

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    
    #if l == 0: 
    #    sub.set_ylim([-40000.0, 5000.0])
    #elif l == 2:
    #    sub.set_ylim([-100.0, 2000.0])
    sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(l)+"}(k)}$", fontsize=30)
    
    if l == 0: 
        sub.legend(loc='lower right', scatterpoints = 1)
    elif l == 2: 
        sub.legend(loc='upper right', scatterpoints = 1)
    elif l == 4: 
        sub.legend(loc='lower right', scatterpoints = 1)
    fig.savefig('qaplot_delP_'+str(l)+'.png', bbox_inches="tight")
    plt.show()
    plt.close()

def plot_delP_extrapolation_test(l, type='normal', fs=1.0, rc=0.4, n_mocks=84):
    '''
    Comparison of delP^corr, delP^uncorr, and P^upw_avg - P^true_avg

    '''
    if l == 0: 
        l_cols = [0, 1]
    elif l == 2: 
        l_cols = [0, 2]
    elif l == 4: 
        l_cols = [0, 4]
    else: 
        raise NotImplementedError()
    
    # Calculate average P(k) for upweight and true cases
    for i_mock in xrange(1, n_mocks+1): 

        # true P_l(k)
        k, Pk_i = np.loadtxt(
                ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/', 
                    'POWER_Q_CutskyN'+str(i_mock)+'.fidcosmo.dat.grid360.P020000.box3600'
                    ]),
                unpack = True, 
                usecols = l_cols 
                ) 
        
        # Upweighted P_l(k)
        k, Pk_upw_i = np.loadtxt(
                ''.join(['/mount/riachuelo1/hahn/power/Nseries/',
                    'POWER_Q_CutskyN'+str(i_mock)+'.fidcosmo.fibcoll.dat.grid360.P020000.box3600'
                    ]),
                unpack = True, 
                usecols = l_cols
                ) 

        if i_mock == 1: 
            Pk_sum = Pk_i
            Pk_upw_sum = Pk_upw_i
        else: 
            Pk_sum += Pk_i
            Pk_upw_sum += Pk_upw_i

    Pk = Pk_sum/np.float(n_mocks)
    Pk_upw = Pk_upw_sum/np.float(n_mocks)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    corrdelP = pickle.load(open('delP'+str(l)+'k_corr_estimated.p', 'rb'))
    corrdelP_noextrap = pickle.load(open('delP'+str(l)+'k_corr_estimated_noextrap.p', 'rb'))
    corrdelP_lp = pickle.load(open('delP'+str(l)+'k_corr_estimated_lpcomponent.p', 'rb'))[0]
    corrdelP_comp = pickle.load(open('delP'+str(l)+'k_corr_estimated_lpcomponent.p', 'rb'))[1]
    corrdelP_noextrap_lp = pickle.load(open('delP'+str(l)+'k_corr_estimated_lpcomponent_noextrap.p', 'rb'))[0]
    corrdelP_noextrap_comp = pickle.load(open('delP'+str(l)+'k_corr_estimated_lpcomponent_noextrap.p', 'rb'))[1]
    #corrdelP_noextrap_firstorder = pickle.load(open('delP'+str(l)+'k_corr_noextrap_firstorder.p', 'rb'))
    uncorrdelP = tophat.delP_uncorr(k, l, fs=fs, rc=rc)
    
    if type == 'normal': 
        sub.plot(
                k, 
                corrdelP, 
                c = pretty_colors[1], 
                lw = 4, 
                label = "Correlated"
                )
        
        sub.scatter(
                k, 
                corrdelP_noextrap, 
                c = pretty_colors[2], 
                label = "Correlated (No extrap.)"
                )

        for i_lp in xrange(len(corrdelP_comp)): 
        
            sub.plot(
                    k, corrdelP_comp[i_lp],
                    c = pretty_colors[7+i_lp],
                    lw = 2, 
                    ls = '--',
                    label = "Correlated: "+r"$l' ="+str(corrdelP_lp[i_lp])+"$"
                    )
        
        for i_lp in xrange(len(corrdelP_noextrap_comp)): 
        
            sub.scatter(
                    k, corrdelP_noextrap_comp[i_lp],
                    c = pretty_colors[7+i_lp],
                    label = "No extrapolation: "+r"$l' ="+str(corrdelP_noextrap_lp[i_lp])+"$"
                    )
    elif type == 'difference': 
        sub.plot(
                k, 
                corrdelP - corrdelP_noextrap, 
                c = pretty_colors[1], 
                lw = 4, 
                )

        for i_lp in xrange(len(corrdelP_comp)): 
        
            sub.plot(
                    k, 
                    corrdelP_comp[i_lp] - corrdelP_noextrap_comp[i_lp],
                    c = pretty_colors[7+i_lp],
                    lw = 2, 
                    ls = '--',
                    label = r"$l' ="+str(corrdelP_lp[i_lp])+"$"
                    )

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    
    #if l == 0: 
    #    sub.set_ylim([-40000.0, 5000.0])
    #elif l == 2:
    #    sub.set_ylim([-100.0, 2000.0])
    if type == 'normal':
        sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(l)+"}^{corr}(k)}$", fontsize=30)
    elif type == 'difference': 
        sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(l)+"}^{corr}(k) - \Delta P_{"+str(l)+"}^{corr}(k)_{no\;extrap}}$", fontsize=30)
    
    if l == 0: 
        sub.legend(loc='lower left', scatterpoints = 1)
    elif l == 2: 
        sub.legend(loc='upper left', scatterpoints = 1)
    elif l == 4: 
        sub.legend(loc='lower left', scatterpoints = 1)
     
    if type == 'normal':
        fig.savefig('qaplot_delPcorr_'+str(l)+'_extrapolationtest.png', bbox_inches="tight")
    if type == 'difference':
        fig.savefig('qaplot_delPcorr_'+str(l)+'_extrapolationtest_diff.png', bbox_inches="tight")
    plt.show()
    plt.close()

def plot_corrected_Pk(l, fs=1.0, rc=0.4):
    '''
    '''

    if l == 0: 
        l_cols = [0, 1]
    else: 
        l_cols = [0, 2]
    
    # true P_l(k)
    k, Pk = np.loadtxt(
            'POWER_Q_CutskyN1.fidcosmo.dat.grid360.P020000.box3600',
            unpack = True, 
            usecols = l_cols
            ) 
    
    # Upweighted P_l(k)
    k, Pk_upw= np.loadtxt(
            'POWER_Q_CutskyN1.fidcosmo.fibcoll.dat.grid360.P020000.box3600',
            unpack = True, 
            usecols = l_cols
            ) 

    corrdelP = pickle.load(open('delP'+str(l)+'k_corr.p', 'rb'))
    uncorrdelP = tophat.delP_uncorr(k, l, fs=fs, rc=rc)

    # FC P_l(k)
    fc_Pk = Pk + uncorrdelP + corrdelP
    
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    sub.plot(k, fc_Pk, color=pretty_colors[3], label='FC TopHat')
    sub.plot(k, Pk_upw, color=pretty_colors[5], label='Upweighted (data)')
    sub.plot(k, Pk, color=pretty_colors[1], label='True (data)')

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)

    sub.set_yscale('log')
    sub.set_ylabel(r"$\mathtt{P_{"+str(l)+"}(k)}$", fontsize=30)

    sub.legend(loc='upper left')
    fig.savefig('qaplot_P_'+str(l)+'_tophat.png', bbox_inches="tight")
    plt.close()

def plot_delP_integrand(l, k_value=0.3, n_mocks=20): 
    '''

    plot f_l,l'

    '''

    for i_mock in xrange(1, n_mocks+1): 
        
        q, P0q_i, P2q_i, P4q_i = np.loadtxt(
                ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/', 
                    'POWER_Q_CutskyN'+str(i_mock)+'.fidcosmo.dat.grid360.P020000.box3600'
                    ]),
                unpack = True, 
                usecols = [0,1,2,4]
                ) 

        if i_mock == 1: 
            P0q_sum = P0q_i
            P2q_sum = P2q_i
            P4q_sum = P4q_i
        else: 
            P0q_sum += P0q_i
            P2q_sum += P2q_i
            P4q_sum += P4q_i

    P0q = P0q_sum/np.float(n_mocks)
    P2q = P2q_sum/np.float(n_mocks)
    P4q = P4q_sum/np.float(n_mocks)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for lp in [0,2,4]:
        #print np.array([tophat.f_l_lp(q_i*0.4, k_value*0.4, l, lp) for q_i in q]) - \
        #        np.array([tophat.f_l_lp(q_i*0.4, k_value*0.4, l, 4) for q_i in q])
        
        if lp == 0: 
            Pq_interp = interp1d(q, P0q, kind='cubic')
            extrap_params = [3345.0, -1.6]  
        elif lp == 2: 
            Pq_interp = interp1d(q, P2q, kind='cubic')
            extrap_params = [400.0, -4.] 
        elif lp == 4: 
            Pq_interp = interp1d(q, P4q, kind='cubic')
            extrap_params = [260.0, -1.0]

        Pk_kplus_extrap = lambda q_or_k: tophat.pk_powerlaw(q_or_k, extrap_params)

        sub.plot(
                q, 
                [
                    tophat.delPq_integrand(
                        q_i, k_value, l, lp, 
                        Pq_interp, Pk_kplus_extrap, 
                        rc=0.4, k_min=q[0], k_max=q[-1], 
                        first_order=True, 
                        w1d=False
                        )
                    for q_i in q], 
                c=pretty_colors[lp], 
                lw=2, 
                ls='-',
                label = r"$\mathtt{l' = "+str(lp)+"}$"
                )
        sub.scatter(
                q, 
                [q_i * Pq_interp(q_i) * tophat.f_l_lp(q_i*0.4, k_value*0.4, l, lp, first_order=True) for q_i in q], 
                c=pretty_colors[lp], 
                s=10,
                label = r"$\mathtt{l' = "+str(lp)+"}$ First order"
                )

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{q}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{f_{l,l'} (q r_{fc}, "+str(round(k_value,1))+"r_{fc}})$", fontsize=30)

    sub.legend(loc='upper left', scatterpoints=1)

    fig_file = ''.join(['qaplot_fllp_l', str(l), '_k', str(k_value), '.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    #plt.show()
    plt.close()

def plot_fllp(l, k_value=0.3, rc=0.4): 
    '''

    plot f_l,l'

    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for lp in [0,2,4]:
        q = np.arange( 0.002, 1.0, 0.005 ) 
        
        start_time = time.time()
        fllp_integrated = [tophat.f_l_lp(q_i*rc, k_value*rc, l, lp, first_order=True) for q_i in q] 
        print 'Integrated f_l,lp takes ', time.time() - start_time, ' seconds'

        start_time = time.time()
        fllp_estimated = [tophat.f_l_lp_est(q_i*rc, k_value*rc, l, lp) for q_i in q]
        print 'Estimated f_l,lp takes ', time.time() - start_time, ' seconds'

        sub.scatter(
                q, 
                fllp_integrated, 
                c=pretty_colors[lp], 
                s=10,
                label = r"$\mathtt{l' = "+str(lp)+"}$ First order"
                )
        sub.plot(
                q, 
                fllp_estimated,
                c=pretty_colors[lp], 
                lw = 4, 
                ls = '--', 
                label = r"theoretical estimate"
                )

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{q}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{f_{l,l'} (q r_{fc}, "+str(round(k_value,1))+"r_{fc}})$", fontsize=30)

    sub.legend(loc='upper left', scatterpoints=1)

    fig_file = ''.join(['qaplot_fl', str(l), 'lp024_k', str(k_value), '.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    plt.close()

if __name__=="__main__": 
    #for k in [0.0025]:#, 0.05, 0.1, 0.3]: 
    #    plot_fllp(0, k_value=k, rc=0.4)
    #    plot_fllp(2, k_value=k, rc=0.4)
    #    plot_fllp(4, k_value=k, rc=0.4)

    plot_delP_extrapolation_test(0, fs=1.0, rc=0.4, n_mocks=20)
    plot_delP_extrapolation_test(2, fs=1.0, rc=0.4, n_mocks=20)
    plot_delP_extrapolation_test(4, fs=1.0, rc=0.4, n_mocks=20)
    plot_delP_extrapolation_test(0, type = 'difference', fs=1.0, rc=0.4, n_mocks=20)
    plot_delP_extrapolation_test(2, type = 'difference', fs=1.0, rc=0.4, n_mocks=20)
    plot_delP_extrapolation_test(4, type = 'difference', fs=1.0, rc=0.4, n_mocks=20)
    #plot_delP_uncorr([0,2])
    #plot_fllp(0, k_value=0.3)
    #for l in [0,2]:
    #    for k in [0.0025, 0.05, 0.1, 0.3]: 
    #        plot_fllp(l, k_value=k)

"""
def plot_delP_uncorr(l, fs=1.0, rc=0.4): 
    '''

    '''
    if not isinstance(l, list): 
        l = [l]

    # import k, Pk data  
    k, P0k, P2k = np.loadtxt(
            'POWER_Q_CutskyN1.fidcosmo.dat.grid360.P020000.box3600',
            unpack = True, 
            usecols = [0,1,2]
            ) 

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for i_l, l_i in enumerate(l):
        delP = tophat.delP_uncorr(k, l_i, fs=fs, rc=rc)

        sub.plot(
                k, 
                delP, 
                c=pretty_colors[i_l+1], 
                lw=4, 
                label=r"$\Delta P_"+str(l_i)+"^{uncorr}(k)$"
                )

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{\Delta P_{i}^{uncorr}(k)}$", fontsize=30)

    sub.legend(loc='lower left')
    fig.savefig('qaplot_delP_uncorr.png', bbox_inches="tight")
    plt.close()

def plot_delP_corr(l): 
    '''

    '''
    if not isinstance(l, list): 
        l = [l]

    # import k, Pk data  
    k, P0k, P2k = np.loadtxt(
            'POWER_Q_CutskyN1.fidcosmo.dat.grid360.P020000.box3600',
            unpack = True, 
            usecols = [0,1,2]
            ) 

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for i_l, l_i in enumerate(l):
        #delP = tophat.delP_uncorr(k, l_i, fs=fs, rc=rc)

        delP = pickle.load(open('delP'+str(l_i)+'k_corr.p', 'rb'))

        sub.plot(
                k, 
                delP, 
                c=pretty_colors[i_l+1], 
                lw=4, 
                label=r"$\Delta P_"+str(l_i)+"^{corr}(k)$"
                )

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{\Delta P_{i}^{corr}(k)}$", fontsize=30)

    sub.legend(loc='upper left')
    fig.savefig('qaplot_delP_corr.png', bbox_inches="tight")
    plt.close()

def plot_fllp_diagonal(l, k_value=0.3, rc=0.4): 
    '''

    plot f_l,l'

    '''
    # import q values
    q = np.loadtxt(
            ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/', 
                'POWER_Q_CutskyN1.fidcosmo.dat.grid360.P020000.box3600'
                ]),
            unpack = True, 
            usecols = [0]
            ) 
    q = np.arange(0.002, 1.0, 0.005)
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    lp = l  # l' = l
        
    sub.scatter(
            q, 
            [tophat.f_l_lp(q_i*rc, k_value*rc, l, lp, first_order=True) for q_i in q], 
            c=pretty_colors[3], 
            s=10,
            label = r"$\mathtt{f_{l = "+str(l)+",l' = "+str(lp)+"}}$ First order"
            )
    
    def fllp_diagonal(q_in): 
        
        if q_in < k_value: 
            return tophat.W_2d(q_in * rc) * (q_in / k_value)**(l+1)
        else: 
            return tophat.W_2d(q_in * rc) * (k_value / q_in)**l

    sub.plot(
            q, 
            [fllp_diagonal(q_i) for q_i in q],
            c=pretty_colors[5], 
            lw = 4, 
            ls = '--', 
            label = r"theoretical estimate"
            )

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{q}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{f_{l,l'} (q r_{fc}, "+str(round(k_value,1))+"r_{fc}})$", fontsize=30)

    sub.legend(loc='upper left', scatterpoints=1)

    fig_file = ''.join(["qaplot_f_l", str(l), "lp", str(l), '_k', str(k_value), '.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    plt.close()

def plot_fllp_offdiagonal(l, k_value=0.3, rc=0.4): 
    '''

    plot f_l,l'

    '''
    # import q values
    q = np.loadtxt(
            ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/', 
                'POWER_Q_CutskyN1.fidcosmo.dat.grid360.P020000.box3600'
                ]),
            unpack = True, 
            usecols = [0]
            ) 
    q = np.arange(0.002, 1.0, 0.005)
    lps = [0,2,4]
    lps.pop(np.where(np.array(lps) == l)[0])
    print lps
    
    def fllp_offdiagonal(q_in, l_in, lp_in): 
        
        if l_in > lp_in: 
            if q_in < k_value: 
                return (2.0 * l_in + 1.)/2.0 * (q_in / k_value) * fllp_poly(l_in, lp_in, q_in / k_value) * tophat.W_2d(q_in * rc)
            else: 
                return 0.0
        else: 
            if q_in > k_value: 
                return (2.0 * l_in + 1.)/2.0 * fllp_poly(lp_in, l_in, k_value / q_in) * tophat.W_2d(q_in * rc)
            else: 
                return 0.0 
        
    for i_lp, lp in enumerate(lps): 
        prettyplot()
        pretty_colors = prettycolors()

        fig = plt.figure(figsize=(10,10))
        sub = fig.add_subplot(111)
    
        sub.scatter(
                q, 
                [tophat.f_l_lp(q_i*rc, k_value*rc, l, lp, first_order=True) for q_i in q], 
                c=pretty_colors[3], 
                s=10,
                label = r"$\mathtt{f_{l = "+str(l)+",l' = "+str(lp)+"}}$ First order"
                )

        sub.plot(
                q, 
                [fllp_offdiagonal(q_i, l, lp) for q_i in q],
                c=pretty_colors[5], 
                lw = 4, 
                ls = '--', 
                label = r"theoretical estimate"
                )

        sub.set_xlim([10**-3,10**0])
        sub.set_xscale("log") 
        sub.set_xlabel(r"$\mathtt{q}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
        sub.set_ylabel(r"$\mathtt{f_{l,l'} (q r_{fc}, "+str(round(k_value,1))+"r_{fc}})$", fontsize=30)

        sub.legend(loc='upper left', scatterpoints=1)

        fig_file = ''.join(["qaplot_f_l", str(l), "lp", str(lp), '_k', str(k_value), '.png'])
        fig.savefig(fig_file, bbox_inches='tight')
        plt.close()

def fllp_poly(l1_in, l2_in, x):  
    
    if (l1_in == 2) and (l2_in == 0): 
        return x**2 - 1.
    elif (l1_in == 4) and (l2_in == 0): 
        return (7./4.)*x**4 - (5./2.)*x**2 + 3./4.
    elif (l1_in == 4) and (l2_in == 2): 
        return x**4 - x**2
    else: 
        raise ValueError


"""
