'''

Plots for xi tohat fiber collision correction methods 

'''

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
    
    corrdelP = pickle.load(open('delP'+str(l)+'k_corr_downsampled.p', 'rb'))
    #corrdelP_noextrap = pickle.load(open('delP'+str(l)+'k_corr_noextrap.p', 'rb'))
    #corrdelP_noextrap_firstorder = pickle.load(open('delP'+str(l)+'k_corr_noextrap_firstorder.p', 'rb'))
    uncorrdelP = tophat.delP_uncorr(k, l, fs=fs, rc=rc)

    sub.plot(
            k, 
            uncorrdelP/(2.0 * np.pi)**0., 
            c=pretty_colors[3], 
            lw=4, 
            label="Uncorrelated"
            )
    nonzero = np.where(corrdelP != 0.0) 
    sub.plot(
            k[nonzero], 
            corrdelP[nonzero] * (2.0 * np.pi)**-3., 
            c = pretty_colors[1], 
            lw = 4, 
            label = "Correlated"
            )
    #sub.scatter(
    #        k, 
    #        corrdelP_noextrap/(2.0 * np.pi)**0., 
    #        c = pretty_colors[2], 
    #        label = "Correlated (No extrap.)"
    #        )
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
        sub.legend(loc='upper left', scatterpoints = 1)
    fig.savefig('qaplot_delP_'+str(l)+'.png', bbox_inches="tight")
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

def plot_fllp(l, k_value=0.3): 
    '''

    plot f_l,l'

    '''
    # import k values
    q, P0q, P2q, P4q = np.loadtxt(
            'POWER_Q_CutskyN1.fidcosmo.dat.grid360.P020000.box3600',
            unpack = True, 
            usecols = [0,1,2,4]
            ) 
    

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for lp in [0,2,4]:
        #print np.array([tophat.f_l_lp(q_i*0.4, k_value*0.4, l, lp) for q_i in q]) - \
        #        np.array([tophat.f_l_lp(q_i*0.4, k_value*0.4, l, 4) for q_i in q])
        
        if lp == 0: 
            Pq_interp = interp1d(q, P0q, kind='cubic')
        elif lp == 2: 
            Pq_interp = interp1d(q, P2q, kind='cubic')
        elif lp == 4: 
            Pq_interp = interp1d(q, P4q, kind='cubic')

        sub.plot(
                q, 
                [q_i * Pq_interp(q_i) * tophat.f_l_lp(q_i*0.4, k_value*0.4, l, lp) for q_i in q], 
                c=pretty_colors[lp], 
                lw=4, 
                ls='--',
                label = r"$\mathtt{l' = "+str(lp)+"}$"
                )

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{q}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{f_{l,l'} (q r_{fc}, "+str(round(k_value,1))+"r_{fc}})$", fontsize=30)

    sub.legend(loc='upper left')
    plt.show()
    plt.close()

if __name__=="__main__": 
    plot_delP(0, fs=1.0, rc=0.4, n_mocks=20)
    plot_delP(2, fs=1.0, rc=0.4, n_mocks=20)
    #plot_delP_uncorr([0,2])
    #plot_fllp(0, k_value=0.3)
    #plot_fllp(2, k_value=0.0025)
    #plot_fllp(2, k_value=0.05)
    #plot_fllp(2, k_value=0.1)
    #plot_fllp(2, k_value=0.3)
    #plot_corrected_Pk(2)

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

"""
