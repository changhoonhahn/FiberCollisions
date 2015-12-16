'''

Plots for xi tohat fiber collision correction methods 

'''
import time
import pickle
import numpy as np 
import scipy as sp

from scipy.interpolate import interp1d

# local ----
from fourier_corr import fourier_corr as tophat
from corr_spec.corr_average import CorrAvgSpec as AvgSpec

# plotting ----
from defutility.plotting import prettyplot 
from defutility.plotting import prettycolors 
"""
def plot_delP(l, fs=1.0, rc=0.4, n_mocks=84, Ngrid=360, extrap_params=[[3345.0, -1.6], [400.0, -4.]], k_fixed=None, **kwargs):
    '''
    Comparison of delP^corr, delP^uncorr, and P^upw_avg - P^true_avg

    '''
    # average P_l(k) and P_l^upw(k)
    true_cat_corr = {
            'catalog': {'name': 'nseries'}, 
            'correction': {'name': 'true'}, 
            'spec': {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid}
            }
    upw_cat_corr = {
            'catalog': {'name': 'nseries'}, 
            'correction': {'name': 'upweighted'}, 
            'spec': {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid}
            }

    k, Pks = avgSpec(n_mocks, 'pallk', true_cat_corr)
    k_upw, Pk_upws = avgSpec(n_mocks, 'pallk', upw_cat_corr)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    # uncorrelated delP
    #uncorrdelP = tophat.delP_uncorr(k, l, fs=fs, rc=rc)
    
    sub.plot(
            k, 
            uncorrdelP, 
            c=pretty_colors[3], 
            lw=4, 
            label="Uncorrelated"
            )
    
    # correlated delP 
    #corrdelP = tophat.delP_corr(k, Pks, l, fs=fs, rc=rc, extrap_params=extrap_params, k_fixed)

    sub.plot(
            k, 
            corrdelP, 
            c = pretty_colors[1], 
            lw = 4, 
            label = "Correlated"
            )
    
    # combined uncorrelated and correlated delP
    sub.plot(
            k, 
            uncorrdelP + corrdelP, 
            c='k', 
            lw=2, 
            label="Combined"
            )

    # delP from data
    if l == 0: 
        l_index = 0
    elif l == 2: 
        l_index = 1
    elif l == 4: 
        l_index = 2
    sub.plot(
            k, 
            Pk_upws[l_index] - Pks[l_index], 
            c = 'k', 
            lw = 4,
            ls = '--', 
            label = r"$\mathtt{P^{upw}(k) - P^{true}(k)}$"
            )
    
    if 'xrange' in kwargs.keys(): 
        sub.set_xlim(kwargs['xrange'])
    else:
        sub.set_xlim([10**-3,10**0])

    if 'yrange' in kwargs.keys(): 
        sub.set_ylim(kwargs['yrange'])

    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(l)+"}(k)}$", fontsize=30)
    
    if l == 0: 
        sub.legend(loc='lower right', scatterpoints = 1)
    elif l == 2: 
        sub.legend(loc='upper right', scatterpoints = 1)
    elif l == 4: 
        sub.legend(loc='lower right', scatterpoints = 1)

    fig_file = ''.join([
        'figure/',
        'plot_delP_', str(l), '_k_fixed', str(round(k_fixed,2)),'_Ngrid', str(Ngrid), '.png'
        ])
    fig.savefig(fig_file, bbox_inches="tight")
    plt.show()
    plt.close()
"""

def plot_delP_lp_component(ell, mock='nseries', n_mocks=84, Ngrid=360):
    '''
    Compare each of the l' components of delP_l^corr. l' components are read from pickle files.
    '''
    # average P_l(k) and P_l^upw(k)
    specdict = {
            'P0': 20000,
            'Lbox': 3600, 
            'Ngrid': Ngrid, 
            'ell': ell 
            }
    true_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'true'},
            'spec': specdict
            }
    true_spec = AvgSpec(n_mocks, 'pk', true_cat_corr)
    true_spec.read()
    k = true_spec.k
    Pk = getattr(true_spec, 'p'+str(ell)+'k')
    upw_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'upweight'},
            'spec': specdict
            }
    upw_spec = AvgSpec(n_mocks, 'pk', upw_cat_corr)
    upw_spec.read()
    k_upw = upw_spec.k
    Pk_upw = getattr(upw_spec, 'p'+str(ell)+'k')

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    for i_mock in range(1,n_mocks+1): 
        if mock == 'nseries': # this needs to be fixed
            pass
            #if Ngrid == 960:
            #    corrdelP = pickle.load(open('delP'+str(l)+'k_corr_k_fixed0.6_kmax0.8_Ngrid'+str(Ngrid)+'.p', 'rb'))
            #elif Ngrid == 720: 
            #    corrdelP = pickle.load(open('delP'+str(l)+'k_corr_estimated_k_fixed0.6_Ngrid'+str(Ngrid)+'.p', 'rb'))
        elif mock == 'nseriesbox': 
            lps = [0, 2, 4, 6]
    
        for i_lp, lp in enumerate(lps):
            pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                'corrdelP2k_lp', str(lp), '_power3600z_BoxN', str(i_mock), '.fourier_tophat.fs1.0.rc0.43.kfit4.3.kfixed4.34.dat.p'
                ])
            k_lpcomp, corrdelP_lpcomp = pickle.load(open(pickle_file, 'rb'))

            sub.plot(
                k_lpcomp, corrdelP_lpcomp,
                c = pretty_colors[i_lp+1],
                lw = 3, 
                ls = '--',
                label = "Correlated: "+r"$l' ="+str(lp)+"$"
                )
            
            try: 
                corrdelP += corrdelP_lpcomp
            except UnboundLocalError: 
                corrdelP = corrdelP_lpcomp

    sub.plot(
            k_lpcomp, 
            corrdelP, 
            c= 'gray',
            lw=2, 
            label='Total'
            )

    sub.plot(
            k, 
            Pk_upw - Pk, 
            c= 'k',
            lw=2, 
            label='data'
            )

    sub.set_xlim([10**-3,10**1])
    if ell == 2: 
        sub.set_ylim([-50.0, 250.])
    elif ell == 4: 
        sub.set_ylim([-50.0, 1000.])

    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(ell)+"}^{corr}(k)}$", fontsize=30)
    
    if ell == 0: 
        sub.legend(loc='lower left', scatterpoints = 1)
    elif ell == 2: 
        sub.legend(loc='upper left', scatterpoints = 1)
    elif ell == 4: 
        sub.legend(loc='lower left', scatterpoints = 1)
     
    #fig.savefig(
    #        'figure/', 
    #        'qaplot_delP_'+str(l)+'corr_lp_components_', str(n_mocks), mock, 'mocks.png', 
    #       bbox_inches="tight")
    plt.show()
    plt.close()


def plot_delP_lp_comp_nseries(ell, n_mocks=1, Ngrid=960):
    '''
    Compare each of the l' components of delP_l^corr. l' components are read from pickle files.
    '''
    # average P_l(k) and P_l^upw(k)
    specdict = {
            'P0': 20000,
            'Lbox': 3600, 
            'Ngrid': Ngrid, 
            'ell': ell 
            }
    true_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'true'},
            'spec': specdict
            }
    true_spec = AvgSpec(20, 'pk', true_cat_corr)
    true_spec.read()
    k = true_spec.k
    Pk = getattr(true_spec, 'p'+str(ell)+'k')
    upw_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'upweight'},
            'spec': specdict
            }
    upw_spec = AvgSpec(20, 'pk', upw_cat_corr)
    upw_spec.read()
    k_upw = upw_spec.k
    Pk_upw = getattr(upw_spec, 'p'+str(ell)+'k')

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(14,8))
    sub = fig.add_subplot(111)
    
    for mock in ['nseriesbox', 'nseries']:
        if mock == 'nseries': # this needs to be fixed
            lps = [0, 2, 4]
            lstyle = '--'
        elif mock == 'nseriesbox': 
            lps = [0, 2, 4, 6]
            lstyle = '-'

        for i_lp, lp in enumerate(lps):
            for i_mock in range(1,n_mocks+1): 

                if mock == 'nseriesbox': 
                    pickle_file = ''.join([
                        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                        'corrdelP', str(ell), 'k_lp', str(lp), '_power3600z_BoxN', str(i_mock), 
                        '.fourier_tophat.fs1.0.rc0.43.kfit4.3.kfixed4.34.dat.p'
                        ])
                    lstyle = '-'
                    label = "Correlated: "+r"$l' ="+str(lp)+"$"

                elif mock == 'nseries': 
                    if ell == 0: 
                        q_str = ''
                    else: 
                        q_str = 'Q_'
                    pickle_file = ''.join([
                        '/mount/riachuelo1/hahn/power/Nseries/', 
                        'corrdelP', str(ell), 'k_lp', str(lp), '_POWER_', q_str, 'CutskyN', str(i_mock), 
                        '.fidcosmo.fourier_tophat.fs1.0.rc0.43.kfit0.7.kfixed0.84.dat.grid960.P020000.box3600.p'
                        ])
                    label = None

                k_lpcomp, corrdelP_lpcomp = pickle.load(open(pickle_file, 'rb'))
                if mock == 'nseriesbox': 
                    print 'lp = ', lp
                    print corrdelP_lpcomp.min(), corrdelP_lpcomp.max()
                
                if i_mock == 1: 
                    avg_corrdelP_lpcomp = corrdelP_lpcomp
                else: 
                    avg_corrdelP_lpcomp += corrdelP_lpcomp
            
            avg_corrdelP_lpcomp /= np.float(n_mocks)

            sub.plot(
                k_lpcomp, avg_corrdelP_lpcomp,
                c = pretty_colors[i_lp+1],
                lw = 3, 
                ls = lstyle,
                label = label
                )
            
            try: 
                corrdelP += avg_corrdelP_lpcomp
            except UnboundLocalError: 
                corrdelP = avg_corrdelP_lpcomp

        # del P^uncorr
        if mock == 'nseries': 
            uncorrdelPk_pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/',
                'uncorrdelP', str(ell), 
                'k_POWER_CutskyN1.fidcosmo.fourier_tophat.fs1.0.rc0.43.kfit0.7.kfixed0.84.dat.grid960.P020000.box3600.p'])
        
        elif mock == 'nseriesbox': 
            uncorrdelPk_pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                'uncorrdelP', str(ell), 
                'k_power3600z_BoxN1.fourier_tophat.fs1.0.rc0.43.kfit4.3.kfixed4.34.dat.p'
                ])

        k, uncorrdelP = pickle.load(open(uncorrdelPk_pickle_file, 'rb'))

        if mock == 'nseriesbox': 
            uncorr_label = 'Uncorrelated'
            tot_label = 'Total'
            tot_lstyle = '-'
        else: 
            uncorr_label = None
            tot_label = None
            tot_lstyle = '--'
        sub.plot(k, uncorrdelP, c= 'red', label = uncorr_label, ls= tot_lstyle, lw=2)
        sub.plot(k_lpcomp, corrdelP + uncorrdelP, 
                c= 'gray',
                ls= tot_lstyle,
                lw=2, 
                label=tot_label
                )
        
        del corrdelP

    sub.plot(
            k, 
            Pk_upw - Pk, 
            c= 'k',
            lw=2, 
            label='data'
            )

    sub.set_xlim([10**-3,10**1])
    if ell == 0: 
        sub.set_ylim([-1000., 100.])
    elif ell == 2: 
        sub.set_ylim([-50.0, 250.])
    elif ell == 4: 
        sub.set_ylim([-50.0, 1000.])

    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(ell)+"}^{corr}(k)}$", fontsize=30)
    
    if ell == 0: 
        sub.legend(loc='lower right', scatterpoints = 1)
    elif ell == 2: 
        sub.legend(loc='upper left', scatterpoints = 1)
    elif ell == 4: 
        sub.legend(loc='lower left', scatterpoints = 1)
     
    fig.savefig(
            ''.join([
                'figure/', 
                'qaplot_delP_'+str(ell)+'corr_lp_components_', str(n_mocks), 'nseriesbox.png']), 
            bbox_inches="tight")
    #plt.show()
    plt.close()


def plot_delP_extrapolation_test(l, type='normal', fs=1.0, rc=0.4, n_mocks=84, Ngrid=360):
    '''
    Comparison of delP^corr, delP^uncorr, and P^upw_avg - P^true_avg

    '''
    k, Pk = pk_extrap.average_Pk(l, n_mocks, Ngrid=Ngrid)
    k_upw, Pk_upw = pk_extrap.average_Pk_upw(l, n_mocks, Ngrid=Ngrid)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    corrdelP = pickle.load(open('delP'+str(l)+'k_corr_estimated_k_fixed0.6_Ngrid'+str(Ngrid)+'.p', 'rb'))
    corrdelP_noextrap = pickle.load(open('delP'+str(l)+'k_corr_estimated_noextrap_k_fixed0.6_Ngrid'+str(Ngrid)+'.p', 'rb'))

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

    elif type == 'difference': 
        sub.plot(
                k, 
                corrdelP - corrdelP_noextrap, 
                c = pretty_colors[1], 
                lw = 4, 
                )

    sub.set_xlim([10**-3,10**0])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    
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
        fig.savefig('qaplot_delPcorr_'+str(l)+'_extrapolationtest_k_fixed0.6_Ngrid'+str(Ngrid)+'.png', bbox_inches="tight")
    if type == 'difference':
        fig.savefig('qaplot_delPcorr_'+str(l)+'_extrapolationtest_diff_k_fixed0.6_Ngrid'+str(Ngrid)+'.png', bbox_inches="tight")
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

def plot_delP_integrand(l, k_value=0.3, k_fixed=0.6, n_mocks=20, Ngrid=360): 
    '''

    plot f_l,l'

    '''

    for i_mock in xrange(1, n_mocks+1): 
        
        q, P0q_i, P2q_i, P4q_i = np.loadtxt(
                ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/', 
                    'POWER_Q_CutskyN'+str(i_mock)+'.fidcosmo.dat.grid', 
                    str(Ngrid), '.P020000.box3600'
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

        Pk_kplus_extrap = lambda q_or_k: pk_extrap.pk_powerlaw(q_or_k, extrap_params, k_fixed=k_fixed)

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

def plot_delP_corr_extrapolations(l, n_mocks=10, k_fixed=0.6, k_max=0.5, Ngrid=360, fs=1.0, rc=0.4, **kwargs):
    '''
    Comparison of delP^corr, delP^uncorr, and P^upw_avg - P^true_avg

    '''

    if isinstance(k_max, list) or isinstance(k_max, np.ndarray):  
        pass
    else: 
        k_max = [k_max]
    
    # average P_l(k) and P_l^upw(k)
    k, Pk = pk_extrap.average_Pk(l, n_mocks, Ngrid=Ngrid)
    k_upw, Pk_upw = pk_extrap.average_Pk_upw(l, n_mocks, Ngrid=Ngrid)
    
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    # uncorrelated delP
    uncorrdelP = tophat.delP_uncorr(k, l, fs=fs, rc=rc)
    
    sub.plot(
            k, 
            uncorrdelP, 
            c=pretty_colors[3], 
            lw=4, 
            label="Uncorrelated"
            )

    # correlated delP 
    for k_max_i in k_max: 

        pickle_file = ''.join([
            'delP', str(l), 'k_corr', 
            '_k_fixed', str(round(k_fixed, 1)),
            '_kmax', str(round(k_max_i, 2)), 
            '_Ngrid', str(Ngrid), '.p'
            ])

        corrdelP = pickle.load(open(pickle_file, 'rb'))

        sub.plot(
                k, 
                corrdelP, 
                lw = 4, 
                ls = '--', 
                label = "Correlated "+r"$\mathtt{k_{fit} ="+str(round(k_max_i, 2))+"}$"
                )

    # delP from data
    sub.plot(
            k, 
            Pk_upw - Pk, 
            c = 'k', 
            lw = 4,
            ls = '--', 
            label = r"$\mathtt{P^{upw}(k) - P^{true}(k)}$"
            )
    
    if 'xrange' in kwargs.keys(): 
        sub.set_xlim(kwargs['xrange'])
    else:
        sub.set_xlim([10**-3,10**0])

    if 'yrange' in kwargs.keys(): 
        sub.set_ylim(kwargs['yrange'])

    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(l)+"}(k)}$", fontsize=30)
    
    if l == 0: 
        sub.legend(loc='lower right', scatterpoints = 1)
    elif l == 2: 
        sub.legend(loc='upper right', scatterpoints = 1)
    elif l == 4: 
        sub.legend(loc='lower right', scatterpoints = 1)

    fig_file = ''.join([
        'qaplot_delP_', str(l), '_k_fixed0.6_Ngrid', str(Ngrid), '_extrapolations.png'
        ])
    fig.savefig(fig_file, bbox_inches="tight")
    plt.show()
    plt.close()

if __name__=="__main__": 
    #for k in [0.0025]:#, 0.05, 0.1, 0.3]: 
    #    plot_fllp(0, k_value=k, rc=0.4)
    #    plot_fllp(2, k_value=k, rc=0.4)
    #    plot_fllp(4, k_value=k, rc=0.4)

    #plot_delP(0, fs=1.0, rc=0.4, n_mocks=20)
    #plot_delP(2, fs=1.0, rc=0.4, n_mocks=20)
    #plot_delP(4, fs=1.0, rc=0.4, n_mocks=20)

    #plot_delP(2, n_mocks=10, Ngrid=960)#, xrange=[0.1, 1.0], yrange=[-100., 1000.])
    #plot_delP(2, n_mocks=10, Ngrid=720)
    #for l_i in [0,2,4]:
    #    plot_delP(l_i, n_mocks=10, Ngrid=720)
    plot_delP_lp_comp_nseries(0, n_mocks=4, Ngrid=960)
    #plot_delP_lp_comp_nseries(2, n_mocks=3, Ngrid=960)
    #plot_delP_lp_component(2, mock='nseriesbox', n_mocks=1, Ngrid=960)
    #plot_delP_lp_component(4, n_mocks=10, Ngrid=720)
    #plot_delP_lp_component(2, n_mocks=10, Ngrid=960)
    #plot_delP_lp_component(4, n_mocks=10, Ngrid=960)
        #plot_delP_extrapolation_test(l_i, type='normal', n_mocks=10, Ngrid=720)
        #plot_delP_extrapolation_test(l_i, type='difference', n_mocks=10, Ngrid=720)

    #plot_delP_uncorr([0,2])
    #plot_fllp(0, k_value=0.3)
    #for l in [0,2]:
    #    for k in [0.0025, 0.05, 0.1, 0.3]: 
    #        plot_fllp(l, k_value=k)
