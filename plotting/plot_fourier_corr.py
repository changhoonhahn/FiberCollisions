'''

Plots for xi tohat fiber collision correction methods 

'''
import time
import pickle
import numpy as np 
import scipy as sp

from scipy.interpolate import interp1d

# local ----
from corr_spec.corr_average import CorrAvgSpec as AvgSpec
from fourier_corr.test_fourier_corr import delPcorr_ktrust_inf_bestfit

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

def plot_delP_lp_comp_nseries(ell, Ngrid=960, noextrap=''):
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

    fig = plt.figure(figsize=(17,8))
    sub = fig.add_subplot(111)
    
    for mock in ['nseriesbox']:#, 'nseries']:
        if mock == 'nseries': # this needs to be fixed
            lps = [0, 2, 4]
            lstyle = '--'
        elif mock == 'nseriesbox': 
            lps = [0, 2, 4, 6, 8, 10]
            lstyle = '-'

        for i_lp, lp in enumerate(lps):

            if mock == 'nseriesbox': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                    'corrdelP', str(ell), 'k_lp', str(lp), noextrap, '_AVG_power3600z_BoxN.dat.p'
                    ])
                lstyle = '-'
                label = "Correlated: "+r"$l' ="+str(lp)+"$"

            elif mock == 'nseries': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/', 
                    'corrdelP', str(ell), 'k_lp', str(lp), noextrap, '_AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600.p'
                    ])
                label = None

            k_lpcomp, corrdelP_lpcomp = pickle.load(open(pickle_file, 'rb'))
            if mock == 'nseriesbox': 
                print 'lp = ', lp
                print pickle_file 
                print corrdelP_lpcomp.min(), corrdelP_lpcomp.max()
            
            sub.plot(
                k_lpcomp, corrdelP_lpcomp,
                c = pretty_colors[i_lp+1],
                lw = 2, 
                ls = lstyle,
                label = label
                )
                
            if i_lp == 0: 
                corrdelP = corrdelP_lpcomp
            else: 
                corrdelP += corrdelP_lpcomp

        # del P^uncorr
        if mock == 'nseries': 
            uncorrdelPk_pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/',
                'uncorrdelP', str(ell), 'k_AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600.p'
                ])
        elif mock == 'nseriesbox': 
            uncorrdelPk_pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                'uncorrdelP', str(ell), 'k_AVG_power3600z_BoxN.dat.p'
                ])

        uncorr_k, uncorrdelP = pickle.load(open(uncorrdelPk_pickle_file, 'rb'))

        if mock == 'nseriesbox': 
            uncorr_label = 'Uncorrelated'
            uncorr_lstyle = '-.'
            tot_label = 'Total'
            tot_lstyle = '-'
        else: 
            uncorr_label = None
            tot_label = None
            tot_lstyle = '--'

        sub.plot(uncorr_k, uncorrdelP, c= 'gray', label = uncorr_label, ls= uncorr_lstyle, lw=2)
        sub.plot(k_lpcomp, corrdelP + uncorrdelP, 
                c= 'black',
                ls= tot_lstyle,
                lw=4, 
                label=tot_label
                )
        
        del corrdelP

    sub.scatter(
            k, 
            Pk_upw - Pk, 
            c= 'k',
            label='data'
            )
    #        lw=2, 
    #        ls = '--', 
    #        )

    if ell == 0: 
        sub.set_ylim([-1000., 100.])
        sub.set_xlim([10**-3,10**1])
    elif ell == 2: 
        sub.set_ylim([-50.0, 350.])
        sub.set_xlim([10**-4,10**1])
    else: 
        raise NotImplementedError
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
                'qaplot_delP_'+str(ell)+'corr_lp_components', noextrap, '_nseriesbox.png']), 
            bbox_inches="tight")
    #plt.show()
    plt.close()

def plot_delPcorr_ktrust(ell, k_trust=0.5, Ngrid=960):
    '''
    Investigate the contribution of delPcorr from q range k_trust to infinity. 
    delPcorr is being divided into two parts:
        - 0 to k_trust, which can be reliably calculated 
        - k_trust to infinity, which cannot be reliably calculated because 
        we do not trust the models or the extrapolations beyound this point.
        However they it may be possible to have a polynomial fit. 
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

    fig = plt.figure(figsize=(17,8))
    sub = fig.add_subplot(111)
    
    # Del P 
    sub.scatter(k, Pk_upw - Pk, c= 'k', 
            label=r"Data: $\mathtt{P_"+str(ell)+"^{upw}(k) - P_"+str(ell)+"^{true}(k)}$")
    
    for mock in ['nseriesbox']:
        if mock == 'nseries': # this needs to be fixed
            lps = [0, 2, 4]
            lstyle = '--'
        elif mock == 'nseriesbox': 
            lps = [0, 2, 4, 6, 8, 10]
            lstyle = '-'

        for i_lp, lp in enumerate(lps):
            # total delP^corr integrated from 0 to infinity 
            if mock == 'nseriesbox': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                    'corrdelP', str(ell), 'k_lp', str(lp), 
                    '_AVG_power3600z_BoxN.dat.p'
                    ])
            elif mock == 'nseries': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/', 
                    'corrdelP', str(ell), 'k_lp', str(lp), 
                    '_AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600.p'
                    ])
            print pickle_file
            k_lp, corrdelP_lp = pickle.load(open(pickle_file, 'rb'))
            
            # delP^corr 0 to k_trust
            if mock == 'nseriesbox': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                    'corrdelP', str(ell), 'k_lp', str(lp), 
                    '_qmax', 
                    str(round(k_trust,2)),
                    '_AVG_power3600z_BoxN.dat.p'
                    ])
            print pickle_file

            k_lp_0_ktrust, corrdelP_lp_0_ktrust = pickle.load(open(pickle_file, 'rb'))

            corrdelP_lp_ktrust_inf = corrdelP_lp - corrdelP_lp_0_ktrust
            lp_label = ''.join([
                r"$\mathtt{l' = ", str(lp), "}; ", 
                "\mathtt{k_{trust} =", str(round(k_trust,2)), "}$ to $\infty$"
                ])
            sub.plot(k_lp, corrdelP_lp_ktrust_inf, 
                    c= pretty_colors[i_lp],
                    ls= '-',
                    lw=4, 
                    label=lp_label)

            if i_lp == 0: 
                corrdelP = corrdelP_lp
                corrdelP_ktrust_inf = corrdelP_lp_ktrust_inf
            else: 
                corrdelP += corrdelP_lp
                corrdelP_ktrust_inf += corrdelP_lp_ktrust_inf

        sub.plot(k_lp, corrdelP, 
                c= 'black',
                ls= '-.',
                lw=2, 
                label="Total, $0$ to $\infty$"
                )

        sub.plot(k_lp, corrdelP_ktrust_inf, 
                c= 'black',
                ls= '-',
                lw=4, 
                label="Total, $k_{trust} ="+str(round(k_trust,2))+"$"+" to $\infty$"
                )

    if ell == 0: 
        sub.set_ylim([-1000., 100.])
        sub.set_xlim([10**-3,10**1])
    elif ell == 2: 
        sub.set_ylim([-50.0, 350.])
        sub.set_xlim([10**-4,10**1])
    else: 
        raise NotImplementedError
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
                'qaplot_delP_'+str(ell)+'corr_lp_components', 
                '_ktrust', str(round(k_trust, 2)), '_inf_nseriesbox.png']), 
            bbox_inches="tight")
    #plt.show()
    plt.close()

def plot_delPcorr_ktrust_poly(ell, k_trust=0.5, Ngrid=960, noextrap=False):
    '''
    Investigate the contribution of delPcorr from q range k_trust to infinity. 
    delPcorr is being divided into two parts:
        - 0 to k_trust, which can be reliably calculated 
        - k_trust to infinity, which cannot be reliably calculated because 
        we do not trust the models or the extrapolations beyound this point.
        However they it may be possible to have a polynomial fit. 
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

    fig = plt.figure(figsize=(17,8))
    sub = fig.add_subplot(111)
    
    for mock in ['nseriesbox']:
        if mock == 'nseries': # this needs to be fixed
            lps = [0, 2, 4]
            lstyle = '--'
        elif mock == 'nseriesbox': 
            lps = [0, 2, 4, 6, 8, 10]
            lstyle = '-'

        for i_lp, lp in enumerate(lps):
            # total delP^corr integrated from 0 to infinity 
            if mock == 'nseriesbox': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                    'corrdelP', str(ell), 'k_lp', str(lp), 
                    '_AVG_power3600z_BoxN.dat.p'
                    ])
            elif mock == 'nseries': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/', 
                    'corrdelP', str(ell), 'k_lp', str(lp), 
                    '_AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600.p'
                    ])
            print pickle_file
            k_lp, corrdelP_lp = pickle.load(open(pickle_file, 'rb'))
            
            # delP^corr 0 to k_trust
            if mock == 'nseriesbox': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                    'corrdelP', str(ell), 'k_lp', str(lp), 
                    '_qmax', 
                    str(round(k_trust,2)),
                    '_AVG_power3600z_BoxN.dat.p'
                    ])
            print pickle_file

            k_lp_0_ktrust, corrdelP_lp_0_ktrust = pickle.load(open(pickle_file, 'rb'))

            corrdelP_lp_ktrust_inf = corrdelP_lp - corrdelP_lp_0_ktrust

            if i_lp == 0: 
                corrdelP = corrdelP_lp
                corrdelP_ktrust_inf = corrdelP_lp_ktrust_inf
            else: 
                corrdelP += corrdelP_lp
                corrdelP_ktrust_inf += corrdelP_lp_ktrust_inf

        sub.plot(k_lp, corrdelP, 
                c= 'black',
                ls= '-.',
                lw=2, 
                label="Total, $0$ to $\infty$"
                )

        sub.plot(k_lp, corrdelP_ktrust_inf, 
                c= 'black',
                ls= '-',
                lw=4, 
                label="Total, $k_{trust} ="+str(round(k_trust,2))+"$"+" to $\infty$"
                )

        if noextrap: 
            extrap_str = '.noextrap'
        else: 
            extrap_str = ''

        coeff_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/Box/',
            'k^n_coeff',
            '.ell', str(ell), 
            '.ktrust', str(round(k_trust,2)), 
            extrap_str, 
            '.dat'])
        n_ks, n_ks_coeff = np.loadtxt(coeff_file, skiprows=1, unpack=True, usecols=[0, 1])

        delp_poly = lambda kk: \
                np.sum([ n_ks_coeff[ii]*kk**n_ks[ii] for ii in range(len(n_ks)) ])

        for i_n_ks in xrange(len(n_ks)): 
            sub.plot(
                    k_lp, 
                    np.array(n_ks_coeff[i_n_ks] * k_lp**n_ks[i_n_ks]), 
                    c=pretty_colors[i_n_ks], 
                    ls='--', 
                    lw=2, 
                    label = r"$\mathtt{k^{"+str(n_ks[i_n_ks])+"}}$ Contrib.")

        sub.plot(k_lp, np.array([delp_poly(k_i_lp) for k_i_lp in k_lp]), 
                c='red', ls='--', lw=4, label="Polynomial (Analytic)")
        
        sub.vlines(0.3, -400., 400., color='gray', linestyle='--', linewidth=3)

    if ell == 0: 
        sub.set_ylim([-1000., 100.])
        sub.set_xlim([10**-3,10**0])
    elif ell == 2: 
        sub.set_ylim([-100.0, 350.])
        sub.set_xlim([10**-3,10**0])
    else: 
        raise NotImplementedError
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(ell)+"}^{corr}(k)}$", fontsize=30)
    
    if ell == 0: 
        sub.legend(loc='lower right', scatterpoints = 1)
    elif ell == 2: 
        sub.legend(loc='upper left', scatterpoints = 1)
     
    fig.savefig(
            ''.join([
                'figure/', 
                'qaplot_delP_'+str(ell)+'corr_lp_components', 
                '_ktrust', str(round(k_trust, 2)), '_inf_nseriesbox_poly', extrap_str, 
                '.png']), 
            bbox_inches="tight")
    #plt.show()
    plt.close()

def plot_delPcorr_ktrust_poly_bestfit(ell, k_trust=0.5):
    '''
    Compare the empirical contribution of delPcorr from q range k_trust to infinity 
    with the polynomial bestfit 
    '''
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(17,8))
    sub = fig.add_subplot(111)

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

    sub.plot(k_val, DelP_tot_k(k_val), 
            c= 'black',
            ls= ':',
            lw=2, 
            label=r"Total Nseries $\mathtt{\Delta P_{"+str(ell)+"}}$")

    sub.plot(k_val, DelP_corr, 
            c= 'black',
            ls= '-.',
            lw=2, 
            label=r"(Nseries $\mathtt{\Delta P_{"+str(ell)+"}) - \Delta P_{"+str(ell)+"}^{uncorr}}$")

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
    
    DelP_corr_ktrust_inf_label = ''.join([ 
        r"$\mathtt{\Delta P_", str(ell), r"^{corr}(k) ", 
        r"\big|_{k_{trust} = ", str(round(k_trust, 2)), "}^{q = q_{max}}}$"])
    
    # plot DelP_corr integrated from q = ktrust - inf (empirical)
    sub.plot(k_val, DelP_corr_ktrust_inf,
            c= 'black',
            ls= '--',
            lw=4, 
            label= DelP_corr_ktrust_inf_label
            )
    
    # plot DelP_corr integrated from q = ktrust - inf (integrated)
    coeff_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/',
        'k^n_coeff',
        '.ell', str(ell), 
        '.ktrust', str(round(k_trust,2)), 
        '.noextrap.dat'])
    n_ks, n_ks_coeff = np.loadtxt(coeff_file, skiprows=1, unpack=True, usecols=[0, 1])
    
    integ_poly_text = r'Integrated : $'
    for i_c, c_integ in enumerate(n_ks_coeff):
        if i_c != 0: 
            if c_integ < 0: 
                integ_poly_text += str(round(c_integ, 2))+'k^{'+str(int(n_ks[i_c]))+'}'
            else:
                integ_poly_text += '+'+str(round(c_integ, 2))+'k^{'+str(int(n_ks[i_c]))+'}'
        else: 
            integ_poly_text += str(round(c_integ, 2))+'k^'+str(int(n_ks[i_c]))

        if i_c == 2: 
            integ_poly_text += '$\n \t$'
    integ_poly_text += '$'

    delp_poly = lambda kk: \
            np.sum(np.vstack(np.array(
                [n_ks_coeff[ii]*kk**n_ks[ii] for ii in range(len(n_ks)) ] 
                )).T, axis=1)

    sub.plot(k_val, delp_poly(k_val), c='red', ls='--', lw=4, label="Polynomial (Integrated)")

    # plot DelP_corr integrated from q = ktrust - inf (bestfit)
    n_ks, coeff_bestfit = delPcorr_ktrust_inf_bestfit(ell, k_trust = k_trust)
    
    bestfit_poly_text = r'Best-fit : $'
    for i_c, c_bestfit in enumerate(coeff_bestfit):
        if i_c != 0: 
            if c_bestfit < 0: 
                bestfit_poly_text += str(round(c_bestfit, 2))+'k^{'+str(int(n_ks[i_c]))+'}'
            else: 
                bestfit_poly_text += '+'+str(round(c_bestfit, 2))+'k^{'+str(int(n_ks[i_c]))+'}'
        else: 
            bestfit_poly_text += str(round(c_bestfit, 2))+'k^'+str(int(n_ks[i_c]))
        
        if i_c == 2: 
            bestfit_poly_text += '$\n \t $'
    bestfit_poly_text += '$'
    
    delp_poly = lambda kk: \
            np.sum(np.vstack(np.array(
                [coeff_bestfit[ii]*kk**n_ks[ii] for ii in range(len(n_ks)) ] 
                )).T, axis=1)
    
    sub.plot(k_val, delp_poly(k_val), c='green', ls='--', lw=4, label="Polynomial (Best-fit)")

    sub.vlines(0.3, -400., 400., color='gray', linestyle='--', linewidth=3)

    sub.text(0.015, 250, bestfit_poly_text, fontsize=20)
    sub.text(0.015, 200, integ_poly_text, fontsize=20)
    if ell == 0: 
        sub.set_ylim([-1000., 100.])
        sub.set_xlim([10**-3,10**0])
    elif ell == 2: 
        sub.set_ylim([-100.0, 350.])
        sub.set_xlim([10**-3,10**0])
    else: 
        raise NotImplementedError
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(ell)+"}^{corr}(k)}$", fontsize=30)
    
    if ell == 0: 
        sub.legend(loc='lower right', scatterpoints = 1)
    elif ell == 2: 
        sub.legend(loc='upper left', scatterpoints = 1)
     
    fig.savefig(
            ''.join([
                'figure/', 
                'qaplot_delP_'+str(ell)+'corr', 
                '.ktrust', str(round(k_trust, 2)), '_inf',  
                '.poly_integrated.poly_bestfit', 
                '.noextrap.png']), 
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

def plot_corrected_Pk(ell, mock='nseriesbox', noextrap='', Ngrid=960):
    '''
    Compare the Top-Hat Convolvution Corrected P_ell(k) with the true P(k) 
    (P_corr/P_true). Compare to delP/P 
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
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
    k_true = true_spec.k
    Pk_true = getattr(true_spec, 'p'+str(ell)+'k')  # P^true
    upw_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'upweight'},
            'spec': specdict
            }
    upw_spec = AvgSpec(20, 'pk', upw_cat_corr)
    upw_spec.read()
    k_upw = upw_spec.k
    Pk_upw = getattr(upw_spec, 'p'+str(ell)+'k')    # P^upw

    fig = plt.figure(figsize=(17,8))
    sub = fig.add_subplot(111)
    
    for mock in ['nseriesbox']:#, 'nseries']:
        if mock == 'nseries': # this needs to be fixed
            lps = [0, 2, 4]
            lstyle = '--'
        elif mock == 'nseriesbox': 
            lps = [0, 2, 4, 6, 8, 10]
            lstyle = '-'

        for i_lp, lp in enumerate(lps):

            if mock == 'nseriesbox': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                    'corrdelP', 
                    str(ell), 
                    'k_lp', 
                    str(lp), 
                    noextrap, 
                    '_AVG_power3600z_BoxN.dat.p'
                    ])
            elif mock == 'nseries': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/', 
                    'corrdelP', 
                    str(ell), 
                    'k_lp', 
                    str(lp), 
                    noextrap, 
                    '_AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600.p'
                    ])

            k_lpcomp, corrdelP_lpcomp = pickle.load(open(pickle_file, 'rb'))
            if mock == 'nseriesbox': 
                print 'lp = ', lp
                print pickle_file 
                print corrdelP_lpcomp.min(), corrdelP_lpcomp.max()
                
            if i_lp == 0: 
                corrdelP = corrdelP_lpcomp
            else: 
                corrdelP += corrdelP_lpcomp

        # del P^uncorr
        if mock == 'nseries': 
            uncorrdelPk_pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/',
                'uncorrdelP', str(ell), 'k_AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600.p'
                ])
        elif mock == 'nseriesbox': 
            uncorrdelPk_pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                'uncorrdelP', str(ell), 'k_AVG_power3600z_BoxN.dat.p'
                ])

        uncorr_k, uncorrdelP = pickle.load(open(uncorrdelPk_pickle_file, 'rb'))

        if mock == 'nseriesbox': 
            uncorr_label = 'Uncorrelated'
            uncorr_lstyle = '-.'
            tot_lstyle = '-'
        else: 
            uncorr_label = None
            tot_label = None
            tot_lstyle = '--'
            
        corr_label = 'Tophat Conv.'
        delpoverp_label = "Sample Variance"
        
        delP_interp = interp1d(k_lpcomp, corrdelP + uncorrdelP, kind='cubic')
        krange = np.where(k_upw > k_lpcomp.min())
        Pk_corr = Pk_true[krange] + delP_interp(k_upw[krange])
        
        #sub.scatter(k_upw, Pk_true/Pk_upw, c= pretty_colors[2], lw=0, label="True")
        sub.scatter(k_upw, Pk_upw/Pk_true, c= pretty_colors[2], lw=0,
                label="Upweighted")

        #sub.scatter(k_upw[krange], Pk_corr/Pk_upw[krange], 
        #        c= pretty_colors[3],
        #        lw=0, 
        #        label=corr_label
        #        )
        sub.scatter(k_upw[krange], (Pk_upw[krange] - (Pk_corr-Pk_true[krange]))/Pk_true[krange], 
                c= pretty_colors[3],
                lw=0, 
                label=corr_label
                )

        Pk_err = true_spec.stddev()

        sub.plot(k_true, 1+Pk_err/Pk_true, 
                c= 'black',
                ls= '--',
                lw=2,
                label=delpoverp_label
                )
        if ell == 0:
            sub.plot(k_true, 1-Pk_err/Pk_true, c= 'black', ls= '--', lw=2, label=None)
        
        del corrdelP

    if ell == 0: 
        sub.set_ylim([0.6, 1.4])
        sub.set_xlim([10**-3,10**0])
    elif ell == 2: 
        sub.set_ylim([0.0, 2.0])
        sub.set_xlim([10**-3,10**0])
    else: 
        raise NotImplementedError
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(
            r"$\mathtt{P_{"+str(ell)+"}^{corr}(k)/P_{"+str(ell)+"}^{true}(k)}$", 
            fontsize=30)
    
    if ell == 0: 
        sub.legend(loc='lower right', scatterpoints = 1)
    elif ell == 2: 
        sub.legend(loc='upper left', scatterpoints = 1)
    elif ell == 4: 
        sub.legend(loc='lower left', scatterpoints = 1)
     
    fig.savefig(
            ''.join([
                'figure/', 
                'qaplot_P_', 
                str(ell), 
                'corr_over_Ptrue', 
                noextrap, 
                '_', 
                mock, 
                '.png']), 
            bbox_inches="tight")
    plt.close()

def plot_corrected_Pk_qmaxtest(ell, q_max, mock='nseriesbox', noextrap='', Ngrid=960):
    '''
    Compare the Top-Hat Convolvution Corrected P_ell(k) *integrated over 
    specified q_max* with the true P(k) (P_corr/P_true). Compare to delP/P 
    to determine an acceptable q_max that corrects within the sample 
    variance. 
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
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
    k_true = true_spec.k
    Pk_true = getattr(true_spec, 'p'+str(ell)+'k')  # P^true
    upw_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'upweight'},
            'spec': specdict
            }
    upw_spec = AvgSpec(20, 'pk', upw_cat_corr)
    upw_spec.read()
    k_upw = upw_spec.k
    Pk_upw = getattr(upw_spec, 'p'+str(ell)+'k')    # P^upw

    fig = plt.figure(figsize=(17,8))
    sub = fig.add_subplot(111)
    
    for mock in ['nseriesbox']:#, 'nseries']:
        if mock == 'nseries': # this needs to be fixed
            lps = [0, 2, 4]
            lstyle = '--'
        elif mock == 'nseriesbox': 
            lps = [0, 2, 4, 6, 8, 10]
            lstyle = '-'

        for i_lp, lp in enumerate(lps):

            if mock == 'nseriesbox': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                    'corrdelP', 
                    str(ell), 
                    'k_lp', 
                    str(lp), 
                    '_qmax', 
                    str(round(q_max, 2)), 
                    '_AVG_power3600z_BoxN.dat.p'
                    ])
            elif mock == 'nseries': 
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/', 
                    'corrdelP', 
                    str(ell), 
                    'k_lp', 
                    str(lp), 
                    '_qmax', 
                    str(round(q_max, 2)), 
                    '_AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600.p'
                    ])

            k_lpcomp, corrdelP_lpcomp = pickle.load(open(pickle_file, 'rb'))
            if mock == 'nseriesbox': 
                print 'lp = ', lp
                print pickle_file 
                print corrdelP_lpcomp.min(), corrdelP_lpcomp.max()
                
            if i_lp == 0: 
                corrdelP = corrdelP_lpcomp
            else: 
                corrdelP += corrdelP_lpcomp

        # del P^uncorr
        if mock == 'nseries': 
            uncorrdelPk_pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/',
                'uncorrdelP', str(ell), 'k_AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600.p'
                ])
        elif mock == 'nseriesbox': 
            uncorrdelPk_pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                'uncorrdelP', str(ell), 'k_AVG_power3600z_BoxN.dat.p'
                ])

        uncorr_k, uncorrdelP = pickle.load(open(uncorrdelPk_pickle_file, 'rb'))

        if mock == 'nseriesbox': 
            uncorr_label = 'Uncorrelated'
            uncorr_lstyle = '-.'
            tot_lstyle = '-'
        else: 
            uncorr_label = None
            tot_label = None
            tot_lstyle = '--'
            
        corr_label = 'Tophat Conv.'
        delpoverp_label = "Sample Variance"
        
        delP_interp = interp1d(k_lpcomp, corrdelP + uncorrdelP, kind='cubic')
        krange = np.where(k_upw > k_lpcomp.min())
        Pk_corr = Pk_true[krange] + delP_interp(k_upw[krange])
        
        #sub.scatter(k_upw, Pk_true/Pk_upw, c= pretty_colors[2], lw=0, label="True")
        sub.scatter(k_upw, Pk_upw/Pk_true, c= pretty_colors[2], lw=0,
                label="Upweighted")

        #sub.scatter(k_upw[krange], Pk_corr/Pk_upw[krange], 
        #        c= pretty_colors[3],
        #        lw=0, 
        #        label=corr_label
        #        )
        sub.scatter(k_upw[krange], (Pk_upw[krange] - (Pk_corr-Pk_true[krange]))/Pk_true[krange], 
                c= pretty_colors[3],
                lw=0, 
                label=corr_label
                )

        Pk_err = true_spec.stddev()

        sub.plot(k_true, 1+Pk_err/Pk_true, 
                c= 'black',
                ls= '--',
                lw=2,
                label=delpoverp_label
                )
        if ell == 0:
            sub.plot(k_true, 1-Pk_err/Pk_true, c= 'black', ls= '--', lw=2, label=None)
        
        del corrdelP

    if ell == 0: 
        sub.set_ylim([0.6, 1.4])
        sub.set_xlim([10**-3,10**0])
    elif ell == 2: 
        sub.set_ylim([0.0, 2.0])
        sub.set_xlim([10**-3,10**0])
    else: 
        raise NotImplementedError
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(
            r"$\mathtt{P_{"+str(ell)+"}^{corr}(k)/P_{"+str(ell)+"}^{true}(k)}$", 
            fontsize=30)
    
    if ell == 0: 
        sub.legend(loc='lower right', scatterpoints = 1)
    elif ell == 2: 
        sub.legend(loc='upper left', scatterpoints = 1)
    elif ell == 4: 
        sub.legend(loc='lower left', scatterpoints = 1)
     
    fig.savefig(
            ''.join([
                'figure/', 
                'qaplot_P_', 
                str(ell), 
                'corr_over_Ptrue', 
                '_qmax', 
                str(round(q_max, 2)),
                '_', 
                mock, 
                '.png']), 
            bbox_inches="tight")
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

def DelP_corr_tot(mock, ell, Ngrid=960):
    '''
    Write out DelP_corr_tot from mock catalogs 
    '''

    if mock == 'nseries': 
        specdict = {
                'P0': 20000,
                'Lbox': 3600, 
                'Ngrid': Ngrid, 
                'ell': ell 
                }

        # true
        true_cat_corr = {
                'catalog': {'name': 'nseries', 'n_mock': 1}, 
                'correction': {'name': 'true'},
                'spec': specdict
                }
        true_spec = AvgSpec(20, 'pk', true_cat_corr)
        true_spec.read()
        k_true = true_spec.k
        Pk_true = getattr(true_spec, 'p'+str(ell)+'k')
        # upweighted
        upw_cat_corr = {
                'catalog': {'name': 'nseries', 'n_mock': 1}, 
                'correction': {'name': 'upweight'},
                'spec': specdict
                }
        upw_spec = AvgSpec(20, 'pk', upw_cat_corr)
        upw_spec.read()
        k_upw = upw_spec.k
        Pk_upw = getattr(upw_spec, 'p'+str(ell)+'k')

        DelP_tot = Pk_upw - Pk_true     # total Delta P from Nseries 
        
    DelP_tot_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'DelP_', str(ell), 'k_tot.', mock, '.p'])

    pickle.dump([k_upw, DelP_tot], open(DelP_tot_file, 'wb'))
    return None

if __name__=="__main__": 
    plot_delPcorr_ktrust_poly_bestfit(2, k_trust=0.5)

    #DelP_corr_tot('nseries', 2, Ngrid=960)

    #for ktrust in [0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 2.0, 3.0, 4.0]: 
    #    plot_delPcorr_ktrust_poly(2, k_trust=ktrust)
    #    plot_delPcorr_ktrust_poly(2, k_trust=ktrust, noextrap=True)

    #for k_trust in [0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 2.0, 3.0, 4.0]: 
    #    plot_delPcorr_ktrust(2, k_trust=k_trust)

    #for qmax in [0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 2.0, 3.0, 4.0]: 
    #    plot_corrected_Pk_qmaxtest(0, qmax)
    #    plot_corrected_Pk_qmaxtest(2, qmax)
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
    #plot_delP_lp_comp_nseries(0, Ngrid=960, noextrap='')
    #plot_delP_lp_comp_nseries(2, Ngrid=960, noextrap='')
    #plot_delP_lp_comp_nseries(0, Ngrid=960, noextrap='_noextrap')
    #plot_delP_lp_comp_nseries(2, Ngrid=960, noextrap='_noextrap')
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
