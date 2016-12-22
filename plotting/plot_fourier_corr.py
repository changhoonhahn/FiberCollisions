'''

Plots for xi tohat fiber collision correction methods 

'''
import time
import pickle
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

# local ----
from corr_spec.corr_average import CorrAvgSpec as AvgSpec
from fourier_corr.test_fourier_corr import delPcorr_bestfit
from fourier_corr.test_fourier_corr import delPuncorr_bestfit
from fourier_corr.test_fourier_corr import delPcorr_ktrust_inf_bestfit

# plotting ----
from defutility.plotting import prettyplot 
from defutility.plotting import prettycolors 

def plot_2pcf_residual(scale='verysmall', fs=0.7, rc=0.43, proj=False):
    ''' 2PCF plotting for the true and upweighted mock catalogs. 
    Plots include the residual plot (1 + xi^fc)/(1 + xi^true). 
    So far, the 2PCF is computed only for the Nseries mock 
    catalogs. 

    Parameters
    ----------
    scale : string
        String that specifies the scale of the 2PCF. This string
        dictates that contour range used for the plotting. 
    '''
    # scale of the 2PCF 
    if scale == 'large': 
        contour_range = np.arange(-0.05, 0.05, 0.005)
        n_rp, n_pi = 40, 40 
        n_mocks = range(1, 21)
    elif scale == 'small': 
        contour_range = np.arange(-0.1, 0.11, 0.01)
    elif scale == 'smaller': 
        contour_range = np.arange(-0.5, 0.5, 0.05)
    elif scale == 'verysmall': 
        contour_range = 20
        n_rp, n_pi = 20, 20 
        n_mocks = range(1, 21)
    elif scale == '5x5': 
        contour_range = 20
        n_rp, n_pi = 100, 100  
        n_mocks = range(1, 6)
    else: 
        raise NotImplementedError

    for corr in ['true', 'upweighted']:    # for each correction
        for i_mock in n_mocks:  # for each mocks
            corr_file = ''.join([
                '/mount/riachuelo1/hahn/2pcf/corr/',
                'corr_2pcf_CutskyN', 
                str(i_mock), '.', 
                corr, 
                '.cute2pcf.', 
                scale, 'scale.dat'
                ])
            print corr_file 

            tpcf = np.loadtxt(corr_file, unpack=True)
    
            if i_mock == n_mocks[0]: 
                rp_bins = tpcf[0].reshape(n_rp, n_pi)[0]
                pi_bins = tpcf[1].reshape(n_rp, n_pi)[:,0]
                r_p, pi = np.meshgrid(rp_bins, pi_bins)

                twop_corr = tpcf[2].reshape(n_rp, n_pi)

            else: 
                twop_corr += tpcf[2].reshape(n_rp, n_pi)

        twop_corr /= np.float(len(n_mocks))  # average 2PCF

        if corr == 'true': 
            true_corr = twop_corr.T

            continue

        prettyplot()
        pretty_colors = prettycolors()

        fig = plt.figure(figsize=(15,10))
        sub = fig.add_subplot(111)

        # contour of 1 - (1 + xi^fc)/(1+xi^true)
        residual_2pcf = 1.0 - (1.0 + twop_corr.T)/(1.0 + true_corr)
        
        if not proj:
            cont = sub.contourf(
                    r_p, pi, 
                    residual_2pcf, 
                    contour_range, 
                    cmap=plt.cm.afmhot
                    )
            plt.colorbar(cont)

            sub.vlines(0.4, 0.0, np.max(r_p), lw=4, linestyle='--', color='red')

            sub.set_xlabel('$\mathtt{r_{p}}$', fontsize=40)
            sub.set_ylabel('$\pi$', fontsize=40)
            sub.set_xlim([np.min(rp_bins), np.max(rp_bins)])
            sub.set_ylim([np.min(pi_bins), np.max(pi_bins)])
            
            sub.set_title(
                    r"$1 - (1 + \xi^\mathtt{"+corr.upper()+r"})/(1+ \xi^\mathtt{TRUE})$", 
                    fontsize=40
                    )
            proj_str = ''
        else: 
            proj_resid = np.sum(residual_2pcf, axis=0)
            sub.plot(r_p[0], proj_resid/np.float(n_pi), 
                    lw=4, c=pretty_colors[2], label=r'Projected $\xi(r_p, \pi)$') 
            
            # top hat function for comparison  
            rbin = np.linspace(0.0, 50.0, 10000)
            inhat = np.where(rbin < rc)
            tophat = np.repeat(0., len(rbin))
            tophat[inhat] = fs
            sub.plot(rbin, tophat, lw=3, ls='--', c='k', 
                    label=r'$f_s='+str(round(fs, 2))+'$ Top Hat')
            
            # x-axis
            sub.set_xlim([0, np.max(r_p[0])])
            sub.set_xlabel(r'$r_p$', fontsize=30)
            # y-axis
            sub.set_ylim([-0.1, 1.0])
            sub.set_ylabel(r"$1 - (1 + \xi^\mathtt{"+corr.upper()+r"})/(1+ \xi^\mathtt{TRUE})$", fontsize=20)

            sub.legend(loc='upper right')
            
            proj_str = '.rp_proj'
    
        fig_name = ''.join([
            'figure/',
            '2pcf_Nseries_', corr, '_', str(len(n_mocks)), 'mocks.', scale, 
            proj_str, 
            '.tophat_comparison.png'
            ])

        fig.savefig(fig_name, bbox_inches="tight")
        plt.close()

def plot_delPuncorr(ell, fs=1.0, rc=0.43): 
    ''' Plot Del P^{uncorr} with the empirical Del P^{tot} from the Nseries mocks
    '''
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    # average P_l(k) and P_l^upw(k) from Nseries mocks
    k_true, Pk_true, k_upw, Pk_upw, Pk_err = nseries_Pk(ell, n_mock=n_moc, Ngrid=nseries_Ngrid)
    delP_tot = Pk_upw - Pk_true
    sub.plot(k_true, delP_tot, c= 'gray', lw=2, label='Total')

    # Del P^uncorr for given fs and rc (hardcoded for Nseries box) 
    pickle_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/',
        'uncorrdelP', str(ell), 'k', 
        '.fs', str(round(fs,2)), 
        '.rc', str(round(rc,2)), 
        '.p'
        ])
    kt, uncorrdelPk = pickle.load(open(pickle_file, 'rb'))

    sub.plot(kt, uncorrdelPk, c='k', lw=2, label=r'$\Delta P^{uncorr}$')
    
    # x-axis 
    sub.set_xlim([10**-3,10**1])
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    # y-axis
    sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(ell)+"}(k)}$", fontsize=30)
    if ell == 0: 
        sub.set_ylim([-1000., 0])
    elif ell == 2: 
        sub.set_ylim([0., 1500.])
    # legend 
    sub.legend(loc='lower right', scatterpoints = 1)
     
    fig.savefig(
            ''.join(['figure/', 
                'qaplot_delP', str(ell), 'uncorr', 
                '.fs', str(round(fs,2)), 
                '.rc', str(round(rc,2)), '.png']), 
           bbox_inches="tight")
    plt.close()

def plot_delPcorr(ell, fs=1.0, rc=0.43, Ngrid=3600, pkmu=True, n_mock=20, nseries_Ngrid=960,
        fold=10, rebin=20):
    ''' Del P_ell^corr
    '''
    # figure
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(17,8))
    sub = fig.add_subplot(111)

    # Empirical Del P^corr
    k_true, Pk_true, k_upw, Pk_upw, Pk_err = nseries_Pk(ell, n_mock=n_mock, Ngrid=nseries_Ngrid)
    DelP_tot = Pk_upw - Pk_true
    DelP_tot_k = interp1d(k_upw, DelP_tot) # Del P (k)
    # subtract out DelP_uncorr 
    k_uncorr, delP_uncorr = DelP_uncorr(ell, fs=fs, rc=rc, mock='nseriesbox')
    krange = np.where((k_uncorr > k_upw[0]) & (k_uncorr < k_upw[-1]))
    k_val = k_uncorr[krange]

    DelP_corr = DelP_tot_k(k_val) - delP_uncorr[krange]      # Del P^corr (k)
    
    sub.scatter(k_val, DelP_corr, c= 'black', label=r'Empirical from Nseries')
    
    pkmu_str = ''
    if pkmu:    # Del P^corr from P(k,mu) 
        if fold: 
            fold_str = ''.join(['.', str(fold), 'fold.', str(rebin), 'rebin'])
        else: 
            fold_str = ''
        pickle_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/Box/', 
            'delPcorr_pkmu', 
            '.l', str(ell), 
            '.fs', str(round(fs,2)), 
            '.rc', str(round(rc,2)), 
            '.Ngrid', str(Ngrid), 
            fold_str, 
            '.p'])
        delPcorrk = interp1d(k_val, DelP_corr) # Del P (k)
        k_pkmu, delPcorr_pkmu = pickle.load(open(pickle_file, 'rb'))
        sub.plot(k_pkmu, delPcorr_pkmu, c=pretty_colors[2], lw=4, 
                label=r'From $\mathtt{P(k, mu)}$')
        pkmu_str = '.pkmu'

        withinrange = np.where((k_pkmu < k_val[-1]) & (k_pkmu > k_val[0]))
        for i_kk, kk in enumerate(k_pkmu[withinrange]):
            if kk > 0.1: 
                if ell == 0: 
                    print kk, 100*(delPcorrk(kk) - delPcorr_pkmu[i_kk])/delPcorr_pkmu[i_kk]
                else: 
                    print kk, (delPcorrk(kk) - delPcorr_pkmu[i_kk])
    
    # x-axis
    sub.set_xscale("log") 
    sub.set_xlim([10**-3,10**0])
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    # y-axis
    if ell == 0: 
        sub.set_ylim([-600., 200.])
        sub.legend(loc='lower left', scatterpoints = 1)
    elif ell == 2: 
        sub.set_ylim([-100.0, 350.])
        sub.legend(loc='upper left', scatterpoints = 1)
    else: 
        raise NotImplementedError
    sub.set_ylabel(r"$\mathtt{\Delta P_"+str(ell)+"^{corr}(k)}$", fontsize=30)

    fig_file = ''.join(['figure/', 
        'qaplot_delP_'+str(ell)+'corr', 
        '.empirical_Ngrid', str(nseries_Ngrid), 
        '.fs', str(round(fs,2)), 
        '.rc', str(round(rc,2)),
        pkmu_str, 
        fold_str, 
        '.png'
        ])
    fig.savefig(fig_file, bbox_inches="tight", dpi=200)
    plt.close()
    return None

def plot_delP_boxcorr(ell, fs=1.0, rc=0.43):
    ''' Del P_ell^corr
    '''
    # figure
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(17,8))
    sub = fig.add_subplot(111)

    # Empirical Del P^corr
    k_true, Pk_true, k_upw, Pk_upw, Pk_err = nseries_Pk(ell)
    DelP_tot = Pk_upw - Pk_true
    DelP_tot_k = interp1d(k_upw, DelP_tot) # Del P (k)
    # subtract out DelP_uncorr 
    k_uncorr, delP_uncorr = DelP_uncorr(ell, fs=fs, rc=rc, mock='nseriesbox')
    krange = np.where((k_uncorr > k_upw[0]) & (k_uncorr < k_upw[-1]))
    k_val = k_uncorr[krange]

    DelP_corr = DelP_tot_k(k_val) - delP_uncorr[krange]      # Del P^corr (k)
    
    sub.scatter(k_val, DelP_corr, c= 'black', label=r'Empirical from Nseries')
   
    d_DelP_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'd_DelP', str(ell), 'k_boxcorr', 
        '.fs', str(round(fs, 2)), 
        '.rc', str(round(rc, 2)), 
        '.p'])
    kk, d_delP = pickle.load(open(d_DelP_file, 'rb'))
    sub.plot(kk, d_delP, c=pretty_colors[2], label=r'$\mathtt{d(\Delta P_l(k))}$')
    
    # x-axis
    sub.set_xscale("log") 
    sub.set_xlim([10**-3,10**0])
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    # y-axis
    if ell == 0: 
        sub.set_ylim([-600., 200.])
        sub.legend(loc='lower left', scatterpoints = 1)
    elif ell == 2: 
        sub.set_ylim([-100.0, 350.])
        sub.legend(loc='upper left', scatterpoints = 1)
    else: 
        raise NotImplementedError
    sub.set_ylabel(r"$\mathtt{\Delta P_"+str(ell)+"^{corr}(k)}$", fontsize=30)

    fig_file = ''.join(['figure/', 
        'qaplot_d_delP'+str(ell)+'_boxcorr', 
        '.fs', str(round(fs, 2)), 
        '.rc', str(round(rc, 2)), 
        '.png'
        ])
    fig.savefig(fig_file, bbox_inches="tight", dpi=200)
    plt.close()
    return None

def plot_delPcorr_untrusted(ell, k_trust=0.5, poly_exp=True, integrated=False, 
        bestfit_poly=False, bestfit_poly_param=None):
    ''' Compare the untrusted Del P^corr integrated over the q range 
    k_trust to infinity calculated from, 

    the polynomial expansion with integrated coefficients 
    versus 
    the actually integrated values
    '''
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(17,8))
    sub = fig.add_subplot(111)
    
    k_true, Pk_true, k_upw, Pk_upw, Pk_err = nseries_Pk(ell)
    # Del P^total 
    DelP_tot = Pk_upw - Pk_true
    DelP_tot_k = interp1d(k_upw, DelP_tot)
    Pk_err_k = interp1d(k_upw, Pk_err)
    # subtract out DelP_uncorr 
    k_uncorr, DelP_uncorr = DelP_uncorr(ell, mock='nseriesbox')
    # subtract out delP_corr |_q = 0^q = k_trust. 
    k_trust, DelP_corr_trust = DelP_corr_0_ktrust(ell, k_trust, mock='nseriesbox')

    krange = np.where((k_uncorr > k_upw[0]) & (k_uncorr < k_upw[-1]))
    k_val = k_uncorr[krange]

    # Empirical Del P^corr |ktrust to inf
    DelP_corr_ktrust_inf = DelP_tot_k(k_val) - DelP_uncorr[krange] - DelP_corr_trust[krange]
    sub.scatter(k_val, DelP_corr_ktrust_inf, c= 'black', label='Empirical')
    
    if poly_exp:    # Polynomial Expansion DelP_corr |q = ktrust - inf
        coeff_file = ''.join(['/mount/riachuelo1/hahn/power/Nseries/Box/',
            'k^n_coeff', 
            '.ell', str(ell), 
            '.ktrust', str(round(k_trust, 2)), 
            '.noextrap.dat'])
        n_ks, n_ks_coeff = np.loadtxt(coeff_file, skiprows=1, unpack=True, usecols=[0, 1])
        # Integrated Polynomial Expression
        integ_poly_text = r'Integrated : $'
        for i_c, c_integ in enumerate(n_ks_coeff):
            if i_c != 0: 
                if c_integ < 0: 
                    integ_poly_text += str(round(c_integ, 2))+'k^{'+str(int(n_ks[i_c]))+'}'
                else:
                    integ_poly_text += '+'+str(round(c_integ, 2))+'k^{'+str(int(n_ks[i_c]))+'}'
            else: 
                if n_ks[i_c] == 0: 
                    integ_poly_text += str(round(c_integ, 2))
                else:
                    integ_poly_text += str(round(c_integ, 2))+'k^'+str(int(n_ks[i_c]))
            if i_c == 2: 
                integ_poly_text += '$\n \t$'
        integ_poly_text += '$'

        delp_poly = lambda kk: \
                np.sum(np.vstack(np.array(
                    [n_ks_coeff[ii]*kk**n_ks[ii] for ii in range(len(n_ks)) ] 
                    )).T, axis=1)
        for i_n_ks in xrange(len(n_ks)):    # polynomial k^n components 
            sub.plot(k_val, 
                    np.array(n_ks_coeff[i_n_ks] * k_val**n_ks[i_n_ks]), 
                    c=pretty_colors[i_n_ks], 
                    lw=2, 
                    label = r"$\mathtt{k^{"+str(int(n_ks[i_n_ks]))+"}}$")
        # total polynomial
        upto_ktrust = np.where(k_val < ktrust)
        sub.plot(k_val[upto_ktrust], delp_poly(k_val[upto_ktrust]), c='red', ls='-', lw=4, label="Polynomial Expansion")
    
    if integrated:      # Integrated Del P^corr|ktrust to qmax 
        pickle_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/Box/', 
            'corrdelP', str(ell), 'k', 
            '.ktrust', str(round(k_trust, 2)), '_qmax.', 
            'AVG_power3600z_BoxN.dat',
            '.p'])
        print pickle_file
        k_integ, delPcorr = pickle.load(open(pickle_file, 'rb'))

        sub.plot(k_integ, delPcorr, c='blue', ls='-', lw=4, label="Integrated")

    if bestfit_poly:    # Bestfit Polynomial t Del P^corr|ktrust to qmax
        n_ks, coeff_bestfit = delPcorr_bestfit(
                ell, 
                highest_power = bestfit_poly_param['highest_power'], 
                k_fit_max = bestfit_poly_param['k_fit_max'])
    
        bestfit_poly_text = r'Best-fit : $'
        for i_c, c_bestfit in enumerate(coeff_bestfit):
            if i_c != 0: 
                if c_bestfit < 0: 
                    bestfit_poly_text += str(round(c_bestfit, 2))+'k^{'+str(int(n_ks[i_c]))+'}'
                else: 
                    bestfit_poly_text += '+'+str(round(c_bestfit, 2))+'k^{'+str(int(n_ks[i_c]))+'}'
            else: 
                if n_ks[i_c] == 0: 
                    bestfit_poly_text += str(round(c_bestfit, 2))
                else:
                    bestfit_poly_text += str(round(c_bestfit, 2))+'k^'+str(int(n_ks[i_c]))
            if len(coeff_bestfit) > 3: 
                if i_c == 2: 
                    bestfit_poly_text += '$\n \t $'
        bestfit_poly_text += '$'

        delp_poly = lambda kk: \
                np.sum(np.vstack(np.array(
                    [coeff_bestfit[ii]*kk**n_ks[ii] for ii in range(len(n_ks)) ] 
                    )).T, axis=1)
        trust_k = k_val[np.where(k_val <= k_trust)]
        sub.plot(trust_k, delp_poly(trust_k), c='green', ls='--', lw=4, label="Polynomial (Best-fit)")

    if ell == 0: 
        sub.vlines(0.3, -600., 400., color='gray', linestyle='--', linewidth=3)
        sub.text(0.01, -500, integ_poly_text, fontsize=20)
        sub.set_ylim([-600., 200.])
        sub.set_xlim([10**-3,10**0])
        sub.legend(loc='lower left', scatterpoints = 1)
    elif ell == 2: 
        sub.vlines(0.3, -400., 400., color='gray', linestyle='--', linewidth=3)
        sub.text(0.005, 250, integ_poly_text, fontsize=20)
        sub.set_ylim([-100.0, 350.])
        sub.set_xlim([10**-3,10**0])
        sub.legend(loc='upper left', scatterpoints = 1)
    else: 
        raise NotImplementedError
    sub.set_xscale("log") 
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(''.join([ 
        r"$\mathtt{\Delta P_", str(ell), r"^{corr}(k) ", 
        r"\big|_{k_{trust} = ", str(round(k_trust, 2)), "}^{q = q_{max}}}$"]), 
        fontsize=30)
    
    fig.savefig(''.join(['figure/', 
        'qaplot_delP_'+str(ell)+'corr', 
        '.ktrust', str(round(k_trust, 2)), '_qmax',  
        '.empirical.poly_expansion.integrated.png']), 
        bbox_inches="tight")
    #plt.show()
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

def DelP_corr_0_ktrust(ell, k_trust, mock='nseriesbox'):
    ''' Del P^{corr} integrated from 0 to k_trust
    '''
    if mock == 'nseriesbox': 
        pickle_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/Box/', 
            'corrdelP', str(ell), 'k', 
            '_qmax', 
            str(round(k_trust,2)),
            '_AVG_power3600z_BoxN.dat.p'
            ])
    else: 
        raise ValueError
    k, DelP_corr = pickle.load(open(pickle_file, 'rb'))
    return [k, DelP_corr]

def DelP_uncorr(ell, fs=0.7, rc=0.43, mock='nseriesbox'): 
    ''' Del P^{uncorr}
    '''
    if mock == 'nseriesbox': 
        DelP_uncorr_file = ''.join([ 
            '/mount/riachuelo1/hahn/power/Nseries/Box/', 
            'uncorrdelP', str(ell), 'k', 
            '.fs', str(round(fs, 2)), 
            '.rc', str(round(rc, 2)),
            '.p'
            ])
    else: 
        raise ValueError
    
    k_uncorr, DelP_uncorr = pickle.load(open(DelP_uncorr_file, 'rb'))
    return [k_uncorr, DelP_uncorr]

def plot_Pk_comp(ell, residual=False, ratio=False, frac_residual=False): 
    ''' Compare P_ell(k) of Nseries versus Nseries Box
    '''
    # Nseries
    output = nseries_Pk(ell, Ngrid=960, n_mock=20)
    k, Pk = output[0], output[1]
    # Nseries Box
    k_box, Pk_box= nseriesbox_Pk(ell)

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    
    if residual: 
        Pk_box_interp = interp1d(k_box, Pk_box) 

        krange = np.where(k > k_box[0])
        Pk_resid = Pk[krange] - Pk_box_interp(k[krange])
        sub.plot(k[krange], Pk_resid, color=pretty_colors[0])
    elif frac_residual: 
        Pk_box_interp = interp1d(k_box, Pk_box) 

        krange = np.where(k > k_box[0])
        Pk_resid = (Pk[krange] - Pk_box_interp(k[krange]))/np.abs(Pk_box_interp(k[krange]))
        sub.plot(k[krange], Pk_resid, color=pretty_colors[0])
    elif ratio:  
        Pk_box_interp = interp1d(k_box, Pk_box) 

        krange = np.where(k > k_box[0])
        Pk_ratio = Pk[krange] / Pk_box_interp(k[krange])
        sub.plot(k[krange], Pk_ratio, color=pretty_colors[0], lw=3)
    else:
        sub.plot(k_box, np.abs(Pk_box), color=pretty_colors[0], label='Nseries Box')
        sub.scatter(k, np.abs(Pk), color=pretty_colors[3], label='Nseries')
    
    # x-axis
    sub.set_xscale("log") 
    sub.set_xlim([10**-3,10**1])
    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    # y-axis
    if residual: 
        sub.set_xlim([10**-3,10**0])
        sub.set_ylim([-200., 200.])
        sub.set_ylabel(
                r"$\mathtt{P^{Nseries}_{"+str(ell)+"} - P^{NseriesBox}_{"+str(ell)+"}}$", 
                fontsize=30)
        plot_str = '.residual'
    elif frac_residual: 
        sub.set_xlim([10**-3,10**0])
        sub.set_ylim([-0.15, 0.15])
        sub.set_ylabel(
                r"$(\mathtt{P^{Nseries}_{"+str(ell)+"} - P^{NseriesBox}_{"+str(ell)+"}})/|\mathtt{P^{NseriesBox}_{"+str(ell)+"}}|$", 
                fontsize=30)
        plot_str = '.frac_residual'
    elif ratio:  
        sub.set_xlim([10**-3,10**0])
        sub.set_ylim([0.8, 1.2])
        sub.set_ylabel(
                r"$\mathtt{P^{Nseries}_{"+str(ell)+"} / P^{NseriesBox}_{"+str(ell)+"}}$", 
                fontsize=30)
        plot_str = '.ratio'
    else: 
        sub.set_yscale("log")
        sub.set_ylim([10**2, 10**5])
        sub.set_ylabel(r"$\mathtt{|P_{"+str(ell)+"}^{True}(k)|}$", fontsize=30)
    
        sub.legend(loc='upper right', scatterpoints = 1)
        plot_str = ''

    fig_file = ''.join(['figure/'
        'P', str(ell), 'k_comparison.Nseries.NseriesBox', 
        plot_str, 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=200)
    plt.close()
    return None

def nseriesbox_Pk(ell): 
    ''' True P_l(k) from Nseries Box
    '''
    # loop through mock realizations
    n_mock = 7
    for i_mock in xrange(1, n_mock+1):
        true_pk_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/Box/', 
            'power3600z_BoxN', str(i_mock), '.dat'
            ])

        l_index = ell/2 + 1
        if ell == 0: 
            l_index *= -1
    
        k, p_ell_k_i = np.loadtxt(true_pk_file, unpack=True, usecols=[0,l_index])

        if i_mock == 1: 
            p_ell_k = p_ell_k_i
        else: 
            p_ell_k += p_ell_k_i

    p_ell_k *= (2.*np.pi)**3/np.float(n_mock)

    return [k, p_ell_k]

def nseries_Pk(ell, Ngrid=960, n_mock=20):
    ''' Import the true and upweighted P(k) for Nseries mock catalog. 
    This is used in many of the modules above for Del P comparisons
    '''
    # True
    specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell}
    true_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'true'},
            'spec': specdict
            }
    true_spec = AvgSpec(n_mock, 'pk', true_cat_corr)
    true_spec.read()

    k_true = true_spec.k
    Pk_true = getattr(true_spec, 'p'+str(ell)+'k')
    Pk_err = true_spec.stddev() # sample variance
    # upweight
    upw_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'upweight'},
            'spec': specdict
            }
    upw_spec = AvgSpec(n_mock, 'pk', upw_cat_corr)
    upw_spec.read()
    k_upw = upw_spec.k
    Pk_upw = getattr(upw_spec, 'p'+str(ell)+'k')
    
    return [k_true, Pk_true, k_upw, Pk_upw, Pk_err]


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
    def plot_delPcorr_ktrust_poly(ell, k_trust=0.5, Ngrid=960):
        '''
        Investigate the contribution of delPcorr from q range k_trust to infinity. 
        delPcorr is being divided into two parts:
            - 0 to k_trust, which can be reliably calculated 
            - k_trust to infinity, which cannot be reliably calculated because 
            we do not trust the models or the extrapolations beyound this point.
            However they it may be possible to have a polynomial fit. 
        '''
        prettyplot()
        pretty_colors = prettycolors()

        fig = plt.figure(figsize=(17,8))
        sub = fig.add_subplot(111)


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
        Pk_true = getattr(true_spec, 'p'+str(ell)+'k')

        upw_cat_corr = {
                'catalog': {'name': 'nseries', 'n_mock': 1}, 
                'correction': {'name': 'upweight'},
                'spec': specdict
                }
        upw_spec = AvgSpec(20, 'pk', upw_cat_corr)
        upw_spec.read()
        k_upw = upw_spec.k
        Pk_upw = getattr(upw_spec, 'p'+str(ell)+'k')
        
        lps = [0, 2, 4, 6, 8, 10]
        for i_lp, lp in enumerate(lps):
            # total delP^corr integrated from 0 to infinity 
            pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                'corrdelP', str(ell), 'k_lp', str(lp), 
                '_noextrap_AVG_power3600z_BoxN.dat.p'
                ])
            print pickle_file
            k_lp, corrdelP_lp = pickle.load(open(pickle_file, 'rb'))
            
            # delP^corr 0 to k_trust
            pickle_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                'corrdelP', str(ell), 'k_lp', str(lp), 
                '_qmax', str(round(k_trust,2)),
                '_AVG_power3600z_BoxN.dat.p'])
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
                ls= ':',
                lw=2, 
                label="Total, $0$ to $\infty$"
                )

        sub.plot(k_lp, corrdelP_ktrust_inf, 
                c= 'black',
                ls= '-',
                lw=4, 
                label="Total, $k_{trust} ="+str(round(k_trust,2))+"$"+" to $\infty$"
                )

        coeff_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/Box/',
            'k^n_coeff',
            '.ell', str(ell), 
            '.ktrust', str(round(k_trust,2)), 
            '.noextrap.dat'])
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
            sub.set_ylim([-1000., 200.])
            sub.set_xlim([10**-3,10**0])
            sub.legend(loc='lower left', scatterpoints = 1)
        elif ell == 2: 
            sub.set_ylim([-100.0, 350.])
            sub.set_xlim([10**-3,10**0])
            sub.legend(loc='upper left', scatterpoints = 1)
        else: 
            raise NotImplementedError
        sub.set_xscale("log") 
        sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
        sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(ell)+"}^{corr}(k)}$", fontsize=30)
        
        fig.savefig(''.join([
            'figure/', 
            'qaplot_delP_'+str(ell)+'corr_lp_components', 
            '_ktrust', str(round(k_trust, 2)), '_inf_nseriesbox_poly.noextrap', 
            '.png']), 
            bbox_inches="tight")
        plt.close()
    def delP_uncorr(ell, k_fit_max=0.2): 
        '''
        '''
        # true Pg(k)
        true_cat_corr = {
                'catalog': {'name': 'nseries', 'n_mock': 1}, 
                'correction': {'name': 'true'},
                'spec': {'P0': 20000, 'Lbox': 3600, 'Ngrid': 960, 'ell': ell}
                }
        true_spec = AvgSpec(20, 'pk', true_cat_corr)
        true_spec.read()
        k_true = true_spec.k
        Pk_true = getattr(true_spec, 'p'+str(ell)+'k')
        Pk_true_k = interp1d(k_true, Pk_true)

        # Del P_uncorr
        DelP_uncorr_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/Box/', 
            'uncorrdelP', str(ell), 'k_AVG_power3600z_BoxN.dat.p'])
        k_uncorr, DelP_uncorr = pickle.load(open(DelP_uncorr_file, 'rb'))
        krange = np.where((k_uncorr > k_true[0]) & (k_uncorr < k_fit_max))
        k_val = k_uncorr[krange]
        
        prettyplot()
        pretty_colors = prettycolors()
        fig = plt.figure(1)
        sub = fig.add_subplot(111)

        sub.plot(k_uncorr, DelP_uncorr, 
                color=pretty_colors[2], lw=4, label=r"$\mathtt{\Delta P^{uncorr}(k)}$")

        bestfit_param = delPuncorr_bestfit(k_val, DelP_uncorr[krange]/Pk_true_k(k_val))
        sub.plot(k_val, Pk_true_k(k_val) * bestfit_param[0]/k_val**2, 
                c='k', ls='--', lw=3, 
                label=r'Best-fit $\mathtt{\frac{A}{k^2} P^{True}_{'+str(ell)+'}(k)}$')
        sub.text(2.*10**-3, 100, r"A = %.2e" %bestfit_param[0], fontsize=20)
        # x-axis
        sub.set_xscale('log')
        sub.set_xlim([10**-3,10**0])
        sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
        # y-axis 
        #sub.set_yscale('log')
        sub.set_ylim([-1000., 200.])
        sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(ell)+"}^{corr}(k)}$", fontsize=30)
        
        sub.legend(loc='lower right')
        
        fig.savefig(
                ''.join(['figure/', 'delP', str(ell), '_uncorr.nongauss_fit',
                    '.kfitmax', str(round(k_fit_max, 2)), '.png']),
                bbox_inches='tight')
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
    def plot_delP_lp_component(ell, mock='nseries', n_mocks=84, Ngrid=960):
        '''
        Compare each of the l' components of Del P_l^corr. 
        l' components are read from pickle files.
        '''
        # average P_l(k) and P_l^upw(k)
        specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell}
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
                raise NotImplementedError
            elif mock == 'nseriesbox': 
                lps = [0, 2, 4, 6]

            for i_lp, lp in enumerate(lps):
                raise ValueError("Code below is outdated; do not run")
                pickle_file = ''.join([
                    '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                    'corrdelP2k_lp', str(lp), '_power3600z_BoxN', str(i_mock), 
                    '.fourier_tophat.fs1.0.rc0.43.kfit4.3.kfixed4.34.dat.p'
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

        sub.plot(k_lpcomp, corrdelP, c= 'gray', lw=2, label='Total')

        sub.plot(k, Pk_upw - Pk, c= 'k', lw=2, label='data')

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
"""

if __name__=="__main__": 
    #plot_delP_boxcorr(0, fs=0.7, rc=0.43)
    #plot_delP_boxcorr(2, fs=0.7, rc=0.43)
    #plot_2pcf_residual(scale='5x5')
    #plot_delPcorr(0, fs=0.7, rc=0.43, Ngrid=3600, pkmu=True, nseries_Ngrid=480, n_mock=10)
    for ell in [0, 2]:
        plot_delPcorr(ell, fs=0.7, rc=0.43, Ngrid=3600, pkmu=True, nseries_Ngrid=960, n_mock=20,
                fold=False, rebin=False)
        plot_delPcorr(ell, fs=0.6, rc=0.43, Ngrid=3600, pkmu=True, nseries_Ngrid=960, n_mock=20,
                fold=10, rebin=20)
    #DelP_corr_tot('nseries', 0, Ngrid=960)
    
    #for ell in [0]:
    #    for ktrust in [0.3, 0.5]: 
    #        plot_delPcorr_ktrust_poly_integ(ell, k_trust=ktrust)
    #for kt in [0.2, 0.3, 0.4, 0.5]:
    #    for n_high in [2, 4]: 
    #        plot_delPcorr_ktrust_poly_bestfit(0, k_trust=kt, highest_power=n_high, k_fit_max=kt)
    #        plot_delPcorr_ktrust_poly_bestfit(2, k_trust=kt, highest_power=n_high, k_fit_max=kt)
