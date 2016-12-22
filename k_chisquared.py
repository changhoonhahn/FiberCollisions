import numpy as np 
import scipy as sp 
from scipy.interpolate import interp1d

# --- Local ---
from fourier_corr import fourier_corr as four_corr 
from corr_spec.corr_average import CorrAvgSpec as AvgSpec


def k_chi2_CM(ell, ktrust=0.3, max_order=10, cm_poly_order='upto2', chi_goal=1., n_mocks=84, Ngrid=960, rebin=10, use='nseries'): 
    ''' Calculate the k value (k_chisquare) where the delta chi^2 value is 1 for the 
    Convolution Method (CM). Note that this calculation only applies for the Nseries mocks!
    '''
    # convolution method parameters
    fs = 0.6 
    rc = 0.43
    fold = 10
    cm_rebin = 20
    
    k_chi, k_w_chi = [], [] 
    chi_w, chi = [], [] 
    
    n_mock_cat = n_mocks
    catdict = {'name': 'nseries', 'n_mock': 1}
    specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell}
    catcorr = {'catalog': catdict, 'spec': specdict}
        
    tru_cat_corr = catcorr.copy()  # P^true(k)
    tru_cat_corr['correction'] = {'name': 'true'}
    tru_avg_spec = AvgSpec(n_mock_cat, 'pk', tru_cat_corr)
    tru_avg_spec.read(rebin=rebin)

    NN_cat_corr = catcorr.copy()    # P^NN(k)
    NN_cat_corr['correction'] = {'name': 'upweight'}
    NN_avg_spec = AvgSpec(n_mock_cat, 'pk', NN_cat_corr)
    NN_avg_spec.read(rebin=rebin)

    k_arr = tru_avg_spec.k
    tru_avg_pk = getattr(tru_avg_spec, ''.join(['p', str(ell), 'k']))
    NN_avg_pk = getattr(NN_avg_spec, ''.join(['p', str(ell), 'k']))
    counts = NN_avg_spec.count

    # Del P Uncorrelated
    k_uncorr, delP_uncorr = four_corr.DelPuncorr(ell, fs=fs, rc=rc, k_arr=k_arr)

    # Del P Corr Trusted Range
    k_corr_kt, delP_corr_kt = four_corr.DelPcorr_pkmu(
            ell, fs=fs, rc=rc, fold=fold, rebin=cm_rebin, 
            dqdmu=True, dmudq=False, ktrust=ktrust)
    corr_kt_interp = interp1d(k_corr_kt, delP_corr_kt, kind='cubic')
    k_range = np.where((k_arr > k_corr_kt[0]) & (k_arr < k_corr_kt[-1]))

    # Del P Corr Untrusted Range Polynomial estimate 
    delpcorr_untrust_orderlp = four_corr.DelPcorr_untrust_poly(
            k_arr, ell=ell, ktrust=ktrust, order=cm_poly_order, 
            fs=fs, rc=rc, fold=10, rebin=50, max_order=max_order)

    # P(k)^CM = Pl^true + Pl^uncorr + Pl^corr_trusted + Pl^corr_untrusted
    cm_pk = tru_avg_pk[k_range] + delP_uncorr[k_range] + corr_kt_interp(k_arr[k_range]) + delpcorr_untrust_orderlp[k_range]
    
    # covariance matrix
    full_cov = NN_avg_spec.covariance(krange=[k_corr_kt[0], k_corr_kt[-1]], rebin=rebin)
    
    # keep covariance matrix within the k range
    k_arr = k_arr[k_range]
    NN_avg_pk = NN_avg_pk[k_range]
    counts = counts[k_range]
    if ell == 0: 
        mono_offset = 0.# 85. #np.mean((cm_pk - NN_avg_pk)[np.where(k_arr < 0.1)])
        print 'monopole offset = ', mono_offset

    # calculate chi^2(k) 
    i_start = (np.abs(k_arr - 0.05)).argmin()
    once = None
    for ik in range(i_start, len(k_arr)): 
        cov = full_cov[:ik, :ik]
        P_resid = cm_pk[:ik] - NN_avg_pk[:ik]
        if ell == 0: 
            P_resid -= mono_offset

        if np.linalg.cond(cov) < 10**20:    # otherwise the cov matrix can't be inverted
            inv_chi2 = np.dot(P_resid, np.dot(np.linalg.inv(cov), P_resid))     #/np.float(len(P_resid))
            pinv_chi2 = np.dot(P_resid, np.dot(np.linalg.pinv(cov), P_resid))   #/np.float(len(P_resid))
            if pinv_chi2 > chi_goal: 
                if once is None: 
                    k_chi.append(0.5*(k_arr[ik-2] + k_arr[ik-1]))
                    k_w_chi.append(
                            (k_arr[ik-2] * counts[ik-2] + k_arr[ik-1] * counts[ik-1])/(counts[ik-2] + counts[ik-1]))
                    prev_inv_chi2 = np.dot(P_resid[:-1], np.dot(np.linalg.inv(cov[:-1, :-1]), P_resid[:-1])) #/np.float(len(P_resid)-1)
                    chi.append(inv_chi2)
                    chi_w.append(0.5 * (inv_chi2 + prev_inv_chi2))

                    #print cat, cov.shape
                    #print 'k_avg = ',  0.5*(k_arr[ik-2] + k_arr[ik-1])
                    #print 'Chi^2 =', inv_chi2, ' ; Chi^2 =', pinv_chi2, ' ', np.linalg.cond(cov)
                    once = 'twice'
                else: 
                    continue
        else: 
            continue
    print 'chi^2 = ', np.mean(chi)
    print 'weighted chi^2 = ', np.mean(chi_w)
    return np.mean(k_w_chi) 

if __name__=='__main__': 
    for poly_order in ['all', 'upto2', 'upto0']: 
        print poly_order
        for ell in [0, 2]: 
            print k_chi2_CM(ell, ktrust=0.3, cm_poly_order=poly_order, max_order=10, chi_goal=1., 
                    n_mocks=84, Ngrid=960, rebin=6, use='nseries')
