'''

fourier correction for P_l(k). 
Interfaces with spec modules in order to actually calculate P_l(k) for mocks/data.

'''
import numpy as np
import pickle

import pk_extrap
import fourier_corr

from util.direc import direc

def fibcoll_angularscale(redshift='median', omega_m=0.31): 
    '''
    Calculate comoving distance of fiber collision angular scale at specified redshift
    '''
    if redshift == 'median': 
        z = 0.55
    elif redshift == 'min': 
        z = 0.43
    elif redshift == 'max': 
        z = 0.7
    else: 
        raise NotImplementedError

    fibcol_ang = 62.0   # arcseconds
    fibcol_ang *= 1./3600. * 2.*np.pi/360. # radians
    
    cosmo = {} 
    cosmo['omega_M_0'] = omega_m 
    cosmo['omega_lambda_0'] = 1.0 - omega_m 
    cosmo['h'] = 0.676
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 
    
    # tranverse comoving distance
    D_M = cosmos.distance.comoving_distance_transverse(z, **cosmo) * cosmo['h']

    return D_M * fibcol_ang

def fourier_tophat_Pk(cat_corr, pk_file_name, true_pk_file_name):
    '''
    Fourier tophat corrected P(k) from true power spectra. Currently the code only corrects for the quadrupole, 
    which is our main objective. Designed to work with spec.spec.Spec.build method 

    Parameters
    ----------
    cat_corr : 
        catalog correction dictionary

    pk_file_name : 
        File name of the output power spectrum

    true_mock_name : 
        File name of the true galaxy mock catalog
    '''

    catdict = cat_corr['catalog']
    corrdict = cat_corr['correction']
    specdict = cat_corr['spec']
    
    if corrdict['name'] != 'fourier_tophat':    # (sanity check) 
        print cat_corr
        raise ValueError
        
    # read in the true power spectrum
    print true_pk_file_name
    if catdict['name'] == 'nseries': 
        tr_k, tr_p0k, tr_p2k, tr_p4k = np.loadtxt(
                    true_pk_file_name, 
                    unpack = True, 
                    usecols =[0,1,2,3] 
                    )
        tr_specs = [tr_p0k, tr_p2k, tr_p4k]
    elif catdict['name'] == 'nseriesbox': 
        tr_k, tr_p0k, tr_p2k, tr_p4k, tr_p6k = np.loadtxt(
                    true_pk_file_name, 
                    unpack = True, 
                    usecols =[0,6,2,3,4] 
                    )
        print tr_p0k
        notnorm_tr_specs = [tr_p0k, tr_p2k, tr_p4k, tr_p6k]
        tr_specs =  [(2.0*np.pi)**3 * tr_spec for tr_spec in notnorm_tr_specs]
    else: 
        raise NotImplemented

    # calculate extrapolation parameters for the true P0(k), P2(k), P4(k), ... used for the integration.
    ells = (2. * np.arange(len(tr_specs))).astype(int)
    print 'ells = ', ells
    tr_extrap_pars = [] 
    for i_ell, ellp in enumerate(ells):
        tr_extrap_pars.append(
                pk_extrap.pk_powerlaw_bestfit(tr_k, tr_specs[i_ell], k_fit=corrdict['k_fit'], k_fixed=corrdict['k_fixed'])
                )
    
    for ell in [0, 2]: 
        # calculate the correlated portion of delP_2(k)
        lps, corrdelPk_lp = fourier_corr.delP_corr(
                tr_k,               # k 
                tr_specs,           # [p0k, p2k, p4k, ..] 
                ell,                  # ell
                fs=corrdict['fs'],  # fs 
                rc=corrdict['rc'],  # rc 
                extrap_params=tr_extrap_pars, 
                k_fixed=corrdict['k_fixed'], 
                lp_comp=True
                )
    
        corrdelPk = np.zeros(len(tr_k))
        for i_lp, lp in enumerate(lps):     # save l' components to pickle file 
            corrdelPk_pickle_file = ''.join([
                '/'.join(pk_file_name.rsplit('/')[:-1]), '/',
                'corrdelP', str(ell), 'k_lp', str(lp), '_', 
                pk_file_name.rsplit('/')[-1], 
                '.p'
                ])
            pickle.dump([tr_k, corrdelPk_lp[i_lp]], open(corrdelPk_pickle_file, 'wb'))
            corrdelPk += corrdelPk_lp[i_lp]

        corrdelPk_pickle_file = ''.join([
            '/'.join(pk_file_name.rsplit('/')[:-1]), '/',
            'corrdelP', str(ell), 'k_', pk_file_name.rsplit('/')[-1], '.p'])
        pickle.dump([tr_k, corrdelPk], open(corrdelPk_pickle_file, 'wb'))

        # calculate the UNcorrelated portion of delP_2(k)
        uncorrdelPk = fourier_corr.delP_uncorr(
                tr_k,               # k 
                ell,                  # ell 
                fs=corrdict['fs'],  # fs
                rc=corrdict['rc']   # rc
                )
        uncorrdelPk_pickle_file = ''.join([
            '/'.join(pk_file_name.rsplit('/')[:-1]), '/',
            'uncorrdelP', str(ell), 'k_', pk_file_name.rsplit('/')[-1], '.p'])
        pickle.dump([tr_k, uncorrdelPk], open(uncorrdelPk_pickle_file, 'wb'))
        
        if ell == 0: 
            corr_p0k = tr_p0k + corrdelPk + uncorrdelPk
        elif ell == 2: 
            corr_p2k = tr_p2k + corrdelPk + uncorrdelPk       # corrected P2(k)
    
    if catdict['name'] == 'nseries': 
        data_list = [tr_k, corr_p0k, corr_p2k, tr_p4k]
        data_fmts = ['%10.5f', '%10.5f', '%10.5f', '%10.5f']
    elif catdict['name'] == 'nseriesbox': 
        data_list = [tr_k, corr_p0k, corr_p2k, tr_p4k, tr_p6k, tr_p0k, tr_p0k]
        data_fmts = ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f']

    np.savetxt(
            pk_file_name, 
            (np.vstack(np.array(data_list))).T, 
            fmt=data_fmts, 
            delimiter='\t'
            )
    return None

def fourier_tophat_lp_component(mock, fs=1.0, rc=0.43):
    '''
    Fourier tophat del P values calculated from the true average power spectrum.
    Hacked together. 

    Parameters
    ----------
    '''
    # calculate average true power spectrum
    if mock == 'nseries': 
        n_mock = 84
        k_fit = 0.7
        k_fixed = 0.84
        data_dir = '/mount/riachuelo1/hahn/power/Nseries/'
        for i_mock in xrange(1, n_mock+1):
            true_pk_file = ''.join([
                data_dir, 
                'POWER_Q_CutskyN', str(i_mock), '.fidcosmo.dat.grid960.P020000.box3600'])

            tr_k, tr_p0k_i, tr_p2k_i, tr_p4k_i = np.loadtxt(
                        true_pk_file, 
                        unpack = True, 
                        usecols =[0,1,2,3] 
                        )
            if i_mock == 1: 
                tr_p0k = tr_p0k_i
                tr_p2k = tr_p2k_i
                tr_p4k = tr_p4k_i
            else: 
                tr_p0k += tr_p0k_i
                tr_p2k += tr_p2k_i
                tr_p4k += tr_p4k_i 
        tr_p0k /= np.float(n_mock)
        tr_p2k /= np.float(n_mock)
        tr_p4k /= np.float(n_mock)

        tr_specs = [tr_p0k, tr_p2k, tr_p4k]

    elif mock == 'nseriesbox': 
        n_mock = 7
        k_fit = 4.3
        k_fixed = 4.34
        data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
        for i_mock in xrange(1, n_mock+1):
            true_pk_file = ''.join([
                data_dir,
                'power3600z_BoxN', str(i_mock), '.dat'])

            tr_k, tr_p0k_i, tr_p2k_i, tr_p4k_i, tr_p6k_i, tr_p8k_i, tr_p10k_i = np.loadtxt(
                        true_pk_file, 
                        unpack = True, 
                        usecols =[0,-1,2,3,4,5,6] 
                        )

            if i_mock == 1: 
                tr_p0k = tr_p0k_i
                tr_p2k = tr_p2k_i
                tr_p4k = tr_p4k_i
                tr_p6k = tr_p6k_i
                tr_p8k = tr_p8k_i
                tr_p10k = tr_p10k_i
            else: 
                tr_p0k += tr_p0k_i
                tr_p2k += tr_p2k_i
                tr_p4k += tr_p4k_i
                tr_p6k += tr_p6k_i
                tr_p8k += tr_p8k_i
                tr_p10k += tr_p10k_i

        tr_p0k /= np.float(n_mock)
        tr_p2k /= np.float(n_mock)
        tr_p4k /= np.float(n_mock)
        tr_p6k /= np.float(n_mock)
        tr_p8k /= np.float(n_mock)
        tr_p10k /= np.float(n_mock)

        notnorm_tr_specs = [tr_p0k, tr_p2k, tr_p4k, tr_p6k, tr_p8k, tr_p10k]
        tr_specs =  [(2.0*np.pi)**3 * tr_spec for tr_spec in notnorm_tr_specs]
    else: 
        raise NotImplemented

    # calculate extrapolation parameters for the 
    # true P0(k), P2(k), P4(k), etc used for the integration.
    ells = (2. * np.arange(len(tr_specs))).astype(int)
    print 'ells = ', ells
    tr_extrap_pars = [] 
    for i_ell, ellp in enumerate(ells):
        tr_extrap_pars.append(
                pk_extrap.pk_powerlaw_bestfit(tr_k, tr_specs[i_ell], k_fit=k_fit, k_fixed=k_fixed)
                )
    
    for ell in [0, 2]: # calculate the correlated portion of delP_ell(k)
        lps, corrdelPk_lp = fourier_corr.delP_corr(
                tr_k,               # k 
                tr_specs,           # [p0k, p2k, p4k, ..] 
                ell,                # ell
                fs=fs,              # fs 
                rc=rc,              # rc 
                extrap_params=tr_extrap_pars, 
                k_fixed=k_fixed, 
                lp_comp=True
                )
    
        corrdelPk = np.zeros(len(tr_k))
        for i_lp, lp in enumerate(lps):     # save l' components to pickle file 
            if mock == 'nseries': 
                avgpk_name = 'AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600'
            elif mock == 'nseriesbox':
                avgpk_name = 'AVG_power3600z_BoxN.dat'

            corrdelPk_pickle_file = ''.join([
                data_dir, 'corrdelP', str(ell), 'k_lp', str(lp), '_', avgpk_name, '.p'
                ])
            pickle.dump([tr_k, corrdelPk_lp[i_lp]], open(corrdelPk_pickle_file, 'wb'))
            corrdelPk += corrdelPk_lp[i_lp]

        corrdelPk_pickle_file = ''.join([
            data_dir, 'corrdelP', str(ell), 'k_', avgpk_name, '.p'
            ])
        pickle.dump([tr_k, corrdelPk], open(corrdelPk_pickle_file, 'wb'))

        # calculate the UNcorrelated portion of delP_2(k)
        uncorrdelPk = fourier_corr.delP_uncorr(
                tr_k,               # k 
                ell,                # ell 
                fs=fs,              # fs
                rc=rc               # rc
                )
        uncorrdelPk_pickle_file = ''.join([
            data_dir, 'uncorrdelP', str(ell), 'k_', avgpk_name, '.p'
            ])
        pickle.dump([tr_k, uncorrdelPk], open(uncorrdelPk_pickle_file, 'wb'))

    return None


if __name__=='__main__':
    fourier_tophat_lp_component('nseriesbox', fs=1.0, rc=0.43)
