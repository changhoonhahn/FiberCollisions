'''

fourier correction for P_l(k). 
Interfaces with spec modules in order to actually calculate P_l(k) for mocks/data.

'''
import codif
import numpy as np
import pickle
import cosmolopy as cosmos
from scipy.interpolate import interp1d

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

def fourier_tophat_lp_component(mock, fs=1.0, rc=0.43, noextrap=''):
    '''
    Fourier tophat del P values calculated from the true average power spectrum.
    Hacked together. Writes out delP values to pickle files. 

    Calculates delP_0^corr and del P_2^corr

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
        # default parameters for nseries box 
        n_mock = 7
        k_fit = 4.
        k_fixed = 4.34
    
        # loop through mock realizations
        data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
        for i_mock in xrange(1, n_mock+1):
            true_pk_file = ''.join([data_dir,
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
    # true P0(k), P2(k), P4(k), ... etc used for the integration.
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
                lp_comp=True, 
                noextrap=noextrap
                )
    
        corrdelPk = np.zeros(len(tr_k))
        for i_lp, lp in enumerate(lps):     # save l' components to pickle file 
            if mock == 'nseries': 
                avgpk_name = 'AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600'
            elif mock == 'nseriesbox':
                avgpk_name = 'AVG_power3600z_BoxN.dat'

            corrdelPk_pickle_file = ''.join([
                data_dir, 'corrdelP', str(ell), 'k_lp', str(lp), noextrap, '_', avgpk_name, '.p'
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

def d_DelPbox_corr_lp(ell, fs=1.0, rc=0.43):
    ''' Calculate the correction to Del P^corr_box from the difference
    in P_l(k) of the Nseries to the Nseries box. Each l' components are 
    written out to pickle files. 
    
    This calculation is almost certaintly an underestimate because: 
        - only the l' = 0, 2, 4 contributions can be calculated 
        because l' = 6, 8 cannot be calculated for the actual survey
        geometry. 
        - the dq integration required only spans the k range of the 
        Nseries P_l(k) (roughly 0.001 to  0.9). So contributions 
        at smaller scales are ignored. 

    Calculates d(delP_0^corr) and d(del P_2^corr)

    Parameters
    ----------
    ell : int
        Int specifying the ell value (l = 0, 2)

    fs : float
        Float that specifies the survey fraction exposed to fiber collisions. 
        By definition fs < 1.0. For SDSS tiling fs~1.0

    rc : float
        Float that specifies the fiber collision comoving scale 
    '''
    # calculate d(P_l(k)) between Nseries and Nseries box 

    nseries_file = lambda i: ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/', 
        'POWER_Q_CutskyN', str(i), '.fidcosmo.dat.grid960.P020000.box3600'])
    box_file = lambda i: ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'power3600z_BoxN', str(i), '.dat'])
    
    for mock in ['nseries', 'nseriesbox']: 
        if mock == 'nseries': 
            n_mock = 20 
            file_func = nseries_file
        elif mock == 'nseriesbox': 
            n_mock = 7 
            file_func = box_file
    
        tr_specs = []
        for el in [0,2,4]: 
            for i_mock in xrange(1, n_mock+1):
                true_file = file_func(i_mock)

                l_index = el/2+1
                if mock == 'nseriesbox' and el == 0: 
                    l_index = -1

                tr_k, tr_plk_i = np.loadtxt(
                            true_file, 
                            unpack = True, 
                            usecols =[0,l_index] 
                            )
                if i_mock == 1: 
                    tr_plk = tr_plk_i
                else: 
                    tr_plk += tr_plk_i
            tr_plk /= np.float(n_mock)
            tr_specs.append(tr_plk)
        
        if mock == 'nseries': 
            nseries_k = tr_k 
            nseries_specs = tr_specs
        elif mock == 'nseriesbox': 
            box_k = tr_k 
            box_specs = [(2.0*np.pi)**3 * tr_spec for tr_spec in tr_specs]
    
    # dP_l(k) = P^nseries_l(k) - P^box_l(k)
    dP_specs = [] 
    k_range = np.where((nseries_k > box_k[0]) & (nseries_k < box_k[-1]))
    k_val = nseries_k[k_range]
    for ii in range(len(nseries_specs)): 
        box_interp = interp1d(box_k, box_specs[ii])
        dP_specs.append((nseries_specs[ii])[k_range] - box_interp(k_val))

    lps, d_DelP_lp = fourier_corr.delP_corr_qmax(
            k_val, 
            dP_specs,           # [p0k, p2k, p4k, ..] 
            ell,                # ell
            q_max=k_val[-1], 
            fs=fs,              # fs 
            rc=rc,              # rc 
            lp_comp=True
            )

    d_DelP_lp_file = lambda lp: ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'd_DelP', str(ell), 'k_boxcorr', '.lp', str(lp), 
        '.fs', str(round(fs, 2)), 
        '.rc', str(round(rc, 2)), 
        '.p'])

    d_DelP = np.zeros(len(k_val))
    for i_lp, lp in enumerate(lps):     # save l' components to pickle file
        pickle.dump([k_val, d_DelP_lp[i_lp]], open(d_DelP_lp_file(lp), 'wb'))
        d_DelP += d_DelP_lp[i_lp]

    d_DelP_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'd_DelP', str(ell), 'k_boxcorr', 
        '.fs', str(round(fs, 2)), 
        '.rc', str(round(rc, 2)), 
        '.p'])
    pickle.dump([k_val, d_DelP], open(d_DelP_file, 'wb'))

    return None

def fourier_tophat_lp_component_qmaxtest(mock, q_max, fs=1.0, rc=0.43):
    '''
    Calculate Fourier tophat del P values from the true average power spectrum
    for specified q_max value. This is the next step test after 'no extrapolation'
    in order to determine what q_max value is necessary to calculate a delP, within
    sample variance. 

    delP_uncorr is constant
    delP_corr = \sum_lp \int\limits_0^{q_max}
    Writes out delP values to pickle files. 

    Calculates delP_0^corr and del P_2^corr

    Parameters
    ----------
    '''
    # calculate average true power spectrum
    if mock == 'nseries': 
        n_mock = 84
        k_fit = 0.7
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
        # default parameters for nseries box 
        n_mock = 7
        k_fit = 4.
    
        # loop through mock realizations
        data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
        for i_mock in xrange(1, n_mock+1):
            true_pk_file = ''.join([data_dir,
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

    if q_max > tr_k[-1]: 
        raise ValueError("For this test, q_max can't be greater than max(k)")
    
    for ell in [0, 2]: # calculate the correlated portion of delP_ell(k)
        lps, corrdelPk_lp = fourier_corr.delP_corr_qmax(
                tr_k,               # k 
                tr_specs,           # [p0k, p2k, p4k, ..] 
                ell,                # ell
                q_max=q_max, 
                fs=fs,              # fs 
                rc=rc,              # rc 
                lp_comp=True
                )
    
        corrdelPk = np.zeros(len(tr_k))
        for i_lp, lp in enumerate(lps):     # save l' components to pickle file 
            if mock == 'nseries': 
                avgpk_name = 'AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600'
            elif mock == 'nseriesbox':
                avgpk_name = 'AVG_power3600z_BoxN.dat'

            corrdelPk_pickle_file = ''.join([
                data_dir, 'corrdelP', str(ell), 'k_lp', str(lp), 
                '_qmax', str(round(q_max,2)), '_', avgpk_name, '.p'
                ])
            pickle.dump([tr_k, corrdelPk_lp[i_lp]], open(corrdelPk_pickle_file, 'wb'))
            corrdelPk += corrdelPk_lp[i_lp]

        corrdelPk_pickle_file = ''.join([
            data_dir, 'corrdelP', str(ell), 'k', 
            '_qmax', str(round(q_max,2)),'_',avgpk_name, '.p'
            ])
        pickle.dump([tr_k, corrdelPk], open(corrdelPk_pickle_file, 'wb'))

    return None

def fourier_tophat_lp_component_ktrust_qmax(mock, ktrust, fs=1.0, rc=0.43):
    '''
    Calculate Fourier tophat Del P^corr values, integrated from q=ktrust to q=qmax 
    using the true average power spectrum. The calculated Del P^corr |ktrust to qmax
    is written out to Pickle Files
    '''
    # calculate average true power spectrum
    if mock == 'nseries': 
        n_mock = 84
        k_fit = 0.7
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
        # default parameters for nseries box 
        n_mock = 7
        # loop through mock realizations
        data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
        for i_mock in xrange(1, n_mock+1):
            true_pk_file = ''.join([data_dir,
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

        notnorm_tr_specs = [tr_p0k, tr_p2k, tr_p4k, tr_p6k, tr_p8k, tr_p10k]
        tr_specs =  [(2.0*np.pi)**3 * tr_spec / np.float(n_mock) for tr_spec in notnorm_tr_specs]
    else: 
        raise NotImplemented
    
    q_max = tr_k[-1]
    if ktrust > q_max: 
        raise ValueError 
    
    for ell in [0, 2]: 
        # calculate Del P^corr integrated over q = ktrust to qmax 
        lps, corrdelPk_lp = fourier_corr.delP_corr_qrange(
                tr_k,               # k 
                tr_specs,           # [p0k, p2k, p4k, ..] 
                ell,                # ell
                q_min=ktrust,       # k_trust
                q_max=q_max,        #
                fs=fs,              # fs 
                rc=rc,              # rc 
                lp_comp=True
                )
    
        corrdelPk = np.zeros(len(tr_k))
        for i_lp, lp in enumerate(lps):     # save l' components to pickle file 
            if mock == 'nseries': 
                avgpk_name = 'AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600'
            elif mock == 'nseriesbox':
                avgpk_name = 'AVG_power3600z_BoxN.dat'

            corrdelPk_pickle_file = ''.join([data_dir, 
                'corrdelP', str(ell), 'k_lp', str(lp), 
                '.ktrust', str(round(ktrust, 2)), '_qmax.', avgpk_name, '.p'])
            pickle.dump([tr_k, corrdelPk_lp[i_lp]], open(corrdelPk_pickle_file, 'wb'))
            corrdelPk += corrdelPk_lp[i_lp]

        corrdelPk_pickle_file = ''.join([data_dir, 
            'corrdelP', str(ell), 'k', 
            '.ktrust', str(round(ktrust, 2)), '_qmax.', avgpk_name, '.p'])
        pickle.dump([tr_k, corrdelPk], open(corrdelPk_pickle_file, 'wb'))
        
    codif.notif(subject='Del P^corr |ktrust='+str(round(ktrust,2))+' to qmax integration')
    return None

def Pk_tophat_conv(cat_corr, pk_file_name, true_pk_file_name):
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



if __name__=='__main__': 
    d_DelPbox_corr_lp(0, fs=0.7, rc=0.43)
    d_DelPbox_corr_lp(2, fs=0.7, rc=0.43)

    #for ktrust in [0.3, 0.4, 0.5, 0.7, 1.0, 3.0]:
    #    fourier_tophat_lp_component_ktrust_qmax('nseriesbox', ktrust, fs=1.0, rc=0.43)
    #for qmax in [0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 2.0, 3.0, 4.0]: 
    #    print 'q_max = ', qmax
    #    fourier_tophat_lp_component_qmaxtest('nseriesbox', qmax, fs=1.0, rc=0.43)
    #print fibcoll_angularscale(redshift='min', omega_m=0.31) 
    #print fibcoll_angularscale(redshift='median', omega_m=0.31)
    #print fibcoll_angularscale(redshift='max', omega_m=0.31)

    #fourier_tophat_lp_component('nseriesbox', fs=1.0, rc=0.43, noextrap='_noextrap')
