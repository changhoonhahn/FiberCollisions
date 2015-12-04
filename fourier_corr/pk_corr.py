'''

fourier correction for P_l(k). 
Interfaces with spec modules in order to actually calculate P_l(k) for mocks/data.

'''
import numpy as np

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

def fourier_tophat_Pk(cat_corr, pk_file_name, true_mock_name):
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
    true_cat_corr = {
            'catalog': catdict, 
            'correction': {'name': 'true'},
            'spec': specdict
            }
    tr_pk_file = ''.join([
        direc('spec', cat_corr),    # spec directory
        'POWER_Q_',
        (true_mock_name).split('/')[-1],
        '.grid', str(specdict['Ngrid']), 
        '.P0', str(specdict['P0']), 
        '.box', str(specdict['Lbox'])
        ])

    tr_k, tr_p0k, tr_p2k, tr_p4k = np.loadtxt(
                tr_pk_file, 
                unpack = True, 
                usecols =[0,1,2,3] 
                )
    tr_specs = [tr_p0k, tr_p2k, tr_p4k]
    
    # calculate extrapolation parameters for the true P0(k), P2(k), P4(k) used for the integration.
    tr_extrap_pars = [] 
    for i_ell, ellp in enumerate([0,2,4]):
        tr_extrap_pars.append(
                pk_extrap.pk_powerlaw_bestfit(tr_k, tr_specs[i_ell], k_fit=corrdict['k_fit'], k_fixed=corrdict['k_fixed'])
                )

    # calculate the correlated portion of delP_2(k)
    corrdelP2k = fourier_corr.delP_corr(
            tr_k,               # k 
            tr_specs,           # [p0k, p2k, p4k] 
            2,                  # ell
            fs=corrdict['fs'],  # fs 
            rc=corrdict['rc'],  # rc 
            extrap_params=tr_extrap_pars, 
            k_fixed=corrdict['k_fixed']
            )

    # calculate the UNcorrelated portion of delP_2(k)
    uncorrdelP2k = fourier_corr.delP_uncorr(
            tr_k,               # k 
            2,                  # ell 
            fs=corrdict['fs'],  # fs
            rc=corrdict['rc']   # rc
            )
    
    corr_p2k = tr_p2k + corrdelP2k + uncorrdelP2k       # corrected P2(k)

    data_list = [tr_k, tr_p0k, corr_p2k, tr_p4k]
    data_fmts = ['%10.5f', '%10.5f', '%10.5f', '%10.5f']

    np.savetxt(
            pk_file_name, 
            (np.vstack(np.array(data_list))).T, 
            fmt=data_fmts, 
            delimiter='\t'
            )
    return None

"""
if __name__=='__main__':
    print fibcoll_angularscale(redshift='median', omega_m=0.31)
    print fibcoll_angularscale(redshift='min', omega_m=0.31)
    print fibcoll_angularscale(redshift='max', omega_m=0.31)
"""