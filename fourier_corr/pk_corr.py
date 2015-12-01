'''

fourier correction for P_l(k). 
Interfaces with spec modules in order to actually calculate P_l(k) for mocks/data.

'''
import pk_extrap
import cosmolopy as cosmos
from spec.spec import Spec

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

if __name__=='__main__':
    print fibcoll_angularscale(redshift='median', omega_m=0.31)
    print fibcoll_angularscale(redshift='min', omega_m=0.31)
    print fibcoll_angularscale(redshift='max', omega_m=0.31)
