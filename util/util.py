'''


Utility codes for FiberCollisions project


Author(s): ChangHoon Hahn 


'''
import cosmolopy as cosmos
import numpy as np 
import scipy as sp

# --- Local ---

def comdis2z(comdis, **cosmo): 
    ''' Given comoving distance and cosmology, determine z 
    using cubic spline

    Notes
    -----
    * Comoving distance *has* to be in Mpc/h
    '''
    z_arr = np.array([0.0+0.05*np.float(i) for i in range(21)]) 
    dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo)*cosmo['h']

    dmz_spline = sp.interpolate.interp1d(dm_arr, z_arr, kind='cubic') 
    
    z = dmz_spline(comdis)

    return z 

def get_fig_dir(): 
    '''
    return figure directory for fibcollision project
    '''
    return '/home/users/hahn/research/figures/boss/fiber_collision/'

def find_nearest(array, value, index=False):
    '''
    Find nearest element in array to value. If index is True then return index
    '''
    idx = (np.abs(array-value)).argmin()
    if index == False: 
        return array[idx]
    else: 
        return idx

def radecz_to_xyz(ra, dec, z, **cosmo):
    ''' Given RA, Dec, redshift AND cosmology, calculate x,y,z in Mpc/h
    '''
    phi = ra 
    theta = 90.0 - dec 
    r = cosmos.distance.comoving_distance(z, **cosmo)*cosmo['h']    # Mpc/h

    x = r * np.cos(np.deg2rad(phi)) * np.sin(np.deg2rad(theta)) 
    y = r * np.sin(np.deg2rad(phi)) * np.sin(np.deg2rad(theta))
    z = r * np.cos(np.deg2rad(theta))
    return (x,y,z)

def ang_sep(ra1, dec1, ra2, dec2): 
    ''' Given a pair of ra and decs in DEGREES gives angular separation in DEGREES
    '''
    # convert to radians 
    ra1 = ra1*np.pi/180.
    dec1 = dec1*np.pi/180.
    ra2 = ra2*np.pi/180.
    dec2 = dec2*np.pi/180.

    x = np.cos(ra1)*np.cos(dec1)*np.cos(ra2)*np.cos(dec2) 
    y = np.sin(ra1)*np.cos(dec1)*np.sin(ra2)*np.cos(dec2) 
    z = np.sin(dec1)*np.sin(dec2)

    rad = np.arccos(x+y+z)
    
    sep = rad
    #sep = np.choose( rad<0.000004848 , ( np.sqrt( (np.cos(dec1)*(ra1-ra2))**2+(dec1-dec2)**2), rad))

    sep = sep*180./np.pi
    return sep
