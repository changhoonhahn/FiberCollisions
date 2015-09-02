'''


Utility codes for FiberCollisions project


Author(s): ChangHoon Hahn 


'''
import cosmolopy as cosmos

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

def get_fibcoll_dir(file_type, **cat_corr): 
    '''
    get data/FFT/power directories given catalog
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    if file_type.lower() not in ('data', 'fft', 'power'): 
        raise NameError('either data, fft, or power') 

    else: 
        if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'): 
            # Lasdamasgeo -----------------------------------------
            
            if file_type.lower() == 'data': 
                file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
            elif file_type.lower() == 'fft': 
                file_dir = '/mount/riachuelo1/hahn/FFT/LasDamas/Geo/'
            else:
                file_dir = '/mount/riachuelo1/hahn/power/LasDamas/Geo/'

        elif catalog['name'].lower() == 'tilingmock': 
            # Tiling mock -------------------------------------------

            if file_type.lower() == 'data': 
                file_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'
            elif file_type.lower() == 'fft': 
                file_dir = '/mount/riachuelo1/hahn/FFT/tiling_mocks/'
            else:
                file_dir = '/mount/riachuelo1/hahn/power/tiling_mocks/'

        elif catalog['name'].lower() == 'qpm':                          # QPM ---------------------------------------------

            if file_type.lower() == 'data': 
                file_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'
            elif file_type.lower() == 'fft': 
                file_dir = '/mount/riachuelo1/hahn/FFT/QPM/dr12d/'
            else:
                file_dir = '/mount/riachuelo1/hahn/power/QPM/dr12d/'

        elif catalog['name'].lower() == 'nseries':                          # N series ---------------------------------------

            if file_type.lower() == 'data': 
                file_dir = '/mount/riachuelo1/hahn/data/Nseries/'
            elif file_type.lower() == 'fft': 
                file_dir = '/mount/riachuelo1/hahn/FFT/Nseries/'
            else:
                file_dir = '/mount/riachuelo1/hahn/power/Nseries/'

        elif catalog['name'].lower() == 'patchy':                       # PATCHY ----------------------------------------
            
            if file_type.lower() == 'data': 
                file_dir = '/mount/riachuelo1/hahn/data/PATCHY/dr12/v6c/'
            elif file_type.lower() == 'fft': 
                file_dir = '/mount/riachuelo1/hahn/FFT/PATCHY/dr12/v6c/'
            else:
                file_dir = '/mount/riachuelo1/hahn/power/PATCHY/dr12/v6c/'

        elif 'bigmd' in catalog['name'].lower():                # Big MD --------------------
            if file_type.lower() == 'data': 
                file_dir = '/mount/riachuelo1/hahn/data/BigMD/'
            elif file_type.lower() == 'fft': 
                file_dir = '/mount/riachuelo1/hahn/FFT/BigMD/'
            else:
                file_dir = '/mount/riachuelo1/hahn/power/BigMD/'

        elif catalog['name'].lower() == 'cmass':                # CMASS ----------------------
            if file_type.lower() == 'data': 
                file_dir = '/mount/riachuelo1/hahn/data/CMASS/'
            elif file_type.lower() == 'fft': 
                file_dir = '/mount/riachuelo1/hahn/FFT/CMASS/'
            else:
                file_dir = '/mount/riachuelo1/hahn/power/CMASS/'
        else: 
            raise NameError('not yet coded')
    return file_dir 

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
