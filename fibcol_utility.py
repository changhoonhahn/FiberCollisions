'''


Utility codes for FiberCollisions project


Author(s): ChangHoon Hahn 


'''



import os.path
import numpy as np
import subprocess
import cosmolopy as cosmos
import fibcol_nbar as fc_nbar
import fibcol_data as fc_data
import fibcol_spec as fc_spec

def fortran_code(fft_power, **cat_corr): 
    ''' Return appropriate FORTRAN code for calculating FFT or powerspectrum

    Parameters
    ----------
    fft_power : 'fft', 'power', 'quadfft', 'quadpower'
    cat_corr : catalog correction dictionary 


    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    spec = cat_corr['spec']
    
    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ----------------------

        # code directory 
        ldg_code_dir = '/home/users/hahn/powercode/FiberCollisions/LasDamas/Geo/' 
       
        if fft_power.lower() == 'fft':                      # FFT -----------------

            # literature
            if correction['name'].lower() == 'floriansn':     # Beutler+2014
                f_code = ldg_code_dir+'FFT_ldg_fkp_w_florian_'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn':  # Gil-Marin+2014
                f_code = ldg_code_dir+'FFT_ldg_fkp_w_hector_'+str(spec['grid'])+'grid.f'
            else: 
                f_code = ldg_code_dir+'FFT_ldg_fkp_w_'+str(spec['grid'])+'grid.f'

        elif fft_power.lower() == 'power':                  # power ----------------

            if correction['name'].lower() in ('true', 'upweight', 'peaknbar'):  
                # FKP estimator
                if spec['grid'] == 360: 
                    f_code = ldg_code_dir+'power_ldg_fkp_360grid_180bin.f'
                elif spec['grid'] == 960: 
                    f_code = ldg_code_dir+'power_ldg_fkp_960grid_480bin.f'

            elif correction['name'].lower() in ('peakshot', 'allpeakshot', 'noweight', 
                    'shotnoise', 'vlospeakshot', 'floriansn', 'hectorsn', 'peakshot_dnn'): 
                # Igal+Irand shot noise incorporated
                if spec['grid'] == 360: 
                    f_code = ldg_code_dir+'power_ldg_fkp_Igal_Iran_360grid_180bin.f'
                elif spec['grid'] == 960: 
                    f_code = ldg_code_dir+'power_ldg_fkp_Igal_Iran_960grid_480bin.f'
            else: 
                raise NameError('what?')
    
        # quadrupole codes ------------------------------------------------
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError('asdflkajsdf') 

    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ----------------

        # code directory 
        code_dir = '/home/users/hahn/powercode/FiberCollisions/TilingMock/'
    
        if fft_power.lower() == 'fft':
            if correction['name'].lower() == 'floriansn': 
                f_code = code_dir+'FFT-fkp-tm-w-nbar-florian-'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn': 
                f_code = code_dir+'FFT-fkp-tm-w-nbar-hector-'+str(spec['grid'])+'grid.f'
            else: 
                f_code = code_dir+'FFT-fkp-tm-w-nbar-'+str(spec['grid'])+'grid.f'

        elif fft_power.lower() == 'power': 
            if correction['name'].lower() in ('peakshot', 'allpeakshot', 'shotnoise', 'floriansn', 'hectorsn', 'vlospeakshot', 'peakshot_dnn'):
                if spec['grid'] == 360: 
                    f_code = code_dir+'power-fkp-tm-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    f_code = code_dir+'power-fkp-tm-w-nbar-Igal-Irand-960grid-480bin.f'
            else: 
                if spec['grid'] ==360: 
                    f_code = code_dir+'power-fkp-tm-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960:
                    f_code = code_dir+'power-fkp-tm-w-nbar-960grid-480bin.f'

        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError('asdlkfjalksdjfklasjf')

    elif catalog['name'].lower() == 'qpm':                      # QPM -----------------------

        code_dir = '/home/users/hahn/powercode/FiberCollisions/QPM/dr12d/'
        
        if fft_power.lower() == 'fft': 
            # FFT code
            if correction['name'].lower() == 'floriansn': 
                f_code = code_dir+'FFT-qpm-fkp-w-nbar-florian-'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn': 
                f_code = code_dir+'FFT-qpm-fkp-w-nbar-hector-'+str(spec['grid'])+'grid.f'
            else: 
                f_code = code_dir+'FFT-qpm-fkp-w-nbar-'+str(spec['grid'])+'grid.f'
    
        elif fft_power.lower() == 'power': 
            if correction['name'].lower() in ('true', 'upweight', 'peaknbar'):
                # normal FKP shot noise correction
                if spec['grid'] == 360: 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-960grid-480bin.f'

            elif correction['name'].lower() in \
                    ('peakshot', 'shotnoise', 'floriansn', 'noweight', 
                            'hectorsn', 'vlospeakshot', 'peakshot_dnn'): 
                if spec['grid'] == 360: 
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-Igal-Irand-960grid-480bin.f'

        # quadrupole codes --------------------------------------------
        # regardess of catalog or correction TEMPORARILY HARDCODED HERE FOR TEST RUN 
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError("not Yet coded") 
    
    elif catalog['name'].lower() == 'nseries':                  # N series ------------------

        code_dir = 'Nseries/'
        
        if fft_power.lower() == 'fft':          # FFT

            if correction['name'].lower() == 'floriansn': 
                f_code = code_dir+'FFT-nseries-fkp-w-nbar-florian-'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn': 
                f_code = code_dir+'FFT-nseries-fkp-w-nbar-hector-'+str(spec['grid'])+'grid.f'
            else: 
                f_code = code_dir+'FFT-nseries-fkp-w-nbar-'+str(spec['grid'])+'grid.f'
    
        elif fft_power.lower() == 'power': 
            if correction['name'].lower() in ('true', 'upweight', 'peaknbar'):
                # normal FKP shot noise correction
                if spec['grid'] == 360: 
                    f_code = code_dir+'power-nseries-fkp-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    f_code = code_dir+'power-nseries-fkp-w-nbar-960grid-480bin.f'

            elif correction['name'].lower() in \
                    ('peakshot', 'shotnoise', 'floriansn', 'noweight', 'hectorsn', 'peakshot_dnn'): 
                if spec['grid'] == 360: 
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-nseries-fkp-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-nseries-fkp-w-nbar-Igal-Irand-960grid-480bin.f'

            elif 'scratch' in correction['name'].lower(): 
                if spec['grid'] == 360: 
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-nseries-fkp-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-nseries-fkp-w-nbar-Igal-Irand-960grid-480bin.f'

        # quadrupole codes --------------------------------------------
        # regardess of catalog or correction TEMPORARILY HARDCODED HERE FOR TEST RUN 
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError("not Yet coded") 

    elif catalog['name'].lower() == 'patchy':                       # PATCHY --------------------

        code_dir = '/home/users/hahn/powercode/FiberCollisions/PATCHY/dr12/v6c/'
        
        if fft_power.lower() == 'fft': 
            # FFT code
            if correction['name'].lower() == 'floriansn': 
                f_code = code_dir+'FFT-patchy-fkp-w-nbar-florian-'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn': 
                f_code = code_dir+'FFT-patchy-fkp-w-nbar-hector-'+str(spec['grid'])+'grid.f'
            else: 
                f_code = code_dir+'FFT-patchy-fkp-w-nbar-'+str(spec['grid'])+'grid.f'
    
        elif fft_power.lower() == 'power': 
            if correction['name'].lower() in ('true', 'upweight', 'peaknbar'):
                # normal FKP shot noise correction
                if spec['grid'] == 360: 
                    f_code = code_dir+'power-patchy-fkp-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    pass
                    f_code = code_dir+'power-qpm-fkp-w-nbar-960grid-480bin.f'

            elif correction['name'].lower() in \
                    ('peakshot', 'shotnoise', 'floriansn', 
                            'hectorsn', 'vlospeakshot', 'peakshot_dnn'): 
                if spec['grid'] == 360: 
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    pass
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-Igal-Irand-960grid-480bin.f'

        # quadrupole codes --------------------------------------------
        # regardess of catalog or correction TEMPORARILY HARDCODED HERE FOR TEST RUN 
        '''
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError("not Yet coded") 
        '''
    
    else: 
        raise NaemError('Not coded!') 

    return f_code

def fortran_code2exe(code): 
    '''
    get .exe file based on fortran code file name  
    '''
    code_dir = '/'.join(code.split('/')[0:-1])+'/' 
    code_file = code.split('/')[-1]

    fort_exe = code_dir+'exe/'+'.'.join(code_file.rsplit('.')[0:-1])+'.exe'

    return fort_exe

def compile_fortran_code(code): 
    '''
    compiles fortran code (very simple, may not work) 
    '''
    # get executable file 
    fort_exe = fortran_code2exe(code) 

    # compile command
    if code == '/home/users/hahn/powercode/FiberCollisions/FFT_FKP_BOSS_cic_il4_v3.f': 
        compile_cmd = ' '.join(['ifort -fast -o', fort_exe, code, '-L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw -lm'])
    elif code == '/home/users/hahn/powercode/FiberCollisions/power_FKP_SDSS_BOSS_v3.f': 
        compile_cmd = ' '.join(['ifort -fast -o', fort_exe, code])
    else: 
        compile_cmd = ' '.join(['ifort -O3 -o', fort_exe, code, '-L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw'])
    print compile_cmd

    # call compile command 
    subprocess.call(compile_cmd.split())

def get_fibcoll_dir(file_type, **cat_corr): 
    '''
    get data/FFT/power directories given catalog
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    if file_type.lower() not in ('data', 'fft', 'power'): 
        raise NameError('either data, fft, or power') 

    else: 
        if catalog['name'].lower() == 'lasdamasgeo': 
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

        else: 
            raise NameError('not yet coded')
    return file_dir 

def get_fig_dir(): 
    '''
    return figure directory for fibcollision project
    '''
    return '/home/users/hahn/research/figures/boss/fiber_collision/'

def find_nearest(array,value, index=False):
    '''
    Find nearest element in array to value. If index is True then return index
    '''
    idx = (np.abs(array-value)).argmin()
    if index == False: 
        return array[idx]
    else: 
        return idx

def fibcol_file_update(): 
    '''
    hardcoded routine in order to update file names as I use more correction schemes without re-running 
    '''
    
    for file_type in ['data', 'fft', 'power']:
        ###### change names of files
        ###### change names of files
        ###### change names of files
        ###### change names of files
        ###### change names of files
        ###### change names of files
        ###### change names of files
        ######
        ######
        ######
        ######
        ######
        ######
        ######
        ######
        pass

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

