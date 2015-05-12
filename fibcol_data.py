'''

Deal with galaxy data for FiberCollisions 

Author(s): ChangHoon Hahn 


'''

import numpy as np
import time 
import os.path
import subprocess
import cosmolopy as cosmos
import warnings 
import matplotlib.pyplot as plt

# --- Local ----
import fibcol_dlos as fc_dlos
import galaxy_environment as genv
import pyspherematch as pysph

# Classes ------------------------------------------------------------
class galaxy_data: 
    def __init__(self, DorR, clobber=False, **cat_corr): 
        ''' Given cat_corr dictionary read in/ build the data or random 
        file and store all the appropriate values  

        Parameters
        ----------
        DorR : 'data' or 'random'
        cat_corr = Catalog correction Dictionary 

        '''
        self.cat_corr = cat_corr                    # save catalog and correction metadata 
        catalog = cat_corr['catalog'] 
        correction = cat_corr['correction'] 

        file_name = get_galaxy_data_file(DorR, **cat_corr)              # get file name 
        self.file_name = file_name          # save file name to class 

        if catalog['name'].lower() == 'lasdamasgeo':    # LasDamas Geo -----------------------

            omega_m = 0.25  # cosmology

            if DorR == 'data':                          # Data -------------------------
               
                # columns that this catalog data will have  
                catalog_columns = ['ra', 'dec', 'z', 'weight']  

                self.columns = catalog_columns

                if (os.path.isfile(file_name) == False) or (clobber == True): 
                    # if file does not exists, make file  
                    print 'Building', file_name

                    if correction['name'].lower() == 'true': # True (No FiberCollisions) 
                        build_true(**cat_corr)
                
                    elif correction['name'].lower() in ('upweight', 'fibcol', 
                            'shotnoise', 'hectorsn', 'floriansn'): 
                        # Fibercollisions Corrected by Upweight correction  
                        build_fibercollided(**cat_corr)
                
                    elif correction['name'].lower() in ('peak', 'peaknbar', 
                            'peakshot', 'peaktest'): 
                        # Fibercollisions Corrected by peak correction 
                        build_peakcorrected_fibcol(sanitycheck=True, **cat_corr)
                        
                    elif correction['name'].lower() in ('peakshot_dnn'): 
                        ''' Peak Correction (with dNN env) + Shotnoise 
                        Correction needs to specify: fit, nth NN, and sigma 
                        '''
                        build_peak_fpeak_dNN( NN=correction['NN'], **cat_corr) 

                    elif correction['name'].lower() in ('noweight'): 
                        # No weight
                        build_noweight(**cat_corr) 

                    else: 
                        raise NameError('Correction Name Unknown') 

                # read data (ra, dec, z, weights) 
                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         

                for i_col, catalog_column in enumerate(catalog_columns): 
                    column_data = file_data[i_col]
                    setattr(self, catalog_column, column_data) 

            elif DorR.lower() == 'random':          # Random Catalogs -------------------------

                # columns of random catalog (NOTE CZ IS CONVERTED TO Z) 
                catalog_columns = ['ra', 'dec', 'z']        

                self.columns = catalog_columns

                # Read data 
                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])         

                for i_col, catalog_column in enumerate(catalog_columns): 
                    if catalog_column == 'z': 
                        column_data = file_data[i_col]/299800.0
                    else: 
                        column_data = file_data[i_col]
                    # assign to class
                    setattr(self, catalog_column, column_data)

        elif catalog['name'].lower() == 'tilingmock':   
            # TILING MOCKS ------------------------------------------------ 

            omega_m = 0.274         # survey cosmology 

            if DorR == 'data':              # for mocks 

                catalog_columns = ['ra', 'dec', 'z', 'weight']       # columns that this catalog data will have  
                self.columns = catalog_columns

                if (os.path.isfile(file_name) == False) or (clobber == True): # if file does not exists, make file  
                    print 'Constructing ', file_name 

                    if correction['name'].lower() == 'true':    # true mocks
                        # all weights = 1 (fibercollisions *not* imposed) 
                        build_true(**cat_corr) 
                
                    elif correction['name'].lower() in ('upweight', 'fibcol', 'shotnoise', 'floriansn', 'hectorsn'):
                        # upweighted mocks
                        build_fibercollided(**cat_corr)         # build-in fibercollisions using spherematch idl code
                
                    elif correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot'): 
                        # Correction methods that have both peak and tail contributions   
                        # peak/peaknbar = peak + tail correction 
                        # peaktest = fpeak peak correction + remove rest 
                        # peakshot = fpeak peak correction + shot noise for rest
                        build_peakcorrected_fibcol(sanitycheck=True, **cat_corr)
                
                    elif correction['name'].lower() in ('allpeak', 'allpeakshot'):
                        # all peak corrected 
                        build_peakcorrected_fibcol(sanitycheck=True, **cat_corr)

                # Read data  
                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weight

                # assign to data columns class
                for i_col, catalog_column in enumerate(catalog_columns): 
                    column_data = file_data[i_col]
                    # assign to class
                    setattr(self, catalog_column, column_data) 

            elif DorR.lower() == 'random':              # Randoms -----------------------------------------------------

                catalog_columns = ['ra', 'dec', 'z']        # columns of catalog 
                self.columns = catalog_columns

                if (os.path.isfile(file_name) == False) or (clobber == True): # if file does not exists, make file  
                    print 'Constructing ', file_name
                    build_corrected_randoms(sanitycheck=False, **cat_corr)       # impose redshift limit

                print 'Reading ', file_name                     # just because it takes forever
                t0 = time.time() 
                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])      # read random file
                print 'took ', (time.time()-t0)/60.0, ' mins'       # print time 

                for i_col, catalog_column in enumerate(catalog_columns): 
                    column_data = file_data[i_col]
                    # assign data column to class
                    setattr(self, catalog_column, column_data)
                
        elif catalog['name'].lower() == 'qpm':          # QPM -------------------------------

            omega_m = 0.31  # survey cosmology 

            if DorR == 'data':                          # Data ------------------------------
                
                # catalog columns 
                catalog_columns = ['ra', 'dec', 'z', 'wfkp', 'wfc', 'comp'] 
                self.columns = catalog_columns

                if (os.path.isfile(file_name) == False) or (clobber == True):
                    # File does not exist or Clobber = True!

                    print 'Constructing ', file_name 
                    if correction['name'].lower() == 'true':                    
                        # true mocks 
                        # all weights = 1 (fibercollisions *not* imposed) 
                        build_true(**cat_corr) 

                    elif correction['name'].lower() in ('upweight', 'shotnoise', 'floriansn', 'hectorsn'): 
                        # upweighted mocks
                        build_fibercollided(**cat_corr) 

                    elif correction['name'].lower() in ('peaknbar', 'peakshot'): 
                        # peak corrected mocks 
                        build_peakcorrected_fibcol(**cat_corr)  # build peak corrected file 

                    elif correction['name'].lower() in ('peakshot_dnn'):
                        # peak + dLOS env correct mocks 
                        build_peak_fpeak_dNN(NN=correction['NN'], **cat_corr) 

                    elif correction['name'].lower() in ('tailupw'):         
                        # tail upweight correction 
                        build_tailupweight_fibcol(**cat_corr)  # build peak corrected file 

                    elif correction['name'].lower() in ('noweight'): 
                        # n oweight 
                        build_noweight(**cat_corr) 
                    else: 
                        raise NotImplementedError() 

                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5])         # ra, dec, z, wfkp, wfc, comp

                # assign to data columns class
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col]) 
                        
            
            elif DorR == 'random':              # Random ------------------------------------

                catalog_columns = ['ra', 'dec', 'z', 'wfkp']    # catalog columns 

                if (os.path.isfile(file_name) == False) or (clobber == True):
                    print 'Constructing ', file_name 
                    build_random(**cat_corr) 

                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])             # ra, dec, z, wfkp

                # assign to object data columns
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col])
        
        elif catalog['name'].lower() == 'patchy': 
            # PATCHY Mocks --------------------------------------------------
            omega_m = 0.31              # survey cosmology 

            if DorR == 'data':                      
                # mocks ------------------------------------------------------

                catalog_columns = ['ra', 'dec', 'z', 'nbar', 'wfc']         # catalog columns 
                self.columns = catalog_columns

                if (os.path.isfile(file_name) == False) or (clobber == True):
                    # File does not exist or Clobber = True!
                    print 'Constructing ', file_name 

                    if correction['name'].lower() == 'true': 
                        # true mocks 
                        build_true(**cat_corr) 

                    elif correction['name'].lower() in ('upweight', 'shotnoise', 
                            'floriansn', 'hectorsn'): 
                        # upweighted mocks
                        build_fibercollided(**cat_corr) 

                    elif correction['name'].lower() in ('peaknbar', 'peakshot'): 
                        # peak corrected mocks 
                        build_peakcorrected_fibcol(**cat_corr)  # build peak corrected file 

                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4])

                # assign to data columns class
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col]) 
            
            elif DorR == 'random':
                # Random ---------------------------------------------------

                catalog_columns = ['ra', 'dec', 'z', 'nbar']    # catalog columns 

                if (os.path.isfile(file_name) == False) or (clobber == True):
                    print 'Constructing ', file_name 
                    build_random(**cat_corr) 

                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])   # ra, dec, z, nbar

                for i_col, catalog_column in enumerate(catalog_columns):    
                    # assign to object data columns
                    setattr(self, catalog_column, file_data[i_col])

        else: 
            raise NameError('not yet coded') 
        ''' COMMENTED OUT FOR NOW 
        # CMASS --------------------------------------------------------------------------------------------------------
        elif catalog['name'].lower() == 'cmass':
            omega_m = 0.274
            # Only coded for dr12v4 
            # NO CORRECTION IMPOSED YET
            # Data ------------------------------------------------------------------------------------------------------
            if DorR == 'data': 
                catalog_columns = ['ra', 'dec', 'z', 'wsys', 'wnoz', 'wfc', 'nbar', 'comp']
                self.columns = catalog_columns

                data_dir = '/mount/riachuelo1/hahn/data/'              # data directory

                file_name = ''.join([data_dir, 'cmass-dr12v4-N-Reid-weights-zlim.dat'])
                if readdata == True: 
                    file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5,6,7])   #ra,dec,z,wsys,wnoz,wfc,nbar,comp

                    # assign to data columns class
                    for i_col, catalog_column in enumerate(catalog_columns): 
                        setattr(self, catalog_column, file_data[i_col]) 

            elif DorR == 'random': 
                catalog_columns = ['ra', 'dec', 'z', 'nbar', 'comp']
                self.columns = catalog_colums
                
                file_name = ''.join([data_dir, 'cmass-dr12v4-N-Reid-weights-zlim.ran.dat'])
                if readdata == True: 
                    file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4])   #ra,dec,z,nbar,comp

                    # assign to data columns class
                    for i_col, catalog_column in enumerate(catalog_columns): 
                        setattr(self, catalog_column, file_data[i_col]) 
        '''
        
        # survey cosmology metadata 
        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.7 
        cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        self.cosmo = cosmo 

def get_galaxy_data_file(DorR, **cat_corr): 
    ''' Returns file name of Data or Random catalog 
    given cat_corr dictionary 
    
    Parameters
    ----------
    DorR : 'data' or 'random'
    cat_corr : catalog + correction dictionary 

    Returns
    -------
    file_name : D/R file name 
    '''
    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction'] 
    
    if catalog['name'].lower() == 'lasdamasgeo':    # LasDamasGeo -----------------------
        if DorR.lower() == 'data':                  # data

            if correction['name'].lower() == 'true':            # true

                file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/',
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.dat']) 

            elif correction['name'].lower() in ('upweight', 'fibcol', 
                    'shotnoise', 'hectorsn', 'floriansn'): 
                ''' LDG mocks with fiber collision weights 
                '''
                file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/', 
                    'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog['n_mock']), 
                    catalog['letter'], '_no.rdcz.fibcoll.dat']) 

            elif correction['name'].lower() in ('peak', 'peaknbar', 'peakshot', 'peaktest'): 

                if correction['name'].lower()  == 'peak':   # to correct for poor naming convention 
                    correction['name'] = 'peaknbar'
                  
                # specify peak correction fit (expon, gauss, true) 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                elif correction['fit'].lower() in ('true'):
                    corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.fpeak', str(correction['fpeak'])])
                else: 
                    raise NameError('peak fit has to be specified: gauss, expon, true') 
                
                file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/', 
                    'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog['n_mock']), 
                    catalog['letter'], '_no.rdcz.fibcoll.dat'+corr_str]) 
            elif correction['name'].lower() in ('peakshot_dnn'): 
                ''' Peak Correction (with dNN env) + Shotnoise 
                Correction needs to specify: fit, nth NN, and sigma 
                '''

                # specify peak correction fit (expon, gauss) 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    corr_str = ''.join([ '.', correction['fit'].lower(), '.peakshot_d', str(correction['NN']), 'NN', 
                        '.sigma', str(correction['sigma']) ])
                else: 
                    raise NameError('peak fit has to be specified: gauss, expon') 
                
                file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/', 
                    'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog['n_mock']), 
                    catalog['letter'], '_no.rdcz.fibcoll.dat'+corr_str]) 

            elif correction['name'].lower() in ('noweight'):    # noweight

                file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/', 
                    'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog['n_mock']), 
                    catalog['letter'], '_no.rdcz.noweight.dat']) 

            else: 
                raise NotImplementedError('asdfasdf') 

        if DorR.lower() == 'random': 
            file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat'

    elif catalog['name'].lower() == 'tilingmock':   
        # Tiling Mock ---------------------------------
        if DorR == 'data': 
            data_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'      # data directory

            if correction['name'].lower() == 'true': 
                # all weights = 1 (fibercollisions *not* imposed) 
                file_name = data_dir+'cmass-boss5003sector-icoll012.zlim.dat'
                        
            elif correction['name'].lower() in ('upweight', 'fibcol', 'shotnoise', 'floriansn', 'hectorsn'): 
                file_name = data_dir+'cmass-boss5003sector-icoll012.zlim.fibcoll.dat'

            elif correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot'): 

                if correction['name'].lower() == 'peak': 
                    # Correct for poor naming conventions 
                    correction['name'] = 'peaknbar'
                
                # specify fit (expon or gauss) to peak 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                   pass 
                else: 
                    raise NameError('peak fit has to be specified: gauss or expon') 

                # correction string in file name 
                corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                    '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                file_name = data_dir+'cmass-boss5003sector-icoll012.zlim.dat'+corr_str 
                
            # Fibercollisions Corrected by all peak correction -----------------------------------------------------------
            elif correction['name'].lower() in ('allpeak', 'allpeakshot'):

                # specify peak correction fit (expon or gauss) 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    pass 
                else: 
                    raise NameError('peak fit has to be specified: gauss or expon') 

                corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                    '.sigma', str(correction['sigma'])]) 

                file_name = data_dir+'cmass-boss5003sector-icoll012.zlim.dat'+corr_str 
                
            else: 
                raise NameError('Correction Name Unknown') 

        elif DorR.lower() == 'random':              # Randoms --------------------------------------------------------
            file_name = ''.join(['/mount/riachuelo1/hahn/data/tiling_mocks/', 
                'randoms-boss5003-icoll012-vetoed.zlim.dat']) 

    elif catalog['name'].lower() == 'qpm':                  # QPM -------------------- 
        if DorR == 'data':                                  # data  
            data_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'              

            if correction['name'].lower() == 'true': 
                file_name = ''.join([data_dir, 
                    'a0.6452_', str("%04d" % catalog['n_mock']), 
                    '.dr12d_cmass_ngc.vetoed.dat']) 
                
            elif correction['name'].lower() in ('upweight', 'shotnoise', 
                    'floriansn', 'hectorsn'): 
                file_name = ''.join([data_dir, 
                    'a0.6452_', str("%04d" % catalog['n_mock']), 
                    '.dr12d_cmass_ngc.vetoed.fibcoll.dat']) 
                    
            elif correction['name'].lower() in ('peaknbar', 'peakshot'): 
                if correction['name'].lower() == 'peaknbar': 
                    pass
                    # should not happen 
                    #warnings.warn('peaknbar requires corrected nbar(z) and randoms, not coded for QPM') 

                # specify fit (expon or gauss) to peak 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    # correction specifier string 
                    corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])]) 
                elif correction['fit'].lower() in ('true'): 
                    # correction specifier string 
                    corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.fpeak', str(correction['fpeak'])]) 
                else: 
                    raise NameError('peak fit has to be specified: gauss or expon') 

                file_name = ''.join([data_dir, 'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.fibcoll', corr_str, '.dat']) 

            elif correction['name'].lower() in ('peakshot_dnn'):
                
                # specify peak correction fit (expon, gauss) 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    corr_str = ''.join([ '.', correction['fit'].lower(), '.peakshot_d', str(correction['NN']), 'NN', 
                        '.sigma', str(correction['sigma']) ])
                else: 
                    raise NameError('peak fit has to be specified: gauss, expon') 

                file_name = ''.join([data_dir, 
                    'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.fibcoll', corr_str, '.dat']) 

            elif correction['name'].lower() in ('tailupw'):     
                # only upweight uncorrelated chance alignment fc pairs 
                corr_str = 'tailupw' 
                file_name = ''.join([data_dir, 
                    'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.fibcoll.', corr_str, '.dat'
                    ]) 
            elif correction['name'].lower() in ('noweight'):    # no weight

                file_name = ''.join([data_dir, 
                    'a0.6452_', str("%04d" % catalog['n_mock']), 
                    '.dr12d_cmass_ngc.vetoed.noweight.dat']) 

            else: 
                raise NotImplementedError('not yet coded') 
                
        elif DorR == 'random':              # Random 
            data_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/' 
                
            file_name = ''.join([data_dir, 'a0.6452_rand50x.dr12d_cmass_ngc.vetoed.dat'])             # hardcoded to 50x so it does'nt take forever

    elif catalog['name'].lower() == 'patchy': 
        # PATHCY mocks ------------------------------------------------
        data_dir = '/mount/riachuelo1/hahn/data/PATCHY/dr12/v6c/'   # data directory

        if DorR == 'data': 
            # mock catalogs  

            if correction['name'].lower() == 'true': 
                # true mocks
                file_name = ''.join([data_dir, 
                    'Patchy-Mocks-DR12CMASS-N-V6C-Portsmouth-mass_', 
                    str("%04d" % catalog['n_mock']), '.vetoed.dat']) 
                
            elif correction['name'].lower() in ('upweight', 'shotnoise', 
                    'floriansn', 'hectorsn'): 
                # upweighted mocks 
                file_name = ''.join([data_dir, 
                    'Patchy-Mocks-DR12CMASS-N-V6C-Portsmouth-mass_', 
                    str("%04d" % catalog['n_mock']), '.vetoed.fibcoll.dat']) 

            elif correction['name'].lower() in ('peaknbar', 'peakshot'): 
                # peak corrected mocks 
                if correction['name'].lower() == 'peaknbar': 
                    pass

                # specify fit (expon or gauss) to peak 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    # correction specifier string 
                    corr_str = ''.join(['.', correction['fit'].lower(), 
                        '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), 
                        '.fpeak', str(correction['fpeak'])]) 

                elif correction['fit'].lower() in ('true'): 
                    # correction specifier string 
                    corr_str = ''.join(['.', correction['fit'].lower(), 
                        '.', correction['name'].lower(), 
                        '.fpeak', str(correction['fpeak'])]) 
                else: 
                    raise NameError('peak fit has to be specified: gauss or expon') 

                file_name = ''.join([data_dir, 
                    'Patchy-Mocks-DR12CMASS-N-V6C-Portsmouth-mass_', 
                    str("%04d" % catalog['n_mock']), '.vetoed.fibcoll', corr_str, '.dat'
                    ]) 

            elif correction['name'].lower() in ('tailupw'):     
                # only upweight uncorrelated chance alignment fc pairs 
                corr_str = 'tailupw' 
                file_name = ''.join([data_dir, 
                    'Patchy-Mocks-DR12CMASS-N-V6C-Portsmouth-mass_', 
                    str("%04d" % catalog['n_mock']), '.vetoed.fibcoll', corr_str, '.dat']) 
            else: 
                raise NameError('not yet coded') 

        elif DorR == 'random': 
            # random catalog 

            file_name = ''.join([data_dir, 'Random-DR12CMASS-N-V6C-x50.vetoed.dat'])

    return file_name 

# ------------------------------------------------------------------------
# Build galaxy data  

def build_true(**cat_corr): 
    ''' Adjust original data for convenience purposes: 
    
    * Tiling Mock : Reads in true mock catalogs and imposes predetermined redshift limits on them and attaches a flag mostly hardcoded since it's a simple procedure 
    * QPM: Handles messy weights 
    * PATCHY: Returns mocks with only necessary columns and w_fc = 1

    Parameters
    ----------
    cat_corr : catalog correction dictionary 

    '''
    catalog = cat_corr['catalog']
    
    if catalog['name'].lower() == 'tilingmock': 
        # Tiling Mock ------------------------------------------------
        # import original true data 
        orig_true_data = np.loadtxt('/mount/riachuelo1/hahn/data/tiling_mocks/cmass-boss5003sector-icoll012.dat') 
        orig_ra = orig_true_data[:,0]
        orig_dec = orig_true_data[:,1]
        orig_z = orig_true_data[:,2]
        orig_w = orig_true_data[:,3]
    
        zlimit = (orig_z > 0.43) & (orig_z < 0.7)           # tiling mock redshift limit (not universal) 

        true_zlim_file = galaxy_data('data', readdata=False, **cat_corr)
        np.savetxt(true_zlim_file.file_name, np.c_[orig_ra[zlimit], orig_dec[zlimit], orig_z[zlimit], orig_w[zlimit]], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    elif catalog['name'].lower() == 'qpm': 
        # QPM -----------------------------------------------------
        P0 = 20000.             # hardcoded P0 value
        # import original true data 
        orig_true_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
            'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 
        orig_true_data = np.loadtxt(orig_true_file) 

        orig_true_info_file = orig_true_file+'.info'
        orig_true_info = np.loadtxt(orig_true_info_file)    # gal_id, comp, z_real, z_red, mass_halo, flag_sta, id_halo
        
        # consistency issue with #46
        if catalog['n_mock'] in (46, 52, 53, 54, 56, 61, 707, 756, 794, 819, 831, 835, 838): 
            orig_true_veto_file = ''.join(['/mount/riachuelo1/hahn/data/QPM/dr12d/a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
        else:
            orig_true_veto_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
                'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
        orig_true_veto = np.loadtxt(orig_true_veto_file) 
        n_gal = len(orig_true_veto)
       
        # assign RA, Dec, z, and w_veto 
        true_ra = orig_true_data[:,0]
        true_dec = orig_true_data[:,1]
        true_z = orig_true_data[:,2]
        true_wfkp = orig_true_data[:,3]                           # compute nbar(z_i) from w_fkp 
        true_wfc = np.array([1.0 for i in range(n_gal)])          # fiber collisions weights are all 1 for true

        # check to make sure that the redshifts correspond btw rdz file and rdz.info file 
        if np.max(np.abs(orig_true_info[:,3]/true_z-1.0)) > 10**-5: 
            raise NameError('redshifts dont match') 

        true_comp = orig_true_info[:,1]         # completness weights

        # remove veto mask 
        vetomask = (orig_true_veto == 0)            
        # Only keep galaxies with veto = 0 (for veto values in .veto file) 
        
        true_file = galaxy_data('data', readdata=False, **cat_corr)
        np.savetxt(true_file.file_name, np.c_[
            true_ra[vetomask], true_dec[vetomask], true_z[vetomask], 
            true_wfkp[vetomask], true_wfc[vetomask], true_comp[vetomask]], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    
    elif catalog['name'].lower() == 'lasdamasgeo': 
        # Las Damas Geo ------------------------------------------------------
        orig_true_file = ''.join(['/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/', 
            'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.dat']) 
        orig_true_data = np.loadtxt(orig_true_file, unpack=True, usecols=[0,1,2])         # ra, dec, ***CZ***

        true_ra = orig_true_data[0]
        true_dec = orig_true_data[1]
        true_z = orig_true_data[2]/299800.0         # convert cz to z

        true_weight = np.array([1.0 for j in range(len(true_z))])   # no weights for true (all 1) 

        true_file = galaxy_data('data', readdata=False, **cat_corr)

        np.savetxt(true_file.file_name, np.c_[
            true_ra, true_dec, true_z, true_weight],
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    elif catalog['name'].lower() == 'patchy': 
        # PATCHY mocks ------------------------------------------ 
        
        # read original mock data 
        orig_file = ''.join(['/mount/riachuelo1/hahn/data/PATCHY/dr12/v6c/', 
            'Patchy-Mocks-DR12CMASS-N-V6C-Portsmouth-mass_', 
            str("%04d" % catalog['n_mock']), '.dat']) 

        # ra, dec, z, nbar, wfc, veto 
        orig_ra, orig_dec, orig_z, orig_nbar, orig_veto = np.genfromtxt(orig_file, 
                unpack=True, usecols=[0, 1, 2, 4, 6]) 
        n_gal = len(orig_ra) 
        
        new_wfc = np.array([1.0 for i in range(n_gal)])     # w_fc = 1.0 for true 

        vetomask = (orig_veto == 1)            # only keep galaxies with w_veto = 1
        
        true_file = get_galaxy_data_file('data', **cat_corr)    # file name 
        np.savetxt(true_file, 
                np.c_[
                    orig_ra[vetomask], orig_dec[vetomask], orig_z[vetomask], 
                    orig_nbar[vetomask], new_wfc[vetomask]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f'], 
                delimiter='\t') 
    else: 
        raise NameError('not yet coded') 

def build_random(**cat_corr): 
    ''' Build the random catalogs from original data 

    Paramter
    --------
    cat_corr : Catalog and Correction dictionary
    '''
    catalog = cat_corr['catalog']

    if catalog['name'].lower() == 'qpm':    # QPM
        # read original random catalog  
        orig_true_random = np.loadtxt('/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/randoms/a0.6452_rand50x.dr12d_cmass_ngc.rdz')             # ra, dec, z, wfkp
        orig_true_random_info = np.loadtxt('/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/randoms/a0.6452_rand50x.dr12d_cmass_ngc.rdz.info')   # galid, comp?
        orig_true_random_veto = np.loadtxt('/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/randoms/a0.6452_rand50x.dr12d_cmass_ngc.veto')       # veto  

        vetomask = (orig_true_random_veto == 0)
        #true_random = galaxy_data('random', readdata=False, **{'catalog':{'name':'qpm'}, 'correction':{'name':'true'}})
        true_random_file = get_galaxy_data_file('random', **{'catalog':{'name':'qpm'}, 'correction':{'name':'true'}})
        
        np.savetxt(true_random_file, np.c_[(orig_true_random[:,0])[vetomask], (orig_true_random[:,1])[vetomask], 
            (orig_true_random[:,2])[vetomask], (orig_true_random[:,3])[vetomask], (orig_true_random_info[:,1])[vetomask]], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    elif catalog['name'].lower() == 'patchy':       # PATCHY

        orig_file = ''.join(['/mount/riachuelo1/hahn/data/PATCHY/dr12/v6c/', 
            'Random-DR12CMASS-N-V6C-x50.dat']) 
        orig_ra, orig_dec, orig_z, orig_nbar, orig_veto = np.genfromtxt(orig_file, 
                unpack=True, usecols=[0, 1, 2, 3, 5]) 
        
        vetomask = (orig_veto == 1)     # only keep veto = 1
        
        vetoed_file = get_galaxy_data_file('random', **cat_corr) 
        
        np.savetxt(vetoed_file, 
                np.c_[
                    orig_ra[vetomask], orig_dec[vetomask], orig_z[vetomask], orig_nbar[vetomask]
                    ], fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e'], delimiter='\t') 
    else:
        raise NotImplementedError('asdfasdfasdfasdfadf') 

def build_noweight(**cat_corr): 
    ''' Build mock/random data with no fiber collision weights. 
    In other words, galaxies with w_fc = 0 are *not* included in 
    the sample. 

    Parameters
    ----------
    cat_corr : catalog correction dictionary 
    
    Notes 
    -----
    * Tiling Mock : Reads in true mock catalogs and imposes predetermined redshift limits on them and attaches a flag mostly hardcoded since it's a simple procedure 
    * QPM: Handles messy weights 
    * LasDamasGeo : Reads in fibercollided (wcp assigned) mocks then assigns wcp = 1 to all 
    galaxies that have wcp > 0. Remove rest from file
    * PATCHY: Returns mocks with only necessary columns and w_fc = 1

    '''
    catalog = cat_corr['catalog']
    
    if catalog['name'].lower() == 'qpm': 
        # QPM -----------------------------------------------------
        
        # import original true data 
        orig_true_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
            'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 
        orig_true_data = np.loadtxt(orig_true_file) 

        orig_true_info_file = orig_true_file+'.info'
        orig_true_info = np.loadtxt(orig_true_info_file)    # gal_id, comp, z_real, z_red, mass_halo, flag_sta, id_halo
        
        # consistency issue with #46
        if catalog['n_mock'] in (46, 52, 53, 54, 56, 61, 707, 756, 794, 819, 831, 835, 838): 
            orig_true_veto_file = ''.join(['/mount/riachuelo1/hahn/data/QPM/dr12d/a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
        else:
            orig_true_veto_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
                'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
        orig_true_veto = np.loadtxt(orig_true_veto_file) 
        n_gal = len(orig_true_veto)
       
        # assign RA, Dec, z, and w_veto 
        orig_ra = orig_true_data[:,0]
        orig_dec = orig_true_data[:,1]
        orig_z = orig_true_data[:,2]
        orig_wfkp = orig_true_data[:,3]         # compute nbar(z_i) from w_fkp 
        orig_wfc = orig_true_data[:,4]          # fiber collisions weights are all 1 for true

        # check to make sure that the redshifts correspond btw rdz file and rdz.info file 
        if np.max(np.abs(orig_true_info[:,3]/orig_z-1.0)) > 10**-5: 
            raise NameError('redshifts dont match') 

        orig_comp = orig_true_info[:,1]         # completness weights

        # remove veto mask 
        vetomask = (orig_true_veto == 0) & (orig_wfc >= 1) 
        # Only keep galaxies with veto = 0 (for veto values in .veto file) and
        # wfc >= 1 (remove all fiber collided pairs that aren't included in the actual catlaog 

        vetoed_ra = orig_ra[vetomask]
        vetoed_dec = orig_dec[vetomask]
        vetoed_z = orig_z[vetomask]
        vetoed_wfkp = orig_wfkp[vetomask]
        vetoed_comp = orig_comp[vetomask]
        n_veto = len(vetoed_ra) 
        vetoed_wfc = np.array([1.0 for i in range(n_veto)]) 

        noweight_file = get_galaxy_data_file('data', **cat_corr) 
        np.savetxt(noweight_file, np.c_[
            vetoed_ra, vetoed_dec, vetoed_z, vetoed_wfkp, vetoed_wfc, vetoed_comp], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    
    elif catalog['name'].lower() == 'lasdamasgeo': 
        # Las Damas Geo ------------------------------------------------------
        ldg_cat_corr = {'catalog': catalog, 'correction': {'name': 'upweight'}}
        fibcol_file = get_galaxy_data_file('data', **ldg_cat_corr) 
        fc_ra, fc_dec, fc_z, fc_w = np.loadtxt(fibcol_file, unpack=True, usecols=[0,1,2,3])
        
        hasz = fc_w > 0.0 
        n_hasz = len(fc_ra[hasz]) 
        
        weights = np.array([1.0 for i in range(n_hasz)]) 
    
        now_file = get_galaxy_data_file('data', **cat_corr)
        np.savetxt(now_file, np.c_[fc_ra[hasz], fc_dec[hasz], fc_z[hasz], weights],
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    else: 
        raise NotImplementedError('not yet coded') 

def build_fibercollided(**cat_corr): 
    ''' Build Fibercollided mock catalogs 
    using specific idl routines

    Parameters
    ----------
    cat_corr : catalog correction dictionary 

    Notes
    -----
    '''
    catalog = cat_corr['catalog']
    
    if catalog['name'].lower() == 'lasdamasgeo':        # Las Damas Geo 
        
        fibcollided_cmd = 'idl -e "ldg_fibcollmock_wcp_assign,'+str(catalog['n_mock'])+", '"+\
                str(catalog['letter'])+"'"+'"'
        print fibcollided_cmd
        os.system(fibcollided_cmd)  # call IDL code 

    elif catalog['name'].lower() == 'tilingmock':       # Tiling Mock 
        fibcollided_cmd = ' '.join(['idl', '-e', '"', "build_wcp_assign, 'tilingmock'", '"'])
        os.system(fibcollided_cmd) 

    elif catalog['name'].lower() == 'qpm': 
        # QPM ------------------------------------------------------------
        orig_true_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
            'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 
        orig_true_data = np.loadtxt(orig_true_file) 
        
        orig_true_info_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
            'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz.info']) 
        # gal_id, comp, z_real, z_red, mass_halo, flag_sta, id_halo
        orig_true_info_data = np.loadtxt(orig_true_info_file)    

        if catalog['n_mock'] in (44, 46, 52, 53, 54, 56, 61, 707, 756, 794, 819, 831, 835, 838): 
            orig_true_veto_file = ''.join(['/mount/riachuelo1/hahn/data/QPM/dr12d/a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
        else:
            orig_true_veto_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
                'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
        orig_true_veto = np.loadtxt(orig_true_veto_file) 
        n_gal = len(orig_true_veto)

        if len(orig_true_data[:,0]) != n_gal: 
            print orig_true_file
            print orig_true_veto_file 
            raise ValueError('veto mask doesnt match') 
       
        # assign RA, Dec, z, and w_veto 
        true_ra = orig_true_data[:,0]
        true_dec = orig_true_data[:,1]
        true_z = orig_true_data[:,2]
        true_wfkp = orig_true_data[:,3]
        true_wfc = orig_true_data[:,4] 

        true_comp = orig_true_info_data[:,1]

        vetomask = (orig_true_veto == 0)

        fc_file = galaxy_data('data', readdata=False, **cat_corr)
        np.savetxt(fc_file.file_name, np.c_[
            true_ra[vetomask], true_dec[vetomask], true_z[vetomask], 
            true_wfkp[vetomask], true_wfc[vetomask], true_comp[vetomask]], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

        fibcollided_cmd = ''

    elif catalog['name'].lower() == 'patchy': 
        # PATCHY mocks ----------------------------------------

        # read original mock data 
        orig_file = ''.join(['/mount/riachuelo1/hahn/data/PATCHY/dr12/v6c/', 
            'Patchy-Mocks-DR12CMASS-N-V6C-Portsmouth-mass_', 
            str("%04d" % catalog['n_mock']), '.dat']) 

        # ra, dec, z, nbar, wfc, veto 
        orig_ra, orig_dec, orig_z, orig_nbar, orig_wfc, orig_veto = np.genfromtxt(orig_file, 
                unpack=True, usecols=[0, 1, 2, 4, 7, 6]) 
        n_gal = len(orig_ra) 

        vetomask = (orig_veto == 1)            # only keep galaxies with w_veto = 1
        
        fc_file = get_galaxy_data_file('data', **cat_corr)    # file name 
        np.savetxt(fc_file, 
                np.c_[
                    orig_ra[vetomask], orig_dec[vetomask], orig_z[vetomask], 
                    orig_nbar[vetomask], orig_wfc[vetomask]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f'], 
                delimiter='\t') 
        
        fibcollided_cmd = ''
    else: 
        raise NameError('not yet coded') 

    return fibcollided_cmd 

def build_peakcorrected_fibcol(sanitycheck=False, **cat_corr): 
    ''' Build peak corrected fibercollided mock catalogs (using cosmolopy) 

    Parameters
    ----------
    cat_corr : Catalog + Correction dictionary
    sanitycheck : testing flag 

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # fit functions (using lambda) 
    if correction['fit'].lower() == 'gauss': 
        fit_func = lambda x, sig: np.exp(-0.5*x**2/sig**2)
    elif correction['fit'].lower() == 'expon': 
        fit_func = lambda x, sig: np.exp(-1.0*np.abs(x)/sig)
    elif correction['fit'].lower() == 'true': 
        pass 
    else: 
        raise NameError('correction fit has to be specified as gauss or expon') 

    if catalog['name'].lower() == 'lasdamasgeo': 
        # Las Damas Geo -----------------------------------------------------------

        n_mocks = 160   # total number of mocks

        survey_zmin, survey_zmax = 0.16, 0.44
        
        # only one catalog parameter is defined
        if (isinstance(catalog['n_mock'], int) == False) or \
                (isinstance(catalog['letter'], str) == False): 
            raise NameError('only one mock can be corrected at a time') 
        
        # read in fiber collided mock  
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr) 
        cosmo = fibcoll_mock.cosmo           # survey comoslogy 

        survey_comdis_min = \
                cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
        survey_comdis_max = \
                cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']

        # set peak fraction         
        if correction['name'].lower() in ('allpeak', 'allpeakshot'): 
            f_peak = 1.0
        else: 
            f_peak = correction['fpeak'] 

        appended_ra, appended_dec, appended_z, appended_weight = [], [], [], []
        upweight_again = []
            
        if sanitycheck == True:     # check that the peak p(r) is generated properly
            pr_test = [] 
    
        for i_mock in range(len(fibcoll_mock.weight)): 
            # go through every galaxy in fibercollided mock catalog

            while fibcoll_mock.weight[i_mock] > 1: 
                # for galaxies with wcp > 1
                fibcoll_mock.weight[i_mock] = fibcoll_mock.weight[i_mock]-1.0

                # LOS comoving distance of the galaxy 
                comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], 
                        **cosmo)*cosmo['h']
                
                rand_num = np.random.random(1) 
                
                if rand_num <= f_peak:          # if in the peak 
                    # keep ra and dec
                    appended_ra.append(fibcoll_mock.ra[i_mock])
                    appended_dec.append(fibcoll_mock.dec[i_mock])

                    if correction['name'].lower() in ('allpeak', 'allpeakshot'): 
                        # appended galaxy has weight of 1.0 
                        appended_weight.append(correction['fpeak'])
                    else: 
                        # appended galaxy has weight of 1.0 
                        appended_weight.append(1.0)
                    if correction['fit'].lower() in ('gauss', 'expon'):   
                        # compute the displacement within peak using best-fit --------------
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*correction['sigma']
                            peakpofr = fit_func(rand2, correction['sigma']) 
                        #--------------------------------------------------------------------- 

                    elif correction['fit'].lower() == 'true': 
                        # compute the displacement within peak using actual distribution   
                        dlos_comb_peak_file = ((fibcoll_mock.file_name).rsplit('/', 1))[0]+'/DLOS_norm_peak_dist_'+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined.dat'
                        dlos_mid, dlos_dist = np.loadtxt(dlos_comb_peak_file, unpack=True, usecols=[0,1])

                        dlos_cdf = dlos_dist.cumsum() 

                        rand1 = np.random.random(1) 
                        
                        cdf_closest_index = min(range(len(dlos_cdf)), key = lambda i: abs(dlos_cdf[i]-rand1[0])) 
                        closest_dlos = dlos_mid[cdf_closest_index] 
                       
                        try: 
                            closest_dloses
                        except NameError:
                            closest_dloses = [closest_dlos]
                        else: 
                            closest_dloses.append(closest_dlos)

                        rand2 = np.array([closest_dlos])
                        #--------------------------------------------------------------------- 

                    # in case the displacement falls out of bound (may general large scale issues)
                    if (comdis_imock+rand2 > survey_comdis_max) or \
                            (comdis_imock+rand2 < survey_comdis_min): 
                        collided_z = comdis2z(comdis_imock-rand2, **cosmo)
                    else: 
                        collided_z = comdis2z(comdis_imock+rand2, **cosmo)

                    appended_z.append(collided_z[0]) 

                    if correction['name'] == 'allpeakshot': 
                        # need to re-upweight the fibercollided galaxies by 1-fpeak 
                        upweight_again.append(i_mock)       # save index to be re-upweighted 

                    if sanitycheck == True: 
                        pr_test.append(rand2) 
                else:                           # if not in the peak 
                    if correction['name'] == 'peaktest': 
                        # do nothing 
                        pass
                    elif correction['name'] in ('peak', 'peaknbar'): 
                        # randomly displace 
                        appended_ra.append(fibcoll_mock.ra[i_mock])
                        appended_dec.append(fibcoll_mock.dec[i_mock])
                        appended_weight.append(1.0) 
                        
                        rand1 = np.random.random(1) 
                        appended_z.append(0.16+rand1[0]*0.28)

                        if sanitycheck == True: 
                            original_dm = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], **cosmo)*cosmo['h']
                            corrected_dm = cosmos.distance.comoving_distance(0.16+rand1[0]*0.28, **cosmo)*cosmo['h']
                            pr_test.append(corrected_dm-original_dm)

                    elif (correction['name'] == 'allpeak') or (correction['name'] == 'allpeakshot'): 
                        # it should all be in the peak!
                        raise NameError('should not happen')

                    elif correction['name'] == 'peakshot': 
                        upweight_again.append(i_mock)
        
        if correction['name'] == 'peakshot': 
            # re-upweighting for peak+shotnoise correction 
            for i_upweightagain in upweight_again: 
                fibcoll_mock.weight[i_upweightagain] = fibcoll_mock.weight[i_upweightagain]+1.

        elif correction['name'] == 'allpeakshot': 
            # re-upweighting for all peak+shotnoise correction 
            for i_upweightagain in upweight_again: 
                fibcoll_mock.weight[i_upweightagain] = fibcoll_mock.weight[i_upweightagain]+(1.0-correction['fpeak']) 

        print len(appended_ra), ' galaxies were peak corrected'
        fibcoll_mock.ra = np.concatenate([fibcoll_mock.ra, appended_ra])
        fibcoll_mock.dec = np.concatenate([fibcoll_mock.dec, appended_dec])
        fibcoll_mock.weight = np.concatenate([fibcoll_mock.weight, appended_weight])
        fibcoll_mock.z = np.concatenate([fibcoll_mock.z, appended_z])

    elif catalog['name'].lower() in ('tilingmock', 'qpm', 'patchy'): 
        # Tiling mock, QPM, PATCHY 
        
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits

        if catalog['name'].lower() == 'qpm':
            n_mocks = 100
        
        # read in mock with fibercollisions imposed
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr)  
        
        cosmo = fibcoll_mock.cosmo      # survey cosmology 

        survey_comdis_min = cosmos.distance.comoving_distance(survey_zmin, 
                **cosmo)*cosmo['h']         # in units of Mpc/h
        survey_comdis_max = cosmos.distance.comoving_distance(survey_zmax, 
                **cosmo)*cosmo['h']         # in units of Mpc/h

        if catalog['name'].lower() in ('qpm', 'patchy'): 
            # only use fiber collision weights
            fibcoll_mock.weight = fibcoll_mock.wfc            
                        
        # read in the true galaxies for tail portion  
        # the tail portion of the peak corrections will be generated similar to the mksample procedures for 
        # assigning the random catalog redshift 
        true_cat_corr = {'catalog':catalog, 'correction': {'name': 'true'}}
        true_data = galaxy_data('data', **true_cat_corr)
        
        if catalog['name'].lower() in ('qpm',  'patchy'): 
            true_weight = true_data.wfc
        else: 
            true_weight = true_data.weight
        true_z = true_data.z

        true_weight_cum = true_weight.cumsum()
        true_weight_max = true_weight_cum.max()
        i_true = np.arange(0,len(true_weight_cum))

        # set peak fraction (f_peak in correction direction reflects the correct fpeak, 
        # but for all peak correction methods this is simply used as weight 
        if correction['name'].lower() in ('allpeakshot', 'allpeak'): 
            f_peak = 1.0
        else: 
            f_peak = correction['fpeak'] 
    
        append_ra, append_dec, append_z, append_dlos, append_weight = [], [], [], [], []

        if catalog['name'].lower() == 'qpm': 
            append_wfkp, append_comp = [], [] 
        elif catalog['name'].lower() == 'patchy': 
            append_nbar = [] 

        reupw = []

        if sanitycheck == True:     # check peak p(r) 
            pr_test = [] 
    
        for i_mock in range(len(fibcoll_mock.weight)):  # go through each galaxy 

            while fibcoll_mock.weight[i_mock] > 1:      # for galaxies with wcp > 1

                fibcoll_mock.weight[i_mock] -= 1.0

                # LOS comoving distance of the galaxy 
                comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], 
                        **cosmo)*cosmo['h']
                
                rand_num = np.random.random(1)  # determine whether galaxy is in the peak
               
                if rand_num < f_peak:          # in the peak 

                    append_ra.append(fibcoll_mock.ra[i_mock])    # keep ra and dec
                    append_dec.append(fibcoll_mock.dec[i_mock])

                    if catalog['name'].lower() == 'qpm': 
                        append_wfkp.append(fibcoll_mock.wfkp[i_mock]) 
                        append_comp.append(fibcoll_mock.comp[i_mock])
                    elif catalog['name'].lower() == 'patchy': 
                        append_nbar.append(fibcoll_mock.nbar[i_mock]) 

                    append_weight.append(1.0)       # appended galaxy has weight=1.0 
                    
                    if correction['fit'].lower() in ('gauss', 'expon'): 
                        # compute the displacement within peak ------------------                        
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*correction['sigma']
                            peakpofr = fit_func(rand2, correction['sigma']) 
                        #----------------------------------- --------------------------------- 

                    elif correction['fit'].lower() == 'true': 
                        # compute displacements within peak using true distribution ------
                        dlos_comb_peak_file = ((fibcoll_mock.file_name).rsplit('/', 1))[0]+'/DLOS_norm_peak_dist_'+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined.dat'
                        dlos_mid, dlos_dist = np.loadtxt(dlos_comb_peak_file, unpack=True, usecols=[0,1])

                        dlos_cdf = dlos_dist.cumsum() 

                        rand1 = np.sum(dlos_dist)*np.random.random(1) 
                        
                        cdf_closest_index = min(range(len(dlos_cdf)), key = lambda i: abs(dlos_cdf[i]-rand1[0])) 
                        if rand1 > dlos_cdf[cdf_closest_index]: 
                            cdf_smaller_index = cdf_closest_index
                            cdf_larger_index = cdf_closest_index+1
                        else: 
                            cdf_smaller_index = cdf_closest_index-1
                            cdf_larger_index = cdf_closest_index

                        # linear interpolations
                        closest_dlos = dlos_mid[cdf_smaller_index] + \
                                (dlos_mid[cdf_larger_index] - dlos_mid[cdf_smaller_index])*(rand1[0]-dlos_cdf[cdf_smaller_index])/(dlos_cdf[cdf_larger_index]-dlos_cdf[cdf_smaller_index])
                       
                        try: 
                            closest_dloses
                        except NameError:
                            closest_dloses = [closest_dlos]
                        else: 
                            closest_dloses.append(closest_dlos)

                        rand2 = np.array([closest_dlos])
                        #--------------------------------------------------------------------- 
                    else: 
                        raise NotImplementedError('does not work') 

                    # in case the displacement falls out of bound 
                    # NOTE: THIS IS NOT CORRECT, BUT IMPLEMENTED FOR SIMPLICITY 
                    if (comdis_imock+rand2 > survey_comdis_max) or \
                            (comdis_imock+rand2 < survey_comdis_min): 
                        collided_z = comdis2z(comdis_imock - rand2, **cosmo)
                    else: 
                        collided_z = comdis2z(comdis_imock + rand2, **cosmo)

                    append_z.append(collided_z[0]) 
                    append_dlos.append(rand2)
                    
                    if correction['name'] == 'allpeakshot': 
                        # for allpeak + shot noise correction 
                        # need to re-upweight the fibercollided galaxies by 1-fpeak 
                        reupw.append(i_mock)       # save index to be re-upweighted 

                    if sanitycheck == True:             
                        # append calculated LOS displacement 
                        pr_test.append(rand2) 

                else:                                           # Not in peak 
                    if correction['name'] == 'peaktest': 
                        # do nothing. 
                        # This effectively discards the galaxies that fall into the tail 
                        pass

                    elif (correction['name'] == 'peak') or (correction['name'] == 'peaknbar'): 
                        '''
                        generate random z based on redshift distribution of galaxies 
                        '''
                        # RA, Dec remain the same 
                        # weight = 1
                        append_ra.append(fibcoll_mock.ra[i_mock])
                        append_dec.append(fibcoll_mock.dec[i_mock])
                        append_weight.append(1.0) 
                        
                        wtarg = np.random.random(1)*true_weight_max
                        zindx = np.floor(np.interp(wtarg, true_weight_cum, i_true)).astype(int)+1
                        
                        # fail safe
                        qqq = np.where(wtarg < true_weight_cum[0])[0]
                        zindx[qqq] = 0 

                        # assign redshift 
                        append_z.append(true_z[zindx]) 

                    elif correction['name'] == 'peakshot': 
                        # peak + shot noise correction 
                        # record galaxy index that need to be upweighted again 
                        reupw.append(i_mock)
                    
                    elif correction['name'] in ('allpeakshot', 'allpeak'): 
                        # for all peak correction methods, all galaxies should be in the peak!
                        raise NameError('should not happen')

        if correction['name'] == 'peakshot': 
            # re-upweighting for peak+shotnoise correction 
            for i_reupw in reupw: 
                fibcoll_mock.weight[i_reupw] += 1.0

        elif correction['name'] == 'allpeakshot': 
            # re-upweighting for all peak+shotnoise correction 
            for i_reupw in reupw: 
                fibcoll_mock.weight[i_reupw] += 1.0-correction['fpeak']

        print len(append_ra), ' galaxies were peak corrected'
    
        # append the "append" galaxies to the end of the arrays
        fibcoll_mock.ra = np.concatenate([fibcoll_mock.ra, append_ra])
        fibcoll_mock.dec = np.concatenate([fibcoll_mock.dec, append_dec])
        fibcoll_mock.weight = np.concatenate([fibcoll_mock.weight, append_weight])

        if catalog['name'].lower() == 'qpm': 
            fibcoll_mock.wfkp = np.concatenate([fibcoll_mock.wfkp, append_wfkp])
            fibcoll_mock.comp = np.concatenate([fibcoll_mock.comp, append_comp])
        elif catalog['name'].lower() == 'patchy': 
            fibcoll_mock.nbar = np.concatenate([fibcoll_mock.nbar, append_nbar])

        fibcoll_mock.z = np.concatenate([fibcoll_mock.z, appended_z])
        
        if sanitycheck == True:         # output sanitycheck 
            sanitycheck_file = ''.join([
                'peak_corrected_dlos_values_', correction['name'].lower(), '_', 
                correction['fit'].lower(), '.dat']) 
            np.savetxt(sanitycheck_file,
                    np.c_[append_dlos], fmt=['%10.5f']) 
    else: 
        raise NameError('Error here') 

    peakcorr_file = get_galaxy_data_file('data', **cat_corr) 
    
    if catalog['name'].lower() == 'qpm': 
        np.savetxt(peakcorr_file, 
                np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, 
                    fibcoll_mock.wfkp, fibcoll_mock.weight, fibcoll_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 

    elif catalog['name'].lower() == 'patchy': 
        np.savetxt(peakcorr_file, 
                np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, 
                    fibcoll_mock.nbar, fibcoll_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 

    elif catalog['name'].lower() == 'tilingmock': 
        np.savetxt(peakcorr_file, 
                np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    if sanitycheck == True: 
        np.savetxt(peakcorr_file+'.sanitycheck', np.c_[pr_test], fmt=['%10.5f'], delimiter='\t') 

def build_corrected_randoms(sanitycheck=False, **cat_corr): 
    ''' 
    Construct corrected randoms based on corrected data
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # Las Damas Geo ------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        # read series of corrected data in order to minimize single catalog dependent variance?
        # especially for las damas, which is divided into 4 parts 
        print 'Combining Mocks'  
        for i_mock in range(1,41): 
            for letter in ['a', 'b', 'c', 'd']: 
                i_cat = {'name':catalog['name'], 'n_mock': i_mock, 'letter': letter} 
                i_cat_corr = {'catalog': i_cat, 'correction': correction} 
                corr_data = galaxy_data('data', **i_cat_corr)
                
                print corr_data.file_name
            
                # target weights based on the corrected weights
                if (i_mock == 1) and (letter == 'a'): 
                    targ_weight = corr_data.weight 
                    targ_z = corr_data.z
                else: 
                    targ_weight = np.concatenate((targ_weight, corr_data.weight))
                    targ_z = np.concatenate((targ_z, corr_data.z))
    # Tiling Mocks ----------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'tilingmock': 
        # build corrected randoms for tiling mocsk 
        # read corrected data (Only one mock available) 
        corr_data = galaxy_data('data', **cat_corr)
        print 'Generating random catalog based on : ', corr_data.file_name
            
        # target weights based on the corrected weights
        targ_weight = corr_data.weight 
        targ_z = corr_data.z

    # QPM Mocks --------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'qpm': 
        # read series of corrected data in order to minimize single catalog dependent variance?
        # especially for las damas, which is divided into 4 parts 
        print 'Combining Mocks'  
        for i_mock in range(1,41):                                     #### CURRENTLY ONLY COMBINING 100 BECAUSE 1000 IS TOO FREAKING MANY!
            i_cat = {'name':catalog['name'], 'n_mock': i_mock} 
            i_cat_corr = {'catalog': i_cat, 'correction': correction} 
            corr_data = galaxy_data('data', **i_cat_corr)
                
            print 'combining ... ', corr_data.file_name
            
            # target weights based on the corrected weights
            if i_mock == 1: 
                targ_weight = corr_data.wfc                         # target weight IGNORES veto weights which the randoms take into account separately 
                targ_z = corr_data.z
            else: 
                targ_weight = np.concatenate((targ_weight, corr_data.wfc))
                targ_z = np.concatenate((targ_z, corr_data.z))
    else: 
        # not yet coded 
        raise NameError('not coded yet') 

    targ_weight_cum = targ_weight.cumsum()            # cumulative weights
    targ_weight_max = targ_weight_cum.max()           # max weight (a.k.a. total weight) 
    i_targ = np.arange(0,len(targ_weight_cum))
    print 'Corrected Galaxy Weight Sum = ', targ_weight_max
    
    if catalog['name'].lower() == 'tilingmock': 
        # for tiling mock import original randoms 
        rand_true = {'catalog': catalog, 'correction': {'name': 'true'}}
        true_random = galaxy_data('random', readdata=False, **rand_true)
        
        true_random_data = np.loadtxt('/mount/riachuelo1/hahn/data/tiling_mocks/randoms-boss5003-icoll012-vetoed.dat', unpack=True, usecols=[0,1]) 
        true_random.ra = true_random_data[0]
        true_random.dec = true_random_data[1]
    else:
        # read true random
        rand_true = {'catalog': catalog, 'correction': {'name': 'true'}}
        true_random = galaxy_data('random', **rand_true)

    Nran = len(true_random.ra)  # number of random galaxies
    print 'N_random = ', Nran

    corr_random_z = np.zeros(Nran) 
    # divide into chunks like make_catalog_z.py
    nchunk = 50 
    for i_chunk in range(nchunk): 
        # chunk indicies
        istart = Nran/(nchunk)*i_chunk
        iend = Nran/(nchunk)*(i_chunk+1)
        if i_chunk == nchunk-1: 
            iend = Nran 
        
        # randomly sample the sum of the target weights 
        wtarg = np.random.random(iend-istart)*targ_weight_max

        # get zindex of the target 
        zindx = np.floor(np.interp(wtarg, targ_weight_cum, i_targ)).astype(int)+1       # identical to make_catalog_z.py
        qqq = np.where(wtarg < targ_weight_cum[0])[0]
        zindx[qqq]  = 0 

        # assign redshift 
        corr_random_z[istart:iend] = targ_z[zindx]

    if len(corr_random_z) != Nran: 
        raise NameError("NOT ENOUGH CORRECTED RANDOMS!") 

    ## write Corrected Random File  
    corr_random = galaxy_data('random', readdata=False, **cat_corr)    
    corr_random_filename = corr_random.file_name 

    # write corrected randoms
    if catalog['name'].lower() == 'qpm': 
        # QPM has different data format so needs to be separate
        np.savetxt(corr_random_filename, np.c_[true_random.ra, true_random.dec, corr_random_z, true_random.nbar, true_random.wveto], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%1.5e', '%10.5f'], delimiter='\t') 
    else: 
        np.savetxt(corr_random_filename, np.c_[true_random.ra, true_random.dec, corr_random_z], 
                fmt=['%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def build_fc_info_catalog(**cat_corr): 
    '''
    build mock catalogs with fiber-collision information (upweighted neighbor info) using specific idl routines
    '''
    catalog = cat_corr['catalog']
    
    # QPM --------------------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'qpm': 
        upw_cat_corr = {'catalog': catalog, 'correction': {'name': 'upweight'}}     # upweight cat_corr dict
        fc_upw = galaxy_data('data', readdata=False, **upw_cat_corr)                # get upweight name 
    
        fc_upw_info_file_name = '.fcinfo.'.join((fc_upw.file_name).rsplit('.', 1)) 
        
        fc_info_cat_cmd = ''.join(['idl -e "', 
            "build_fibcoll_catalog, 'qpm', ", 
            "'"+fc_upw.file_name+"'", ", '"+fc_upw_info_file_name+"'"+'"'])
        
        print fc_info_cat_cmd
        os.system(fc_info_cat_cmd) 
    else: 
        raise NameError('not yet coded') 

    return fc_info_cat_cmd

def build_tailupweight_fibcol(**cat_corr): 
    ''' Only upweight galaxies in the tail but keep everything else true 
    tail definition is where abs(dLOS) > 20 Mpc/h
    '''
    catalog = cat_corr['catalog']

    if catalog['name'].lower() == 'qpm': 
        upw_cat_corr = {'catalog': catalog, 'correction': {'name': 'upweight'}}     # upweight cat_corr dict
        fc_upw = galaxy_data('data', readdata=False, **upw_cat_corr)                # get upweight name 
    
        fc_upw_info_file_name = '.fcinfo.'.join((fc_upw.file_name).rsplit('.', 1)) 
    
        if os.path.isfile(fc_upw_info_file_name) == False:                          # if info file doesn't exist
            build_fc_info_catalog(**cat_corr) 

        # read ra, dec, z, wfkp, wcp, wcomp, d_fc from fiber-collision info catalog
        mock_ra, mock_dec, mock_z, wfkp, wcp, wcomp, dfc, fc_index = np.loadtxt(
                fc_upw_info_file_name, unpack=True, usecols=[0,1,2,3,4,5,6,10] 
                ) 

        # run through each line
        for i_row in range(len(mock_ra)): 
            if dfc[i_row] != 0.0:       # if there is some sort of fiber-collision
                if np.abs(dfc[i_row]) > 20.0:    # if in the tail then keep upweights
                    pass
                else:                            # undo upweighting scheme 
                    if wcp[i_row] != 0.0:       
                        # only down-weighted fc-pairs should have dLOS_fc != 0.0 so catch errors
                        raise NameError('something went wrong') 

                    wcp[i_row] = 1.0 
                    wcp[fc_index[i_row]] = wcp[fc_index[i_row]]-1.0         # subtract from fc weight 
    

        tailupw = galaxy_data('data', readdata=False, **cat_corr)
        np.savetxt(tailupw.file_name, 
                np.c_[
                    mock_ra, mock_dec, mock_z, wfkp, wcp, wcomp
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    else: 
        raise NameError('not yet coded')

def build_peak_fpeak_dNN(NN=3, **cat_corr): 
    ''' Build corrected mock catalogs using peak correction with *env dependent fpeak* 

    - Calculates kth nearest neighbor distances for upweighted galaxies 
    - Use bestfit fpeak(d_NN) to calculate the fpeak
    - Peak correct using fpeak 

    Input: 
    **cat_corr dictionary 
    
    (using cosmolopy) 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # lambda function 
    if correction['fit'].lower() == 'gauss': 
        fit_func = lambda x, sig: np.exp(-0.5*x**2/sig**2)
    elif correction['fit'].lower() == 'expon': 
        fit_func = lambda x, sig: np.exp(-1.0*np.abs(x)/sig)
    elif correction['fit'].lower() == 'true': 
        pass 
    else: 
        raise NameError('correction fit has to be specified as gauss or expon') 

    # LDG -------------------------------------------------------------------- 
    if catalog['name'].lower() == 'lasdamasgeo': 

        survey_zmin, survey_zmax = 0.16, 0.44      # survey redshift limit
        n_mocks = 160

        # LOS comoving distance for redshift limit
        #comdis_lo = cosmos.distance.comoving_distance(0.16, **cosmo)*cosmo['h']
        #comdis_hi = cosmos.distance.comoving_distance(0.44, **cosmo)*cosmo['h'] 
        
        # in case more than one catalog parameter is defined
        if (isinstance(catalog['n_mock'], int) == False) or \
                (isinstance(catalog['letter'], str) == False): 
            raise NameError('only one mock can be corrected at a time') 
        
        # read in fiber collided galaxies
        fc_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fc_mock = galaxy_data('data', **fibcoll_cat_corr) 
        n_gal = len(fc_mock.weight)        # Ngal 

        cosmo = fc_mock.cosmo                   # import cosmology from galaxy class 
        # survey comoving distnace limits 
        survey_comdis_min = \
                cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
        survey_comdis_max = \
                cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']

        upw_index = np.arange(n_gal)[ fc_mock.weight > 1 ]  # index of upweighted galaxy 

        upw_dNN = genv.d_NN(fc_mock.ra[upw_index], fc_mock.dec[upw_index], fc_mock.z[upw_index], 
                n=NN, **fc_cat_corr)       # upweight kth Nearest Neighbor distance
        
        fake_ra, fake_dec, fake_z, fake_weight, re_upw = [], [], [], [], []
        for ii_upw, i_upw in enumerate(upw_index):     # go through upweighted galaxy catalog 

            while fc_mock.weight[i_upw] > 1: # for galaxies with wcp > 1

                fc_mock.weight[i_upw] = fc_mock.weight[i_upw] - 1.0

                # LOS comoving distance of the galaxy 
                comdis_iupw = cosmos.distance.comoving_distance(\
                        fc_mock.z[i_upw], **cosmo) * cosmo['h']

                fpeak_iupw = fc_dlos.fpeak_dNN( upw_dNN[ii_upw], n=NN, **fc_cat_corr ) 
                
                rand_num = np.random.random(1)      # random # sampled for the peak 
                if rand_num <= fpeak_iupw:          # IN THE PEAK  

                    # place fake galaxy a dLOS apart from the upweighted galaxy
                    fake_ra.append(fc_mock.ra[i_upw])
                    fake_dec.append(fc_mock.dec[i_upw])
                    fake_weight.append(1.0)

                    if correction['fit'].lower() in ('gauss', 'expon'):   
                        # compute the displacement within peak using best-fit 
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*correction['sigma']
                            peakpofr = fit_func(rand2, correction['sigma']) 

                    # if fake galaxy out of bound (may generate large scale issues)
                    # put it in the opposite sign 
                    # **THINK ABOUT THIS OMRE *** 
                    if (comdis_iupw + rand2 > survey_comdis_min) or \
                            (comdis_iupw + rand2 < survey_comdis_max): 
                        collided_z = comdis2z(comdis_iupw - rand2, **cosmo)
                    else: 
                        collided_z = comdis2z(comdis_imock + rand2, **cosmo)

                    fake_z.append(collided_z[0])        # append collided z

                else:                       # NOT IN THE PEAK 
                    if correction['name'] == 'peakshot': 
                        re_upw.append(i_upw)
                    else: 
                        raise NotImplementedError
        
        if correction['name'] == 'peakshot': 
            # re-upweighting for peak+shotnoise correction 
            for i_reupw in re_upw: 
                fc_mock.weight[i_reupw] += 1.

        print len(fake_ra), ' galaxies were peak corrected'
        # add the fake galaxies to the catalog 
        fc_mock.ra = np.concatenate([fc_mock.ra, fake_ra])
        fc_mock.dec = np.concatenate([fc_mock.dec, fake_dec])
        fc_mock.z = np.concatenate([fc_mock.z, fake_z])
        fc_mock.weight = np.concatenate([fc_mock.weight, fake_weight])

    # QPM / Tiling Mock ------------------------------------------------
    elif catalog['name'].lower() in ('tilingmock', 'qpm'): 

        survey_zmin, survey_zmax = 0.43, 0.7        # survey redshift limit 
        
        if catalog['name'].lower() == 'qpm':
            n_mocks = 100
        
        # read in fibercollided mocks
        fc_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fc_mock = galaxy_data('data', **fc_cat_corr) 

        cosmo = fc_mock.cosmo                       # set comoslogy 
        # redshift limits comoving distance  
        survey_comdis_min = \
                cosmos.distance.comoving_distance(survey_zmin, **cosmo)*cosmo['h']  # in units of Mpc/h
        survey_comdis_max = \
                cosmos.distance.comoving_distance(survey_zmax, **cosmo)*cosmo['h'] 

        if catalog['name'].lower() == 'qpm': 
            # fibercollisions weights are ultimately the weights I will be using here
            fc_mock.weight = fc_mock.wfc            
        n_gal = len(fc_mock.weight)            # Ngal 
        
        upw_index = np.arange(n_gal)[ fc_mock.weight > 1 ]  # index of upweighted galaxy 

        upw_dNN = genv.d_NN(fc_mock.ra[upw_index], fc_mock.dec[upw_index], fc_mock.z[upw_index], 
                n=NN, **fc_cat_corr)       # upweight kth Nearest Neighbor distance
        print 'dNN min= ', np.min(upw_dNN), ', max =', np.max(upw_dNN)
        
        fake_ra, fake_dec, fake_z, fake_weight, re_upw = [], [], [], [], []
        if catalog['name'].lower() == 'qpm': 
            fake_wfkp = [] 
            fake_comp = [] 

        for ii_upw, i_upw in enumerate(upw_index):     # go through upweighted galaxy catalog 

            # for galaxies with wcp > 1
            while fc_mock.weight[i_upw] > 1: 
                fc_mock.weight[i_upw] = fc_mock.weight[i_upw]-1.0

                # LOS comoving distance of the galaxy 
                comdis_iupw = \
                        cosmos.distance.comoving_distance(fc_mock.z[i_upw], **cosmo)*cosmo['h']
                
                # fpeak given dNN of galaxy 
                fpeak_iupw = fc_dlos.fpeak_dNN( upw_dNN[ii_upw], n=NN, **fc_cat_corr )      

                rand_num = np.random.random(1)          
                if rand_num < fpeak_iupw:          # if in the peak 
                
                    # keep RA/Dec w = 1
                    fake_ra.append(fc_mock.ra[i_upw])
                    fake_dec.append(fc_mock.dec[i_upw])
                    fake_weight.append(1.0)

                    if catalog['name'].lower() == 'qpm': 
                        fake_wfkp.append(fc_mock.wfkp[i_upw]) 
                        fake_comp.append(fc_mock.comp[i_upw])
                    
                    if correction['fit'].lower() in ('gauss', 'expon'): 
                        # compute the displacement within peak
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*correction['sigma']
                            peakpofr = fit_func(rand2, correction['sigma']) 

                    # in case the displacement falls out of bound 
                    if (comdis_iupw + rand2 > survey_comdis_min) or \
                            (comdis_iupw + rand2 < survey_comdis_max): 

                        collided_z = comdis2z(comdis_iupw - rand2, **cosmo)
                    else: 
                        collided_z = comdis2z(comdis_iupw + rand2, **cosmo)
                    
                    fake_z.append(collided_z[0]) 
                    
                else:       # Tail ------------------------------------------------
                    if correction['name'] in ('peakshot_dnn'): 
                        # peak + shot noise correction 
                        # record galaxy index that need to be upweighted again 
                        re_upw.append(i_upw)
                    else: 
                        raise NotImplementedError

        if correction['name'] in ('peakshot_dnn'): 
            # re-upweight for peak+shotnoise correction 
            for i_reupw in re_upw: 
                fc_mock.weight[i_reupw] = fc_mock.weight[i_reupw]+1.
        else:
            raise NotImplementedError

        print len(fake_ra), ' galaxies were peak corrected'
        fc_mock.ra = np.concatenate([fc_mock.ra, fake_ra])
        fc_mock.dec = np.concatenate([fc_mock.dec, fake_dec])
        fc_mock.weight = np.concatenate([fc_mock.weight, fake_weight])
        if catalog['name'].lower() == 'qpm': 
            fc_mock.wfkp = np.concatenate([fc_mock.wfkp, fake_wfkp])
            fc_mock.comp = np.concatenate([fc_mock.comp, fake_comp])
        fc_mock.z = np.concatenate([fc_mock.z, fake_z])

    else: 
        raise NotImplementedError('Error here') 

    peakcorr = galaxy_data('data', readdata=False, **cat_corr)
    peakcorr_file = peakcorr.file_name 
    
    if catalog['name'].lower() == 'qpm': 
        # QPM has specific formatting
        np.savetxt(peakcorr_file, 
                np.c_[fc_mock.ra, fc_mock.dec, fc_mock.z, 
                    fc_mock.wfkp, fc_mock.weight, fc_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 
    else: 
        np.savetxt(peakcorr_file, np.c_[fc_mock.ra, fc_mock.dec, fc_mock.z, fc_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def comdis2z(comdis, **cosmo): 
    '''
    Given comoving distance and cosmology, determine z 
    Comoving distance *has* to be in Mpc/h
    '''
    z_arr = np.array([0.0+0.05*np.float(i) for i in range(21)]) 
    dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo)*cosmo['h']

    # use numpy interpolate
    z = np.interp(comdis, dm_arr, z_arr) 
    return z 
