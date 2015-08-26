'''

Code to handle galaxy data for FiberCollisions 

Author(s): ChangHoon Hahn 


'''

import numpy as np
import scipy as sp 
import time 
import random
import os.path
import subprocess
import cosmolopy as cosmos
import warnings 
import matplotlib.pyplot as plt

# --- Local ----
import fibcol_dlos as fc_dlos
import fibcol_nbar as fc_nbar
import photoz as photoz 
import galaxy_environment as genv
import pyspherematch as pysph
from utility.fitstables import mrdfits

# Classes ------------------------------------------------------------
class galaxy_data: 
    def __init__(self, DorR, clobber=False, **cat_corr): 
        ''' Given cat_corr dictionary read in/ build the data or random 
        file and store all the appropriate values  

        Parameters
        ----------
        DorR : 'data' or 'random'
        cat_corr :  Catalog correction Dictionary 
        cosmology : 'fiducial' uses fiducial cosmology; 'boxcosmo' uses box cosmology
        clobber : 

        '''
        self.cat_corr = cat_corr                    # save catalog and correction metadata 
        catalog = cat_corr['catalog'] 
        correction = cat_corr['correction'] 

        file_name = get_galaxy_data_file(DorR, **cat_corr)              # get file name 
        self.file_name = file_name          # save file name to class 

        if catalog['name'].lower() == 'lasdamasgeo':    # LasDamas Geo -----------------------
            omega_m = 0.25  # cosmology

            if DorR == 'data':                      # Data -------------------------
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

                    elif 'scratch' in correction['name'].lower():           
                        # scratch pad for different methods 

                        build_ldg_scratch(**cat_corr) 

                    elif correction['name'].lower() in ('bigfc'): 
                        # fiber collided with bigger angular scale
                        build_ldg_bigfc(**cat_corr)

                    elif correction['name'].lower() in ('bigfc_peakshot'): 
                        # fiber collided with bigger angular scale
                        build_peakcorrected_fibcol(doublecheck=True, **cat_corr)

                    else: 
                        raise NameError('Correction Name Unknown') 

                # read data (ra, dec, z, weights) 
                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         

                for i_col, catalog_column in enumerate(catalog_columns): 
                    column_data = file_data[i_col]
                    setattr(self, catalog_column, column_data) 

            elif DorR.lower() == 'random':          # Random  -------------------------
                # columns of random catalog (NOTE CZ IS CONVERTED TO Z) 
                catalog_columns = ['ra', 'dec', 'z']        

                self.columns = catalog_columns
                
                if 'down_nz' in correction['name'].lower(): 
                    # true rnadom catalog downsampled by nbar(z) 
                    build_ldg_nz_down('random', **cat_corr)
                    
                    # Read data 
                    file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])         

                    for i_col, catalog_column in enumerate(catalog_columns): 
                        column_data = file_data[i_col]
                        # assign to class
                        setattr(self, catalog_column, column_data)

                else: 
                    # Read data 
                    file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])         

                    for i_col, catalog_column in enumerate(catalog_columns): 
                        if catalog_column == 'z': 
                            column_data = file_data[i_col]/299800.0
                        else: 
                            column_data = file_data[i_col]
                        # assign to class
                        setattr(self, catalog_column, column_data)
        
        elif catalog['name'].lower() == 'ldgdownnz':    # LasDamasGeo downsampled ------------ 
            omega_m = 0.25  # cosmology

            if DorR == 'data':                          # Data -------------------------
                # columns that this catalog data will have  
                catalog_columns = ['ra', 'dec', 'z', 'weight']  

                self.columns = catalog_columns

                if (os.path.isfile(file_name) == False) or (clobber == True): 
                    # if file does not exists, make file  
                    print 'Building', file_name

                    if correction['name'].lower() == 'true':    # True
                        # true mock catalog downsampled by nbar(z) 
                        build_ldg_nz_down('data', **cat_corr)
                    elif correction['name'].lower() in ('upweight'): 
                        # nearest neighbor weights assigned
                        build_fibercollided(**cat_corr)
                    elif correction['name'].lower() in ('bigfc'): 
                        # fiber collided with bigger angular scale
                        build_ldgdownnz_bigfc(**cat_corr)
                    elif correction['name'].lower() in ('peakshot'): 
                        # peak + shotnoise correction for monopole
                        build_peakcorrected_fibcol(doublecheck=True, **cat_corr)

                    elif correction['name'].lower() in ('bigfc_peakshot'): 
                        # fiber collided with bigger angular scale
                        build_peakcorrected_fibcol(doublecheck=True, **cat_corr)
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
                
                if (not os.path.isfile(file_name)) or clobber: 
                    # true rnadom catalog downsampled by nbar(z) 
                    build_ldg_nz_down('random', **cat_corr)
                else: 
                    pass 
                    
                # Read data 
                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])         

                for i_col, catalog_column in enumerate(catalog_columns): 
                    column_data = file_data[i_col]
                    # assign to class
                    setattr(self, catalog_column, column_data)

        elif catalog['name'].lower() == 'tilingmock':       # Tiling Mock ----------------------
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
                
        elif catalog['name'].lower() == 'qpm':                      # QPM --------------------
            omega_m = 0.31  # survey cosmology 

            if DorR == 'data':                  # Data ------------------------------
                # catalog columns 
                catalog_columns = ['ra', 'dec', 'z', 'wfc', 'comp'] 
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
                        build_peakcorrected_fibcol(doublecheck=True, **cat_corr)  # build peak corrected file 

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

                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4])   # ra, dec, z, wfc, comp

                # assign to data columns class
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col]) 
            
            elif DorR == 'random':              # Random ------------------------------------

                catalog_columns = ['ra', 'dec', 'z', 'comp']    # catalog columns 

                if (os.path.isfile(file_name) == False) or (clobber == True):
                    print 'Constructing ', file_name 
                    build_random(**cat_corr) 

                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])   # ra, dec, z, comp 

                # assign to object data columns
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col])
        
        elif catalog['name'].lower() == 'nseries':                  # N-series ----------------
            try: 
                if catalog['cosmology'] == 'fiducial': 
                    omega_m = 0.31 # fiducial cosmology  
                else: 
                    omega_m = 0.286  # survey cosmology 

            except KeyError: 
                omega_m = 0.31 # fiducial cosmology  

            if DorR == 'data':                          # Data ------------------------------
                # catalog columns 
                # ra, dec, z, wfc, comp
                catalog_columns = ['ra', 'dec', 'z', 'wfc', 'comp'] 
                self.columns = catalog_columns
                column_indices = [0,1,2,3,4]
                dtypes = None

                if not os.path.isfile(file_name) or clobber:
                    # File does not exist or Clobber = True!
                    print 'Constructing ', file_name 

                    if correction['name'].lower() == 'true':                    
                        # true mocks 
                        # all weights = 1 (fibercollisions *not* imposed) 
                        build_true(**cat_corr) 

                    elif correction['name'].lower() in (
                            'upweight', 'shotnoise', 'floriansn', 'hectorsn'): 
                        # upweighted mocks
                        build_fibercollided(**cat_corr) 

                    elif correction['name'].lower() in ('peaknbar', 'peakshot'): 
                        # peak corrected mocks 
                        build_peakcorrected_fibcol(doublecheck=True, **cat_corr)  # build peak corrected file 
                        #elif correction['name'].lower() in ('noweight'): 
                        #    # n oweight 
                        #    build_noweight(**cat_corr) 
                    elif 'scratch' in correction['name'].lower():           
                        # scratch pad for different methods 
                        build_nseries_scratch(**cat_corr) 
                    elif correction['name'].lower() == 'photoz': 
                        # photoz assigned fiber collided mock 
                        photoz.build_fibcol_assign_photoz(qaplot=False, **cat_corr) 
                    elif correction['name'].lower() == 'photozpeakshot': 
                        # Peak Shot correction using photometric redshift information
                        build_photoz_peakcorrected_fibcol(doublecheck=True, **cat_corr)
                    else: 
                        raise NotImplementedError() 

                if 'photoz' in correction['name'].lower(): 
                    # corrections to column assignment for methods involve photometry
                    catalog_columns = ['ra', 'dec', 'z', 'wfc', 'comp', 
                            'zupw', 'upw_index', 'z_photo'] 
                    self.columns = catalog_columns
                    column_indices = [0,1,2,3,4,5,6,7]
                    dtypes = {
                            'names': ('ra', 'dec', 'z', 
                                'wfc', 'comp', 'zupw', 'upw_index', 'z_photo'), 
                            'formats': (np.float, np.float, np.float, 
                                np.float, np.float, np.float, np.int, np.float)
                            }

                file_data = np.loadtxt(file_name, 
                        unpack=True, usecols=column_indices, dtype=dtypes)         

                # assign to data columns class
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col]) 
            
            elif DorR == 'random':              # Random ------------------------------------

                catalog_columns = ['ra', 'dec', 'z', 'comp']    # catalog columns 
                
                if not os.path.isfile(file_name) or clobber:
                    print 'Constructing ', file_name 
                    build_random(**cat_corr) 

                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])   # ra, dec, z, comp

                # assign to object data columns
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col])

        elif catalog['name'].lower() == 'patchy':                   # PATCHY ------------------
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
    
        elif 'bigmd' in catalog['name'].lower():              # Big MultiDark ----------------
            try: 
                if catalog['cosmology'] == 'fiducial': 
                    omega_m = 0.31      # (fiducial cosmology) 
                else:
                    raise NotImplementedError('Cosmology Type is not Implemented Yet') 

            except KeyError: 
                catalog['cosmology'] = 'fiducial'
                omega_m = 0.31      # (fiducial cosmology) 

            if DorR == 'data':                          # Data ------------------------------
                # catalog columns 
                catalog_columns = ['ra', 'dec', 'z', 'nbar', 'wfc'] 
                self.columns = catalog_columns

                if not os.path.isfile(file_name) or clobber:
                    # File does not exist or Clobber = True!
                    print 'Constructing ', file_name 

                    if correction['name'].lower() == 'true':  
                        # true mocks 
                        # all weights = 1 (fibercollisions *not* imposed) 
                        build_true(**cat_corr) 

                    elif correction['name'].lower() in ('upweight'): 
                        # upweighted mocks
                        build_fibercollided(**cat_corr) 

                    elif correction['name'].lower() in ('peakshot'): 
                        # peak corrected
                        build_peakcorrected_fibcol(doublecheck=True, **cat_corr)

                    else: 
                        raise NotImplementedError() 

                # ra, dec, z, nbar, wfc
                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4])

                # assign to data columns class
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col]) 
            
            elif DorR == 'random':              # Random ------------------------------------
                catalog_columns = ['ra', 'dec', 'z']    # catalog columns 
                
                if not os.path.isfile(file_name) or clobber:
                    print 'Constructing ', file_name 
                    build_random(**cat_corr) 

                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])   # ra, dec, z, comp

                # assign to object data columns
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col])

        elif 'cmass' in catalog['name'].lower():                    # CMASS -----------------
            # and all CMASS like catalogs including 
            # 'cmasslowz_high' and 'cmasslowz_low' from CMASS LOWZ combined sample 
            try: 
                if catalog['cosmology'] == 'fiducial': 
                    # Firm shift in cosmology use to OmegaM = 0.31 in WG 
                    omega_m = 0.31      # (fiducial cosmology) 
                else:
                    raise NotImplementedError('You should use fiducial cosmology') 

            except KeyError: 
                catalog['cosmology'] = 'fiducial'
                omega_m = 0.31      # (fiducial cosmology) 

            if DorR == 'data':              # data --------------------
                catalog_columns = ['ra', 'dec', 'z', 'nbar', 'wsys', 'wnoz', 'wfc', 'comp']
                self.columns = catalog_columns
                
                if not os.path.isfile(file_name) or clobber:
                    # Construct file 
                    print 'Constructing ', file_name 
                    
                    if correction['name'].lower() in ('upweight'): 
                        # upweighted
                        build_fibercollided(**cat_corr) 
                    elif correction['name'].lower() in ('peakshot'): 
                        # peak correction  
                        build_peakcorrected_fibcol(doublecheck=True, **cat_corr) 
                    else: 
                        raise NotImplementedError('Only upweight works for now') 

                #ra,dec,z,wsys,wnoz,wfc,nbar,comp
                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5,6,7])   

                # assign to data columns class
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col]) 

            elif DorR == 'random': 
                catalog_columns = ['ra', 'dec', 'z', 'nbar', 'comp']
                self.columns = catalog_columns
                
                if not os.path.isfile(file_name) or clobber: 
                    print 'Constructing ', file_name 
                    build_random(**cat_corr) 

                #ra,dec,z,nbar,comp
                file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4])  

                # assign to data columns class
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col]) 
    
        else: 
            raise NameError('not yet coded') 
        # survey cosmology metadata 
        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.676
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
    
    if catalog['name'].lower() == 'lasdamasgeo':                # LasDamasGeo --------------
        data_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'

        if DorR.lower() == 'data':                  # data

            if correction['name'].lower() == 'true':    # true
                file_name = ''.join([data_dir, 
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], 
                    '_no.rdcz.dat']) 

            elif correction['name'].lower() in ('upweight', 'fibcol', 
                    'shotnoise', 'hectorsn', 'floriansn'): 
                # LDG mocks with fiber collision weights 
                
                file_name = ''.join([data_dir, 
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], 
                    '_no.rdcz.fibcoll.dat']) 

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
                
                file_name = ''.join([data_dir, 'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.fibcoll.dat'+corr_str]) 
            elif correction['name'].lower() in ('bigfc_peakshot'): 

                # specify peak correction fit (expon, gauss, true) 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    corr_str = ''.join([
                        '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                elif correction['fit'].lower() in ('true'):
                    corr_str = ''.join([
                        '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.fpeak', str(correction['fpeak'])])
                else: 
                    raise NameError('peak fit has to be specified: gauss, expon, true') 
                
                file_name = ''.join([data_dir, 
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], 
                    '_no.rdcz.big_fibcoll.dat'+corr_str]) 

            elif correction['name'].lower() in ('peakshot_dnn'): 
                # Peak Correction (with dNN env) + Shotnoise 
                # Correction needs to specify: fit, nth NN, and sigma 

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

            elif 'scratch' in correction['name'].lower(): # easily adjustable methods 
                file_name = ''.join([data_dir, 'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.fibcoll.', 
                    correction['name'].lower(), '.dat']) 
            
            elif correction['name'].lower() in ('bigfc'): # fc with bigger ang scale
                file_name = ''.join([data_dir, 
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], 
                    '_no.rdcz.big_fibcoll.dat'])

            else: 
                raise NotImplementedError('asdfasdf') 

        if DorR.lower() == 'random':        # random 
            file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat'
    
    elif catalog['name'].lower() == 'ldgdownnz':                # LasDamasGeo downsampled ----- 
        data_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'

        if DorR.lower() == 'data':                  # data

            if correction['name'].lower() == 'true':    # true
                # true mock catalog downsampled by redshift 
                file_name = ''.join([data_dir, 
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], 
                    '_no.rdcz.down_nz.dat'])

            elif correction['name'].lower() in ('upweight'):    # upweight
                # Nearest angular neighbor weights assigned
                file_name = ''.join([data_dir, 
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], 
                    '_no.rdcz.fibcoll.down_nz.dat'])
            
            elif correction['name'].lower() in ('bigfc'): # fc with bigger ang scale
                file_name = ''.join([data_dir, 
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], 
                    '_no.rdcz.big_fibcoll.down_nz.dat'])

            elif correction['name'].lower() in ('peakshot'): 
                # peak + shot noise correction 
                
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    # specify peak correction fit (expon, gauss, true) 
                    corr_str = ''.join([
                        '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), 
                        '.fpeak', str(correction['fpeak'])])
                else: 
                    raise NameError('peak fit has to be specified: gauss, expon, true') 
                
                file_name = ''.join([data_dir, 
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], 
                    '_no.rdcz.fibcoll.down_nz.dat'+corr_str]) 

            elif correction['name'].lower() in ('bigfc_peakshot'): 
                # peak + shot noise correction 
                
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    # specify peak correction fit (expon, gauss, true) 
                    corr_str = ''.join([
                        '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), 
                        '.fpeak', str(correction['fpeak'])])
                else: 
                    raise NameError('peak fit has to be specified: gauss, expon, true') 
                
                file_name = ''.join([data_dir, 
                    'sdssmock_gamma_lrgFull_zm_oriana', 
                    str("%02d" % catalog['n_mock']), catalog['letter'], 
                    '_no.rdcz.big_fibcoll.down_nz.dat'+corr_str]) 

            else: 
                raise NotImplementedError('asdfasdf') 

        if DorR.lower() == 'random': 

            file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/',
                'sdssmock_gamma_lrgFull.rand_200x_no.rdcz.down_nz.dat']) 

    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ---------------
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
                corr_str = ''.join([
                    '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                    '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                file_name = data_dir+'cmass-boss5003sector-icoll012.zlim.dat'+corr_str 
                
            # Fibercollisions Corrected by all peak correction ------------------------------
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

        elif DorR.lower() == 'random':      # Randoms 
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

    elif catalog['name'].lower() == 'patchy':               # PATHCY mocks -------------
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
    
    elif catalog['name'].lower() == 'nseries':              # N-series ------------------
        data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
        
        if DorR == 'data':  # mock catalogs 
            try: 
                if catalog['cosmology'] == 'fiducial': 
                    cosmo_str = '_fidcosmo' 
                else: 
                    raise NotImplementeError('lkajsdf') 
            except KeyError: 
                cosmo_str = '_fidcosmo' 

            if correction['name'].lower() in ('true'):  # true mocks
                file_name = ''.join([data_dir, 
                    'CutskyN', str(catalog['n_mock']), '.dat']) 
            elif correction['name'].lower() in ('original'):    # original mock
                file_name = ''.join([data_dir, 
                    'CutskyN', str(catalog['n_mock']), '.rdzwc']) 
            elif correction['name'].lower() in ('wcompfile'):   # wcomp file 
                file_name = ''.join([data_dir, 
                    'CutskyN', str(catalog['n_mock']), '.mask_info']) 
            elif correction['name'].lower() in (
                    'upweight', 'shotnoise', 'floriansn', 'hectorsn'): # upweighted mocks 
                file_name = ''.join([data_dir, 
                    'CutskyN', str(catalog['n_mock']), '.fibcoll.dat']) 
            elif correction['name'].lower() in ('peaknbar', 'peakshot'): # peak corrected mocks 
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
                    'CutskyN', str(catalog['n_mock']), '.fibcoll', corr_str, cosmo_str, '.dat' ]) 
            elif 'scratch' in correction['name'].lower():   # scratch pad mocks 
                
                file_name = ''.join([data_dir, 
                    'CutskyN', str(catalog['n_mock']), '.fibcoll.', correction['name'].lower(), cosmo_str, '.dat' ]) 
            elif correction['name'].lower() == 'photoz':    # photoz assigned fiber collided mock  
                file_name = ''.join([data_dir, 
                    'CutskyN', str(catalog['n_mock']), 
                    '.fibcoll.', correction['name'].lower(), '.dat' ]) 
            elif correction['name'].lower() == 'photozpeakshot':    # peak shot correction utilizing photoz
                
                # specify best fit (expon or gauss) function to peak 
                if correction['fit'].lower() in ('gauss'): 
                    # correction specifier string 
                    corr_str = ''.join(['.photoz.', correction['fit'].lower(), 
                        '.peakshot.sigma', str(correction['sigma']), 
                        '.fpeak', str(correction['fpeak'])]) 
                else: 
                    raise NotImplementedError('peak fit has to be specified: gauss or expon') 

                file_name = ''.join([data_dir, 
                    'CutskyN', str(catalog['n_mock']), '.fibcoll', corr_str, cosmo_str, '.dat' ]) 
            else: 
                raise NameError('not yet coded') 

        elif DorR == 'random': 
            # random catalog 
            file_name = ''.join([data_dir, 'Nseries_cutsky_randoms_50x_redshifts_comp.dat']) 
    
    elif 'bigmd' in catalog['name'].lower():                # Big MD ---------------------
        data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
        
        if DorR == 'data':  # mock catalogs 
            if correction['name'].lower() == 'true':    # true mocks
                if catalog['name'].lower() == 'bigmd': 
                    file_name = ''.join([data_dir, 
                        'bigMD-cmass-dr12v4_vetoed.dat']) # hardcoded
                elif catalog['name'].lower() == 'bigmd1': 
                    file_name = ''.join([data_dir, 
                        'bigMD-cmass-dr12v4-RST-standardHAM_vetoed.dat']) # hardcoded
                elif catalog['name'].lower() == 'bigmd2': 
                    file_name = ''.join([data_dir, 
                        'bigMD-cmass-dr12v4-RST-quadru_vetoed.dat']) # hardcoded
                elif catalog['name'].lower() == 'bigmd3':   
                    # "best" bigMD August 3, 2015 
                    file_name = ''.join([data_dir, 
                        'BigMD-cmass-dr12v4-RST-standHAM-Vpeak_vetoed.dat']) 
                else: 
                    raise NotImplementedError('asdfkj')
    
            elif correction['name'].lower() in ('upweight'):    # upweighted
                if catalog['name'].lower() == 'bigmd':
                    file_name = ''.join([data_dir, 
                        'bigMD-cmass-dr12v4_vetoed.fibcoll.dat']) # hardcoded
                elif catalog['name'].lower() == 'bigmd1': 
                    file_name = ''.join([data_dir, 
                        'bigMD-cmass-dr12v4-RST-standardHAM_vetoed.fibcoll.dat']) 
                elif catalog['name'].lower() == 'bigmd2': 
                    file_name = ''.join([data_dir, 
                        'bigMD-cmass-dr12v4-RST-quadru_vetoed.fibcoll.dat'])                 
                elif catalog['name'].lower() == 'bigmd3': 
                    file_name = ''.join([data_dir, 
                        'BigMD-cmass-dr12v4-RST-standHAM-Vpeak_vetoed.fibcoll.dat']) 
                else: 
                    raise NotImplementedError('asdfkj')
            
            elif correction['name'].lower() in ('peakshot'):    # peak + shotnoise 

                corr_str = ''.join([
                    '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                    '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                if catalog['name'].lower() == 'bigmd':
                    file_name = ''.join([data_dir, 
                        'bigMD-cmass-dr12v4_vetoed.fibcoll.dat', corr_str]) # hardcoded
                elif catalog['name'].lower() == 'bigmd1': 
                    file_name = ''.join([data_dir, 
                        'bigMD-cmass-dr12v4-RST-standardHAM_vetoed.fibcoll.dat', corr_str]) 
                elif catalog['name'].lower() == 'bigmd2': 
                    file_name = ''.join([data_dir, 
                        'bigMD-cmass-dr12v4-RST-quadru_vetoed.fibcoll.dat', corr_str]) 
                elif catalog['name'].lower() == 'bigmd3': 
                    file_name = ''.join([data_dir, 
                        'BigMD-cmass-dr12v4-RST-standHAM-Vpeak_vetoed.fibcoll.dat', corr_str]) 
                else: 
                    raise NotImplementedError('asdfkj')
            else: 
                raise NotImplementedError('not yet coded') 
    
        elif DorR == 'random':                  # random catalog 
            # vetomask-ed random catalog 
            if catalog['name'].lower() == 'bigmd': 
                file_name = ''.join([data_dir, 'bigMD-cmass-dr12v4_vetoed.ran'])
            elif catalog['name'].lower() == 'bigmd1': 
                file_name = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-standardHAM_vetoed.ran'])
            elif catalog['name'].lower() == 'bigmd2': 
                file_name = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-quadru_vetoed.ran'])
            elif catalog['name'].lower() == 'bigmd3': 
                file_name = ''.join([data_dir, 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak_vetoed.ran']) 
            else: 
                raise NotImplementedError('asdlfkj')
    
    elif 'cmass' in catalog['name'].lower():                # CMASS ---------------------
        data_dir = '/mount/riachuelo1/hahn/data/CMASS/'
        
        if DorR == 'data':  # mock catalogs 
            if catalog['name'].lower() == 'cmass': 
                if correction['name'].lower() in ('upweight'):      
                    # upweighted
                    file_name = ''.join([data_dir, 
                        'cmass-dr12v4-N-Reid.dat']) # hardcoded
                elif correction['name'].lower() in ('peakshot'):    
                    # peakshot
                    # correction string in file name 
                    corr_str = ''.join([
                        '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                    file_name = ''.join([data_dir, 
                        'cmass-dr12v4-N-Reid.dat', corr_str]) 
                else: 
                    raise NotImplementedError('Only upweight and peakshot implemented so far') 
            elif catalog['name'].lower() == 'cmasslowz_high': 
                # CMASS LOWZ combined sample high redshift bin 
                if correction['name'].lower() in ('upweight'):      
                    # upweighted
                    file_name = ''.join([data_dir, 
                        'cmasslowz-dr12v4-N-Reid-high.dat']) # hardcoded
                elif correction['name'].lower() in ('peakshot'):    
                    # peakshot
                    # correction string in file name 
                    corr_str = ''.join([
                        '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                    file_name = ''.join([data_dir, 
                        'cmasslowz-dr12v4-N-Reid-high.dat', corr_str]) 
            elif catalog['name'].lower() == 'cmasslowz_low':
                # CMASS LOWZ combined sample low redshift bin 
                if correction['name'].lower() in ('upweight'):      
                    # upweighted
                    file_name = ''.join([data_dir, 
                        'cmasslowz-dr12v4-N-Reid-low.dat']) # hardcoded
                elif correction['name'].lower() in ('peakshot'):    
                    # peakshot
                    # correction string in file name 
                    corr_str = ''.join([
                        '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                    file_name = ''.join([data_dir, 
                        'cmasslowz-dr12v4-N-Reid-low.dat', corr_str]) 
    
        elif DorR == 'random':                  # random catalog 
            if catalog['name'].lower() == 'cmass': 
                file_name = ''.join([data_dir, 
                    'cmass-dr12v4-N-Reid.ran.dat'])
            elif catalog['name'].lower() == 'cmasslowz_high': 
                file_name = ''.join([data_dir, 
                    'cmasslowz-dr12v4-N-Reid-high.ran.dat'])
            elif catalog['name'].lower() == 'cmasslowz_low': 
                file_name = ''.join([data_dir, 
                    'cmasslowz-dr12v4-N-Reid-low.ran.dat'])
            else: 
                raise NotImplementedError('lskdfjaklsdfj')

    return file_name 

# ------------------------------------------------------------------------
# Build galaxy data  
# ------------------------------------------------------------------------

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

    elif catalog['name'].lower() == 'qpm':                                      # QPM ------------------------------------
        P0 = 20000.0    # hardcoded P0 value

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
        true_wfkp = orig_true_data[:,3]
        true_wfc = np.array([1.0 for i in range(n_gal)])    # fiber collisions weights are all 1 for true

        # check to make sure that the redshifts correspond btw rdz file and rdz.info file 
        if np.max(np.abs(orig_true_info[:,3]/true_z-1.0)) > 10**-5: 
            raise NameError('redshifts dont match') 

        true_comp = orig_true_info[:,1]         # completness weights

        # remove veto mask 
        vetomask = (orig_true_veto == 0)            
        # Only keep galaxies with veto = 0 (for veto values in .veto file) 
        
        true_file = get_galaxy_data_file('data', readdata=False, **cat_corr)
        np.savetxt(true_file, np.c_[
            true_ra[vetomask], true_dec[vetomask], true_z[vetomask], 
            true_wfc[vetomask], true_comp[vetomask]], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    elif catalog['name'].lower() == 'nseries':              # N Series ---------------------------------
        P0 = 20000.0

        # read rdzw file 
        data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
        orig_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.rdzwc']) 
        orig_ra, orig_dec, orig_z, orig_wfc = np.loadtxt(orig_file, unpack=True, usecols=[0,1,2,4])
    
        # file with completeness
        mask_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.mask_info']) 
        orig_wcomp = np.loadtxt(mask_file, unpack=True, usecols=[0]) 

        # true wfc 
        true_wfc = np.array([ 1.0 for i in range(len(orig_wfc)) ]) 
        
        # write to file 
        true_file = get_galaxy_data_file('data', **cat_corr) 
        np.savetxt(true_file, 
                np.c_[orig_ra, orig_dec, orig_z, true_wfc, orig_wcomp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    elif catalog['name'].lower() == 'lasdamasgeo':          # Las Damas Geo -----------------------------
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

    elif catalog['name'].lower() == 'patchy':               # PATCHY mocks ------------------------
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
    
    elif 'bigmd' in catalog['name'].lower():                # Big MultiDark ------------
        P0 = 20000.0

        # read rdzw file 
        data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
        if catalog['name'].lower() == 'bigmd':
            orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-wcp-veto.dat'])  # hardcoded
        elif catalog['name'].lower() == 'bigmd1': 
            orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-standardHAM-veto.dat'])
        elif catalog['name'].lower() == 'bigmd2': 
            orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-quadru-veto.dat'])
        elif catalog['name'].lower() == 'bigmd3': 
            orig_file = ''.join([data_dir, 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak-veto.dat']) 
        else: 
            raise NotImplementedError('asdlfkjadf') 

        orig_ra, orig_dec, orig_z, orig_wfkp, orig_veto, orig_wfc = np.loadtxt(orig_file, 
                unpack=True, usecols=[0,1,2,3,4,5])

        # true wfc = 1 for all galaxies 
        true_wfc = np.array([ 1.0 for i in range(len(orig_wfc)) ]) 
        nbar = (1.0/P0) * (1.0/orig_wfkp - 1.0) 

        vetomask = np.where(orig_veto == 1)     # if veto = 1 then keep; otherwise discard 
        
        # write to file 
        true_file = get_galaxy_data_file('data', **cat_corr) 
        np.savetxt(true_file, 
                np.c_[orig_ra[vetomask], orig_dec[vetomask], orig_z[vetomask], 
                    nbar[vetomask], true_wfc[vetomask]], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f'], delimiter='\t') 

    else: 
        raise NameError('not yet coded') 

def build_random(**cat_corr): 
    ''' Build the random catalogs from original data 

    Paramter
    --------
    cat_corr : Catalog and Correction dictionary
    '''
    catalog = cat_corr['catalog']
    if 'cmass' in catalog['name'].lower():          # CMASS -------------------------------- 
        data_dir = '/mount/riachuelo1/hahn/data/CMASS/'

        if catalog['name'].lower() == 'cmass': 
            # random data fits file
            data_file = ''.join([data_dir, 'cmass-dr12v4-N-Reid.ran.fits']) 
            cmass = mrdfits(data_file) 
        
            # mask file 
            mask_file = ''.join([data_dir, 'mask-cmass-dr12v4-N-Reid.fits']) 
            mask = mrdfits(mask_file) 
            ipoly = cmass.ipoly # polygon index
            comp = mask.weight[ipoly]
        
            # redshift limit 
            zlimit = np.where((cmass.z >= 0.43) & (cmass.z <= 0.7))

        elif catalog['name'].lower() == 'cmasslowz_high': 
            # random data fits file
            data_file = ''.join([data_dir, 'cmasslowz-dr12v4-N-Reid.ran.fits']) 
            cmass = mrdfits(data_file) 
        
            # mask file 
            mask_file = ''.join([data_dir, 'mask-cmasslowz-dr12v4-N-Reid.fits']) 
            mask = mrdfits(mask_file) 
            ipoly = cmass.ipoly # polygon index
            comp = mask.weight[ipoly]
        
            # redshift limit 
            zlimit = np.where((cmass.z >= 0.5) & (cmass.z < 0.75))

        elif catalog['name'].lower() == 'cmasslowz_low': 
            # random data fits file
            data_file = ''.join([data_dir, 'cmasslowz-dr12v4-N-Reid.ran.fits']) 
            cmass = mrdfits(data_file) 
        
            # mask file 
            mask_file = ''.join([data_dir, 'mask-cmasslowz-dr12v4-N-Reid.fits']) 
            mask = mrdfits(mask_file) 
            ipoly = cmass.ipoly # polygon index
            comp = mask.weight[ipoly]
        
            # redshift limit 
            zlimit = np.where((cmass.z >= 0.2) & (cmass.z < 0.5))
        else: 
            raise NotImplementedError("Only CMASS and CMASS+LOWZ combined sample implemented") 

        #ra, dec, z, nz, comp 
        random_file = get_galaxy_data_file('random', **cat_corr) 
        np.savetxt(random_file, 
                np.c_[
                    (cmass.ra)[zlimit], (cmass.dec)[zlimit], (cmass.z)[zlimit], 
                    (cmass.nz)[zlimit], comp[zlimit]], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f'], delimiter='\t') 

    elif catalog['name'].lower() == 'qpm':            # QPM ------------------------------
        # read original random catalog  
        data_dir = '/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/randoms/'
        orig_true_random = np.loadtxt(''.join([data_dir, 'a0.6452_rand50x.dr12d_cmass_ngc.rdz']))   # ra, dec, z, wfkp
        orig_true_random_info = np.loadtxt(data_dir+'a0.6452_rand50x.dr12d_cmass_ngc.rdz.info')   # galid, comp?
        orig_true_random_veto = np.loadtxt(data_dir+'a0.6452_rand50x.dr12d_cmass_ngc.veto')       # veto  

        vetomask = (orig_true_random_veto == 0)
        true_random_file = get_galaxy_data_file('random', **{'catalog':{'name':'qpm'}, 'correction':{'name':'true'}})
        
        np.savetxt(true_random_file, 
                np.c_[
                    (orig_true_random[:,0])[vetomask], (orig_true_random[:,1])[vetomask], 
                    (orig_true_random[:,2])[vetomask], (orig_true_random_info[:,1])[vetomask]], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    elif catalog['name'].lower() == 'nseries':      # Nseries ----------------------------
        # read original random catalog 
        data_dir = '/mount/riachuelo1/hahn/data/Nseries/'

        orig_rand_file = ''.join([data_dir, 'Nseries_cutsky_randoms_50x_redshifts.dat']) 
        # RA, Decl, Redhsift, w_fkp
        ra, dec, z = np.loadtxt(orig_rand_file, unpack=True, usecols=[0,1,2]) 
    
        # sector completeness
        orig_comp_file = ''.join([data_dir, 'Nseries_cutsky_randoms_50x_maskinfo.dat'])  
        comp = np.loadtxt(orig_comp_file, unpack=True, usecols=[0])
        
        # save rnadom file 
        true_random_file = get_galaxy_data_file('random', **{'catalog':{'name':'nseries'}, 'correction':{'name':'true'}})
        np.savetxt(true_random_file, 
                np.c_[ra, dec, z, comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

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

    elif 'bigmd' in catalog['name'].lower():        # Big MD ----------------------------
        P0 = 20000.0
        # original random catalog 
        data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
        if catalog['name'].lower() == 'bigmd': 
            orig_rand_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-wcp-veto.ran']) 
        elif catalog['name'].lower() == 'bigmd1': 
            orig_rand_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-standardHAM-veto.ran']) 
        elif catalog['name'].lower() == 'bigmd2': 
            orig_rand_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-quadru-veto.ran']) 
        elif catalog['name'].lower() == 'bigmd3': 
            orig_rand_file = ''.join([data_dir, 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak-veto.ran'])
        else: 
            raise NameError('catalo does not exist') 

        # RA, Decl, Redhsift, veto  
        ra, dec, z, wfkp, veto = np.loadtxt(orig_rand_file, unpack=True, usecols=[0,1,2,3,4]) 
    
        nbar = (1.0 / P0) * (1.0/wfkp - 1.0)    # nbar(z) 

        vetomask = np.where(veto == 1)  # impose vetomask 
    
        # save rnadom file (write RA, Decl, Redshift) 
        true_random_file = get_galaxy_data_file('random', 
                **{'catalog':{'name':catalog['name'].lower()}, 'correction':{'name':'true'}})
        np.savetxt(true_random_file, 
                np.c_[ra[vetomask], dec[vetomask], z[vetomask], nbar[vetomask]], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e'], delimiter='\t') 
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
    ''' Build Fibercollided mock catalogs using specific idl routines or by using the given fiber collision weights

    Parameters
    ----------
    cat_corr : catalog correction dictionary 

    Notes
    -----
    '''
    catalog = cat_corr['catalog']

    if 'cmass' in catalog['name'].lower():                  # CMASS ---------------------
        data_dir = '/mount/riachuelo1/hahn/data/CMASS/'

        if catalog['name'].lower() == 'cmass': 
            data_file = ''.join([data_dir, 'cmass-dr12v4-N-Reid.dat.fits']) # fits file 
            data = mrdfits(data_file)

            zlimit = np.where((data.z >= 0.43) & (data.z <= 0.7))

        elif catalog['name'].lower() == 'cmasslowz_high': 
            # high redshift bin of the CMASS LOWZ combined sample
            # 0.5 < z < 0.75

            # original combined data sample
            data_file = ''.join([data_dir, 
                'cmasslowz-dr12v4-N-Reid.dat.fits']) 
            data = mrdfits(data_file) 

            zlimit = np.where((data.z >= 0.5) & (data.z < 0.7))  # redshift limit

        elif catalog['name'].lower() == 'cmasslowz_low': 
            # high redshift bin of the CMASS LOWZ combined sample
            # 0.2 < z < 0.5

            # original combined data sample
            data_file = ''.join([data_dir, 
                'cmasslowz-dr12v4-N-Reid.dat.fits']) 
            data = mrdfits(data_file) 

            zlimit = np.where((data.z >= 0.2) & (data.z < 0.5)) # redshift limit 
    
        # save to file 
        fc_file = get_galaxy_data_file('data', **cat_corr)
        np.savetxt(fc_file, 
                np.c_[
                    (data.ra)[zlimit], (data.dec)[zlimit], (data.z)[zlimit], (data.nz)[zlimit],
                    (data.weight_systot)[zlimit], (data.weight_noz)[zlimit], 
                    (data.weight_cp)[zlimit], (data.comp)[zlimit]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', 
                    '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 
        fibcollided_cmd = ''

    elif catalog['name'].lower() == 'lasdamasgeo':        # Las Damas Geo 
        fibcollided_cmd = 'idl -e "ldg_fibcollmock_wcp_assign,'+str(catalog['n_mock'])+", '"+\
                str(catalog['letter'])+"'"+'"'
        print fibcollided_cmd
        os.system(fibcollided_cmd)  # call IDL code 

    elif catalog['name'].lower() == 'ldgdownnz': 
        fibcollided_cmd = 'idl -e "ldgdownnz_wcp_assign,'+str(catalog['n_mock'])+", '"+\
                str(catalog['letter'])+"'"+'"'
        print fibcollided_cmd
        os.system(fibcollided_cmd)  # call IDL code 

    elif catalog['name'].lower() == 'tilingmock':       # Tiling Mock ----------------
        fibcollided_cmd = ' '.join(['idl', '-e', '"', "build_wcp_assign, 'tilingmock'", '"'])
        os.system(fibcollided_cmd) 

    elif catalog['name'].lower() == 'qpm':              # QPM ------------------------
        orig_true_file = ''.join([
            '/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
            'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 
        orig_true_data = np.loadtxt(orig_true_file) 
        
        orig_true_info_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
            'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz.info']) 
        # gal_id, comp, z_real, z_red, mass_halo, flag_sta, id_halo
        orig_true_info_data = np.loadtxt(orig_true_info_file)    

        if catalog['n_mock'] in (44, 46, 52, 53, 54, 56, 61, 707, 756, 794, 819, 831, 835, 838):
            orig_true_veto_file = ''.join(['/mount/riachuelo1/hahn/data/QPM/dr12d/a0.6452_', 
                str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
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

        fc_file = get_galaxy_data_file('data', **cat_corr)
        np.savetxt(fc_file, 
                np.c_[
                    true_ra[vetomask], true_dec[vetomask], true_z[vetomask], 
                    true_wfc[vetomask], true_comp[vetomask]], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

        fibcollided_cmd = ''
    
    elif catalog['name'].lower() == 'nseries':          # N-series --------------------
        data_dir = '/mount/riachuelo1/hahn/data/Nseries/'   # directory

        # original file 
        orig_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.rdzwc']) 
        orig_ra, orig_dec, orig_z, orig_wfc = np.loadtxt(orig_file, unpack=True, usecols=[0,1,2,4])
    
        # file with mask completeness
        mask_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.mask_info']) 
        orig_wcomp = np.loadtxt(mask_file, unpack=True, usecols=[0]) 

        # write to file 
        true_file = get_galaxy_data_file('data', **cat_corr) 
        np.savetxt(true_file, 
                np.c_[orig_ra, orig_dec, orig_z, orig_wfc, orig_wcomp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

        fibcollided_cmd = ''

    elif catalog['name'].lower() == 'patchy':           # PATCHY mocks ---------------
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
    
    elif 'bigmd' in catalog['name'].lower():            # Big MD --------------------
        P0 = 20000.0
        # read original random catalog 
        data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
        if catalog['name'].lower() == 'bigmd': 
            orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-wcp-veto.dat']) 
        elif catalog['name'].lower() == 'bigmd1': 
            orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-standardHAM-veto.dat']) 
        elif catalog['name'].lower() == 'bigmd2': 
            orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-quadru-veto.dat']) 
        elif catalog['name'].lower() == 'bigmd3': 
            orig_file = ''.join([data_dir, 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak-veto.dat']) 
        else: 
            raise NameError('catalog does not exit') 

        # RA, Decl, Redhsift, veto  
        ra, dec, z, wfkp, veto, wfc = np.loadtxt(orig_file, unpack=True, usecols=[0,1,2,3,4,5]) 
        n_gal = len(ra) 

        nbar = (1.0 / P0) * (1.0/wfkp - 1.0) 
        print nbar
        print 'min nz', min(nbar)
        print 'max nz', max(nbar)
        vetomask = np.where(veto == 1)  # impose vetomask 
    
        fc_file = get_galaxy_data_file('data', **cat_corr)    # file name 
        np.savetxt(fc_file, 
                np.c_[
                    ra[vetomask], dec[vetomask], z[vetomask], 
                    nbar[vetomask], wfc[vetomask]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f'], delimiter='\t') 
        
        fibcollided_cmd = ''
    
    else: 
        raise NameError('not yet coded') 

    return fibcollided_cmd 
"""
def build_peakcorrected_fibcol_old(sanitycheck=False, **cat_corr): 
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

    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo -----------------

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

        survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
        survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']

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

    elif catalog['name'].lower() in ('tilingmock', 'qpm', 'patchy', 'nseries'): 
        # Tiling mock, QPM, PATCHY 
        
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits

        if catalog['name'].lower() == 'qpm':
            n_mocks = 100
        elif catalog['name'].lower() == 'nseries': 
            n_mocks = 84
        
        # read in mock with fibercollisions imposed
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr)  
        
        cosmo = fibcoll_mock.cosmo  # survey cosmology 

        survey_comdis_min = cosmos.distance.comoving_distance(survey_zmin, **cosmo)*cosmo['h']         # in units of Mpc/h
        survey_comdis_max = cosmos.distance.comoving_distance(survey_zmax, **cosmo)*cosmo['h']         # in units of Mpc/h
    
        if catalog['name'].lower() != 'tilingmock':
            # only use fiber collision weights
            fibcoll_mock.weight = fibcoll_mock.wfc            
                        
        if correction['name'] == 'peaknbar': 
            # read in the true galaxies for tail portion  
            # the tail portion of the peak corrections will be generated 
            # similar to the mksample procedures for assigning the random catalog redshift 
            true_cat_corr = {'catalog':catalog, 'correction': {'name': 'true'}}
            true_data = galaxy_data('data', **true_cat_corr)
            
            if catalog['name'].lower() in ('qpm',  'patchy', 'nseries'): 
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

        if catalog['name'].lower() in ('qpm', 'nseries'): 
            append_comp = []   # save comp

        reupw = []

        if sanitycheck == True:     # check peak p(r) 
            pr_test = [] 
    
        for i_mock in range(len(fibcoll_mock.weight)):  # go through each galaxy 

            while fibcoll_mock.weight[i_mock] > 1:  # for galaxies with wcp > 1

                fibcoll_mock.weight[i_mock] -= 1.0

                # LOS comoving distance of the galaxy 
                comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], **cosmo)*cosmo['h']
                
                rand_num = np.random.random(1)  # determine whether galaxy is in the peak
               
                if rand_num < f_peak:          # in the peak 

                    append_ra.append(fibcoll_mock.ra[i_mock])    # keep ra and dec
                    append_dec.append(fibcoll_mock.dec[i_mock])

                    if catalog['name'].lower() in ('qpm', 'nseries'): 
                        append_comp.append(fibcoll_mock.comp[i_mock])

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
                        #-------------------------------------------------------------------- 

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

        if catalog['name'].lower() in ('qpm', 'nseries'): 
            fibcoll_mock.comp = np.concatenate([fibcoll_mock.comp, append_comp])

        fibcoll_mock.z = np.concatenate([fibcoll_mock.z, append_z])
        
        if sanitycheck == True:         # output sanitycheck 
            sanitycheck_file = ''.join([
                'peak_corrected_dlos_values_', correction['name'].lower(), '_', 
                correction['fit'].lower(), '.dat']) 
            np.savetxt(sanitycheck_file,
                    np.c_[append_dlos], fmt=['%10.5f']) 
    else: 
        raise NameError('Error here') 

    peakcorr_file = get_galaxy_data_file('data', **cat_corr) 
    
    if catalog['name'].lower() in ('qpm', 'nseries'): 
        np.savetxt(peakcorr_file, 
                np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.weight, fibcoll_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
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
"""
def build_peakcorrected_fibcol(doublecheck=False, **cat_corr): 
    ''' Build peak corrected fibercollided mock catalogs (using cosmolopy) 

    Parameters
    ----------
    cat_corr : Catalog + Correction dictionary
    doublecheck : save dlos values to test it's working  

    Notes
    -----
    * Currently supported peak correction methods: peakshot 
    * nbar(z) interpolation implemented for CMASS like samples  

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # fit functions (using lambda) ------------------------------------------------------------
    if correction['fit'].lower() == 'gauss': 
        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)
    elif correction['fit'].lower() == 'expon': 
        fit_func = lambda x, sig: np.exp(-1.0*np.abs(x)/sig)
    elif correction['fit'].lower() == 'true': 
        pass 
    else: 
        raise NameError('correction fit has to be specified as gauss or expon') 
    
    # survey specific parameters 
    if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'):     # LasDamasGeo 
        n_mocks = 160   # total number of mocks
        survey_zmin, survey_zmax = 0.16, 0.44
    elif catalog['name'].lower() in ('tilingmock', 'qpm', 'patchy', 'nseries'): 
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits
        if catalog['name'].lower() == 'qpm':
            n_mocks = 100
        elif catalog['name'].lower() == 'nseries': 
            n_mocks = 1 
    elif catalog['name'].lower() in ('bigmd'):             # CMASS, BigMD
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits
        n_mocks = 1 
    elif 'cmass' in catalog['name'].lower():             # CMASS like catalogs
        if catalog['name'].lower() == 'cmass': 
            survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits
        elif catalog['name'].lower() == 'cmasslowz_high': 
            survey_zmin, survey_zmax = 0.5, 0.75    
        elif catalog['name'].lower() == 'cmasslowz_low': 
            survey_zmin, survey_zmax = 0.2, 0.5
        else: 
            raise NotImplementedError('CMASS or CMASSLOWZ combined sample')
        n_mocks = 1 
    else: 
        raise NotImplementedError('Mock Catalog not included')

    # read in fiber collided mock  
    if correction['name'] == 'bigfc_peakshot': 
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'bigfc'}}
    else: 
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
    fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr) 
    cosmo = fibcoll_mock.cosmo      # survey comoslogy 
    
    # comoving distnace of z_min and z_max 
    survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
    survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
    # fiber collision weights
    if catalog['name'].lower() not in ('tilingmock', 'lasdamasgeo', 'ldgdownnz'):
        fibcoll_mock.weight = fibcoll_mock.wfc            
                        
    if correction['name'] == 'peaknbar': 
        # read in the true galaxies for tail portion  
        # the tail portion of the peak corrections will be generated 
        # similar to the mksample procedures for assigning the random catalog redshift 
        true_cat_corr = {'catalog':catalog, 'correction': {'name': 'true'}}
        true_data = galaxy_data('data', **true_cat_corr)
        
        if catalog['name'].lower() in ('qpm',  'patchy', 'nseries', 'cmass'): 
            true_weight = true_data.wfc
        else: 
            true_weight = true_data.weight
        true_z = true_data.z

        true_weight_cum = true_weight.cumsum()
        true_weight_max = true_weight_cum.max()
        i_true = np.arange(0,len(true_weight_cum))

    # import nbar(z) file for CMASS-like catalogs 
    if 'cmass' in catalog['name'].lower(): 
        # hardcoded nbar(z) files 
        if catalog['name'].lower() == 'cmass': 
            nbar_file = '/mount/riachuelo1/hahn/data/CMASS/nbar-cmass-dr12v4-N-Reid-om0p31_Pfkp10000.dat'
        elif 'cmasslowz' in catalog['name'].lower(): 
            nbar_file = '/mount/riachuelo1/hahn/data/CMASS/nbar-cmasslowz-dr12v4-N-Reid-om0p31_Pfkp10000.dat'
        # read in nbar(z) file 
        nbar_z, nbar_nbar = np.loadtxt(nbar_file, skiprows=2, unpack=True, usecols=[0, 3]) 
        # nbar(z) interpolation function
        nbarofz = lambda zz: sp.interpolate.interp1d(nbar_z, nbar_nbar, kind='cubic')       

    # peak fraction         
    try: 
        f_peak = correction['fpeak']
    except KeyError: # all peak and all peak shot corrections
        f_peak = 1.0 

    appended_ra, appended_dec, appended_z, appended_weight = [], [], [], []
    upweight_again = []

    if catalog['name'].lower() in ('qpm', 'nseries'): 
        appended_comp = []   # save comp
    elif 'cmass' in catalog['name'].lower(): 
        appended_comp, appended_wsys, appended_wnoz, appended_nbar = [], [], [], [] 
            
    if doublecheck:     # check that the peak p(r) is generated properly
        dlos_values = [] 
    
    for i_mock in range(len(fibcoll_mock.weight)):  # go through every galaxy in fibercollided mock catalog
        # there is a smarter way to do this; however it's not implemented

        while fibcoll_mock.weight[i_mock] > 1:      # for galaxies with wcp > 1

            fibcoll_mock.weight[i_mock] -= 1.0

            # LOS comoving distance of the galaxy 
            comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], **cosmo)*cosmo['h']
                
            rand_num = np.random.random(1) 
            if rand_num <= f_peak:          # if in the peak 
                appended_ra.append(fibcoll_mock.ra[i_mock])     # keep ra and dec
                appended_dec.append(fibcoll_mock.dec[i_mock])
                
                # keep extra columns
                if catalog['name'].lower() in ('qpm', 'nseries'):  
                    appended_comp.append(fibcoll_mock.comp[i_mock]) 
                elif 'cmass' in catalog['name'].lower(): 
                    appended_comp.append(fibcoll_mock.comp[i_mock]) 
                    appended_wsys.append(fibcoll_mock.wsys[i_mock]) 
                    appended_wnoz.append(1.0) 

                if correction['name'].lower() in ('allpeak', 'allpeakshot'): 
                    # appended galaxy has weight of peak fraction 
                    appended_weight.append(correction['fpeak'])
                else: 
                    # appended galaxy has weight of 1.0 
                    appended_weight.append(1.0)

                # compute the displacement within peak using best-fit
                if correction['fit'].lower() in ('gauss', 'expon'):   
                    rand1 = np.random.random(1) 
                    rand2 = np.random.random(1) 

                    rand2 = (-3.0+rand2*6.0)*correction['sigma']
                    peakpofr = fit_func(rand2, correction['sigma']) 
                    
                    while peakpofr <= rand1: 
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 

                elif correction['fit'].lower() == 'true': 
                    # compute the displacement within peak using actual distribution   
                    dlos_comb_peak_file = ''.join([
                        ((fibcoll_mock.file_name).rsplit('/', 1))[0], '/', 
                        'DLOS_norm_peak_dist_', catalog['name'].lower(), '_', str(n_mocks), 'mocks_combined.dat'])
                    dlos_mid, dlos_dist = np.loadtxt(dlos_comb_peak_file, unpack=True, usecols=[0,1])

                    dlos_cdf = dlos_dist.cumsum()/dlos_dist.sum()

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
                else: 
                    raise NotImplementedError('asdfasdf')

                # in case the displacement falls out of bound (may general large scale issues)
                if (comdis_imock + rand2 > survey_comdis_max) or (comdis_imock + rand2 < survey_comdis_min): 
                    collided_z = comdis2z(comdis_imock-rand2, **cosmo)
                else: 
                    collided_z = comdis2z(comdis_imock+rand2, **cosmo)

                appended_z.append(collided_z[0]) 
                print fibcoll_mock.z[i_mock], collided_z[0] 
                print fibcoll_mock.nbar[i_mock], nbarofz(collided_z[0])
                appended_nbar.append(nbarofz(collided_z[0]))

                if doublecheck: 
                    dlos_values.append(rand2) 

            else:                           # if not in the peak ----------------------------
                if correction['name'] in ('peak', 'peaknbar'): 
                    '''
                    generate random z based on redshift distribution of galaxies 
                    '''
                    # RA, Dec remain the same 
                    # weight = 1
                    appended_ra.append(fibcoll_mock.ra[i_mock])
                    appended_dec.append(fibcoll_mock.dec[i_mock])
                    appended_weight.append(1.0) 
                    
                    wtarg = np.random.random(1)*true_weight_max
                    zindx = np.floor(np.interp(wtarg, true_weight_cum, i_true)).astype(int)+1
                    
                    # fail safe
                    qqq = np.where(wtarg < true_weight_cum[0])[0]
                    zindx[qqq] = 0 

                    # assign redshift 
                    appended_z.append(true_z[zindx]) 

                elif correction['name'] in ('peakshot', 'bigfc_peakshot'): 
                    upweight_again.append(i_mock)

                else: 
                    raise NotImplementedError('asdfasdf') 
    
    if correction['name'] in ('peakshot', 'bigfc_peakshot'): 
        # re-upweighting for peak+shotnoise correction 
        for i_upweightagain in upweight_again: 
            fibcoll_mock.weight[i_upweightagain] += 1.

    print len(appended_ra), ' galaxies were peak corrected'
    # append artificial galaxies to catalog 
    fibcoll_mock.ra = np.concatenate([fibcoll_mock.ra, appended_ra])
    fibcoll_mock.dec = np.concatenate([fibcoll_mock.dec, appended_dec])
    fibcoll_mock.weight = np.concatenate([fibcoll_mock.weight, appended_weight])
    fibcoll_mock.z = np.concatenate([fibcoll_mock.z, appended_z])
        
    if catalog['name'].lower() in ('qpm', 'nseries'): 
        fibcoll_mock.comp = np.concatenate([fibcoll_mock.comp, appended_comp])
    elif 'cmass' in catalog['name'].lower():
        fibcoll_mock.comp = np.concatenate([fibcoll_mock.comp, appended_comp])
        fibcoll_mock.wsys = np.concatenate([fibcoll_mock.wsys, appended_wsys])
        fibcoll_mock.wnoz = np.concatenate([fibcoll_mock.wnoz, appended_wnoz])
        fibcoll_mock.nbar = np.concatenate([fibcoll_mock.nbar, appended_nbar])

    peakcorr_file = get_galaxy_data_file('data', **cat_corr) 
    
    if catalog['name'].lower() in ('qpm', 'nseries'): 
        np.savetxt(peakcorr_file, 
                np.c_[
                    fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, 
                    fibcoll_mock.weight, fibcoll_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 

    elif 'cmass' in catalog['name'].lower():          # CAMSS
        np.savetxt(peakcorr_file, 
                np.c_[
                    fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.nbar
                    fibcoll_mock.wsys, fibcoll_mock.wnoz, fibcoll_mock.weight, 
                    fibcoll_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 

    elif catalog['name'].lower() in ('tilingmock', 'lasdamasgeo', 'ldgdownnz'): 
        np.savetxt(peakcorr_file, 
                np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    else: 
        raise NotImplementedError('asdfasdf')

    if doublecheck: 
        np.savetxt(peakcorr_file+'.dlosvalues', np.c_[dlos_values], fmt=['%10.5f'], delimiter='\t') 

def build_photoz_peakcorrected_fibcol(doublecheck=False, **cat_corr): 
    ''' Build peak corrected fibercollided mock catalogs using assigned photometric redshifts
    that emulated photometric redshift errors (using cosmolopy) 

    Parameters
    ----------
    cat_corr : Catalog + Correction dictionary
    doublecheck : save dlos values to test it's working  

    Notes
    -----
    * Addresses a number of inefficiences in the previous peakshot code. 
        * No longer appends new peak corrected galaxies; instead the original galaxy's redshift is replaced with 
        peak corrected redshift. 
        * No longer downweights all upweighted galaxies only to re-upweight them later for collided galaxies that reside in 
        the tail fiber collided galaxies
    * Choise of +/- 3 sigma for peak sampled dLOS may need to be reviewed. 

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    if correction['name'].lower() != 'photozpeakshot': 
        raise NameError("Only accepts photozpeakshot correction")

    # fit functions (using lambda) ------------------------------------------------------------
    fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)
    
    # redshift limits 
    if catalog['name'].lower() in ('nseries'): 
        # set up mock catalogs 
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits
        n_mocks = 1 
    else: 
        raise NotImplementedError('Catalog not yet included')

    # read in fiber collided mocks with assigned photometric redshift  
    fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'photoz'}}

    data = galaxy_data('data', **fibcoll_cat_corr) 
    cosmo = data.cosmo      # survey comoslogy 
    
    # comoving distance of min and max redshifts    
    survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
    survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
    if 'weight' not in data.columns:
        # resolve nomenclature issue
        data.weight = data.wfc            

    f_peak = correction['fpeak']        # peak fraction 
    
    fcoll = np.where(data.weight == 0)     # fiber collided
    n_fcoll = len(fcoll[0])                 # Ngal fiber collided
    n_fc_peak = int(f_peak * np.float(n_fcoll))     # Ngal fc peak 
    n_fc_tail = n_fcoll - n_fc_peak                 # Ngal fc tail 
    
    # Comoving distance of upweighted galaxy
    Dc_upw = cosmos.distance.comoving_distance(
            data.zupw[fcoll], **cosmo) * cosmo['h']
    # Comoving distance of fibcollided galaxy photoz
    Dc_zphoto = cosmos.distance.comoving_distance(
            data.z_photo[fcoll], **cosmo) * cosmo['h']
   
    # photometric redshit dLOS 
    LOS_d_photo = Dc_zphoto - Dc_upw
    
    # first, determine the collided galaxies that are definitely in the tail based on 
    # photometric redshift 
    def_tail_dlos = np.where( (LOS_d_photo < -175.0) | (LOS_d_photo > 175.0) ) 
    def_tail = (fcoll[0])[def_tail_dlos]
    n_def_tail = len(def_tail) 
    print 'Ngal fiber collided', n_fcoll
    print 'Ngal fiber collided tail', n_fc_tail 
    print 'Ngal fiber collided definitely in tail', n_def_tail 
    
    # Then the rest of the collided galaxies are not definitely in the tail 
    not_def_tail = list(
            set(list(fcoll[0])) - set(list(def_tail))
            )
    print 'Ngal fiber collided not definitely in tail', len(not_def_tail)
    upw_def_tail = (data.upw_index)[def_tail]
    upw_def_not_tail = (data.upw_index)[not_def_tail]
    
    not_tail_fpeak = 1.0 - np.float(n_fc_tail - n_def_tail)/np.float(n_fcoll - n_def_tail)
    print 'fpeak of not definitely in tail', not_tail_fpeak 

    if doublecheck:     # check that the PDF of dLOS peak is properly generated
        dlos_values = [] 

    n_peakcorrected = 0     # Ngal peak corrected
    for i_mock in not_def_tail: 
        # go through each fibercollided galaxy not definitely in the tail        
        
        # use new fpeak that excludes definitely tail galaxies
        # to deterimine whether galxay is in peak or not 
        rand_num = np.random.random(1)      # random number

        if rand_num <= not_tail_fpeak:     # sampled in the peak 
            # sample a dLOS from the best-fit Gaussian dLOS PDF then place 
            # the collided galaxy with wfc = 1 dLOS away from the upweighted galaxy 
            # and then downweight the upweighted galaxy 

            data.weight[ (data.upw_index)[i_mock] ] -= 1.0  # downweight UPW galaxy 
            data.weight[i_mock] += 1.0  # upweight collided galaxy 
            if data.weight[i_mock] > 1.0:
                raise NameError('something went wrong') 

            # comoving distance of upweighted galaxy
            comdis_upw = cosmos.distance.comoving_distance(
                    data.z[ (data.upw_index)[i_mock] ], **cosmo) * cosmo['h']
        
            # compute the displacement within the peak using best-fit function 
            if correction['fit'].lower() in ('gauss'):  # Gaussian 

                rand1 = np.random.random(1) # random number
                
                # random dLOS +/- 3-sigma of the distribution 
                rand2 = np.random.random(1)
                rand2 = (-3.0 + 6.0 * rand2) * correction['sigma']
                
                peak_pofr = fit_func(rand2, correction['sigma']) # probability distribution
                
                while peak_pofr <= rand1: 
                    rand1 = np.random.random(1) # random number
                    
                    # random dLOS +/- 3-sigma of the distribution 
                    rand2 = np.random.random(1)
                    rand2 = (-3.0 + 6.0 * rand2) * correction['sigma']
                    
                    peak_pofr = fit_func(rand2, correction['sigma']) # probability distribution
                    
            else: 
                NotImplementedError('Not yet implemented') 
                
            # in case the displaced coliided galaxy falls out of bound (may generate large scale issues)
            # this will skew the dLOS displacement slightly at low and high redshift limits 
            if (comdis_upw + rand2 > survey_comdis_max) or (comdis_upw + rand2 < survey_comdis_min): 
                collided_z = comdis2z(comdis_upw-rand2, **cosmo)
            else:
                collided_z = comdis2z(comdis_upw+rand2, **cosmo)

            if doublecheck: 
                # append sample peak dLOS value to doublecheck 
                dlos_values.append(rand2) 

            #print data.z[ (data.upw_index)[i_mock] ], data.z[i_mock], collided_z[0] 
            data.z[i_mock] = collided_z[0]

            n_peakcorrected += 1
    
    print n_peakcorrected, ' galaxies were peak corrected'

    # write to file based on mock catalog  
    corrected_file = get_galaxy_data_file('data', **cat_corr) 
    
    if catalog['name'].lower() in ('nseries'): 
        np.savetxt(corrected_file, 
                np.c_[data.ra, data.dec, data.z, data.weight, data.comp, 
                        data.zupw, data.upw_index, data.z_photo], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%i', '%10.5f'], 
                delimiter='\t') 
    else: 
        raise NotImplementedError('asdfasdf')
    
    if doublecheck: 
        np.savetxt(corrected_file+'.dlosvalues', np.c_[dlos_values], fmt=['%10.5f'], delimiter='\t') 

def build_ldg_scratch(**cat_corr): 
    ''' Quick function to test fiber collision correction methods on LasDamasGeo mocks
        

    Parameters
    ----------
    cat_corr : catalog and correction dictionary

    Notes
    -----
    * scratch_peak_ang : Instead of using the

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
        
    omega_m = 0.25  # cosmology for LDG

    cosmo = {} 
    cosmo['omega_M_0'] = omega_m 
    cosmo['omega_lambda_0'] = 1.0 - omega_m 
    cosmo['h'] = 0.676
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 

    d_peak = 21.0
        
    # read rdzw file 
    data_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
    orig_file = ''.join([data_dir, 'sdssmock_gamma_lrgFull_zm_oriana', 
        str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.fibcoll.dat']) 

    orig_ra, orig_dec, orig_z, orig_wfc, z_upw, upw_index = np.loadtxt(orig_file, unpack=True, usecols=[0,1,2,3,4,5], 
            dtype={'names': ('ra', 'dec', 'z', 'wfc', 'zupw', 'upw_index'), 
                'formats': (np.float64, np.float64, np.float64, np.float64, np.float64, np.int32)})

    if correction['name'].lower() in ('scratch_peakknown'): 
        now_index = np.where(orig_wfc < 1)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc 
        
        peak_index = np.where(np.abs(dLOS) < d_peak)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        
        np.random.shuffle(dLOS[peak_index])     # shuffle the dLOS in the peak
        shuffle_dlos = dLOS[peak_index]
        
        print len(shuffle_dlos), np.float(len(shuffle_dlos))/np.float(len(orig_z[now_index]))

        for j, i_now_peak in enumerate(now_peak_index): 
            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + shuffle_dlos[j], **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0

    elif correction['name'].lower() in ('scratch_peakknown_ang'):

        now_index = np.where(orig_wfc < 1)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc 
        
        peak_index = np.where(np.abs(dLOS) < d_peak)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        np.random.shuffle(dLOS[peak_index])     # shuffle the dLOS in the peak
        shuffle_dlos = dLOS[peak_index]
        
        print len(shuffle_dlos), np.float(len(shuffle_dlos))/np.float(len(orig_z[now_index]))

        for j, i_now_peak in enumerate(now_peak_index): 

            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + shuffle_dlos[j], **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0

    elif correction['name'].lower() in ('scratch_peakknown_gauss'): 

        now_index = np.where(orig_wfc < 1)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc  # dLOS values for the pairs 

        peak_index = np.where(np.abs(dLOS) < d_peak)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now_peak in enumerate(now_peak_index): 

            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -1.0*d_peak + 2.0 * d_peak * rand1
            peakpofr = fit_func(rand2, 6.5) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -1.0*d_peak + 2.0 * d_peak * rand1
                peakpofr = fit_func(rand2, 6.5) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0
        
    elif correction['name'].lower() in ('scratch_allpeak_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        now_indices = (now_index[0])

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now in enumerate(now_indices): 

            orig_ra[i_now] = orig_ra[upw_index[i_now]]
            orig_dec[i_now] = orig_dec[upw_index[i_now]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -12.0 + 24.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -12.0 + 24.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now], **cosmo)*cosmo['h']
            
            if z_upw[i_now] < 0.: 
                print z_upw[i_now], comdis_upw, rand2, comdis_upw + rand2
                continue

            orig_z[i_now] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now] = 1.0
            orig_wfc[upw_index[i_now]] -= 1.0
    
    elif correction['name'].lower() in ('scratch_peakknown_gauss_divide'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc  # dLOS values for the pairs 

        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 
    
        append_ra, append_dec, append_z, append_wfc, append_wcomp = [], [], [], [], [] 
        for j, i_now_peak in enumerate(now_peak_index): 

            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            for jj in range(1,5): 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -15.0 + 30.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
                
                while peakpofr <= rand1: 
                    rand1 = np.random.random(1) 
                    rand2 = np.random.random(1) 

                    rand2 = -15.0 + 30.0 * rand1
                    peakpofr = fit_func(rand2, 3.8) 
                
                comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']
    
                if jj == 1: 
                    orig_z[i_now_peak] = comdis2z(comdis_upw + rand2, **cosmo)
                    orig_wfc[i_now_peak] = 0.25
                else: 
                    append_ra.append(orig_ra[upw_index[i_now_peak]]) 
                    append_dec.append(orig_dec[upw_index[i_now_peak]]) 
                    append_z.append(comdis2z(comdis_upw + rand2, **cosmo))
                    append_wfc.append(0.25)
                    append_wcomp.append(orig_wcomp[i_now_peak])
                    
            orig_wfc[upw_index[i_now_peak]] -= 1.0
        
        orig_ra = np.append( orig_ra, np.array(append_ra) ) 
        orig_dec = np.append( orig_dec, np.array(append_dec) ) 
        orig_z = np.append( orig_z, np.array(append_z) ) 
        orig_wfc = np.append( orig_wfc, np.array(append_wfc) ) 
        orig_wcomp = np.append( orig_wcomp, np.array(append_wcomp) ) 

    elif correction['name'].lower() in ('scratch_dividedw_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        now_indices = (now_index[0])

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now in enumerate(now_indices): 

            orig_ra[i_now] = orig_ra[upw_index[i_now]]
            orig_dec[i_now] = orig_dec[upw_index[i_now]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -12.0 + 24.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -12.0 + 24.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now], **cosmo)*cosmo['h']
            
            if z_upw[i_now] < 0.: 
                print z_upw[i_now], comdis_upw, rand2, comdis_upw + rand2
                continue

            orig_z[i_now] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now] = 0.7
            orig_wfc[upw_index[i_now]] -= 0.7
        
    # save to file 
    scratch_file = get_galaxy_data_file('data', **cat_corr) 
    np.savetxt(scratch_file, 
            np.c_[orig_ra, orig_dec, orig_z, orig_wfc ], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def build_ldg_nz_down(DorR, **cat_corr): 
    ''' Downsample LasDamas mocks in order to give the nbar(z) redshift dependence
    similar to CMASS sample 

    Parameters 
    ----------
    * cat_corr : catalog correction dictionary

    
    Notes
    -----
    * simultaneously generates nbar(z) file for new downsampled nbar(z) 
    * rather than using the actual CMASS sample, Nseries mock nbar(z) ratio is used 
    * since the main goal is to mostly to construct a redshif tdependence, 
    the details are spared
    * wcp are NOT assigned

    '''
    catalog = cat_corr['catalog']
    ldg_catalog = catalog.copy()
    ldg_catalog['name'] = 'lasdamasgeo'
    true_cat_corr = {'catalog': ldg_catalog, 'correction': {'name': 'true'}} # only true can go her 

    if catalog['name'].lower() != 'ldgdownnz': 
        raise NameError('kasdkfjlaskjdf') 
    
    ldg_nbar = 0.0000944233         # default LDG nbar(z) 
    
    # read in LDG file 
    ldg_file = get_galaxy_data_file(DorR, **true_cat_corr) 
    if DorR == 'data': 
        ra, dec, z, wcp = np.loadtxt(ldg_file, unpack=True, usecols=[0,1,2,3])
    elif DorR == 'random': 
        ra, dec, cz = np.loadtxt(ldg_file, unpack=True, usecols=[0,1,2])
        z =cz/299800.0
    ngal_tot = len(z) 

    # read in Nseries nbar(z) ratio file
    nseries_cat = {'name': 'nseries'} 
    nseries_cat_corr = {'catalog': nseries_cat, 'correction': {'name': 'true'}} 

    nbar_file = fc_nbar.get_nbar_file(**nseries_cat_corr) # avg nbar file
    ratio_nbar_file = ''.join([ '/'.join( nbar_file.split('/')[:-1] ), '/', 
        'ratio_', (nbar_file.split('/')[-1])]) 
    print ratio_nbar_file
    
    nseries_zmid, nseries_zlow, nseries_zhigh, nseries_fraction = np.loadtxt(
            ratio_nbar_file, unpack=True, usecols=[0,1,2,3]) 

    zmid = nseries_zmid 
    zlow = nseries_zlow 
    zhigh = nseries_zhigh

    nseries_zmid = nseries_zmid - 0.27  # shift nseries redshift down 
    nseries_zlow = nseries_zlow - 0.27  # shift nseries redshift down 
    nseries_zhigh = nseries_zhigh - 0.27  # shift nseries redshift down 
    
    remove_indices = [] 
    ngal_remove = 0
    if DorR == 'random':
        new_nz = [] 

    # loop through redshift bins and downsample 
    for i_z, n_zmid in enumerate(nseries_zmid): 

        z_bin = np.where( (z >= nseries_zlow[i_z]) & (z < nseries_zhigh[i_z]) ) 

        ngal_bin = len(z[z_bin]) # Ngal for redshift bin 
        
        bin_down_frac = nseries_fraction[i_z] 
    
        if DorR == 'random':    # make downsampled nbar(z) file
            if n_zmid > 0.0: 
                new_nz.append(ldg_nbar * bin_down_frac) 
            
        # number of galaxies to remove 
        ngal_bin_remove = np.float(ngal_bin) * (1.0 - bin_down_frac) 
    
        if ngal_bin_remove > 0:  
            # indices to remove
            ngal_remove += int(np.rint(ngal_bin_remove))
            remove_indices += random.sample( z_bin[0], int(np.rint(ngal_bin_remove)) )
    
    if DorR == 'random': 
        new_nz += [0.0 for i in range(len(zmid) - len(new_nz))]

        new_nz_file = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/', 
            'nbar-lasdamasgeo.down_nz.dat']) 
        print new_nz_file
        np.savetxt(new_nz_file,
                np.c_[zmid, zlow, zhigh, new_nz],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e'], delimiter='\t')

    print ngal_remove, ' Galaxies will be removed from mock realization'
    print 100.0 * np.float(ngal_remove)/np.float(ngal_tot), '% of galaxies will be removed'
    
    remain_indices = list(set(range(ngal_tot)) - set(remove_indices)) # remove 'remove indices'

    if len(remain_indices) + len(remove_indices) == ngal_tot: 
        pass
    else: 
        raise TypeError('Remove indices not working') 
    
    ldg_down_nz_file = get_galaxy_data_file(DorR, **cat_corr) 
    if DorR == 'data': 
        np.savetxt(ldg_down_nz_file, 
                np.c_[
                    ra[remain_indices], dec[remain_indices], 
                    z[remain_indices], wcp[remain_indices]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    elif DorR == 'random': 
        np.savetxt(ldg_down_nz_file, 
                np.c_[
                    ra[remain_indices], dec[remain_indices], 
                    z[remain_indices]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    return 

def build_ldg_bigfc(**cat_corr): 
    ''' Build LasDamas realizations with bigger fiber collision angular scale
    to match fiber collided fraction

    Parameters 
    ----------
    * cat_corr : catalog correction dictionary

    Notes
    -----
    * Calls IDL code ldg_bigfibcoll_wcp_assign.pro 
    * theta_ang : 

    '''
    catalog = cat_corr['catalog']
    
    # IDL command 
    fibcollided_cmd = 'idl -e "ldg_bigfibcoll_wcp_assign,'+str(catalog['n_mock'])+", '"+\
            str(catalog['letter'])+"'"+'"'
    print fibcollided_cmd
    os.system(fibcollided_cmd)  # call IDL code 

    return fibcollided_cmd

def build_ldgdownnz_bigfc(**cat_corr): 
    ''' Build LasDamas realizations with bigger fiber collision angular scale
    to match fiber collided fraction

    Parameters 
    ----------
    * cat_corr : catalog correction dictionary

    Notes
    -----
    * Calls IDL code ldg_bigfibcoll_wcp_assign.pro 
    * theta_ang : 

    '''
    catalog = cat_corr['catalog']
    
    # IDL command 
    fibcollided_cmd = 'idl -e "ldgdownnz_bigfc_wcp_assign,'+str(catalog['n_mock'])+", '"+\
            str(catalog['letter'])+"'"+'"'
    print fibcollided_cmd
    os.system(fibcollided_cmd)  # call IDL code 

    return fibcollided_cmd

def build_nseries_scratch(**cat_corr): 
    ''' Quick function to test fiber collision correction methods on Nseries mocks
        

    Parameters
    ----------
    cat_corr : catalog and correction dictionary

    Notes
    -----
    * scratch_peak_ang : Instead of using the

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
        
    omega_m = 0.31

    cosmo = {} 
    cosmo['omega_M_0'] = omega_m 
    cosmo['omega_lambda_0'] = 1.0 - omega_m 
    cosmo['h'] = 0.676
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        
    # read rdzw file 
    data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
    orig_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.rdzwc']) 
    orig_ra, orig_dec, orig_z, orig_wfc, z_upw, upw_index = np.loadtxt(orig_file, unpack=True, usecols=[0,1,2,4,5,6], 
            dtype={'names': ('ra', 'dec', 'z', 'wfc', 'zupw', 'upw_index'), 
                'formats': (np.float64, np.float64, np.float64, np.float64, np.float64, np.int32)})

    # file with completeness
    mask_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.mask_info']) 
    orig_wcomp = np.loadtxt(mask_file, unpack=True, usecols=[0]) 

    if correction['name'].lower() in ('scratch_peakknown'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc 
        
        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        
        np.random.shuffle(dLOS[peak_index])     # shuffle the dLOS in the peak
        shuffle_dlos = dLOS[peak_index]
        
        #print len(shuffle_dlos), np.float(len(shuffle_dlos))/np.float(len(orig_z[now_index]))

        for j, i_now_peak in enumerate(now_peak_index): 
            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + shuffle_dlos[j], **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0

        scratch_file = get_galaxy_data_file('data', **cat_corr) 
        np.savetxt(scratch_file, 
                np.c_[orig_ra, orig_dec, orig_z, orig_wfc, orig_wcomp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    elif correction['name'].lower() in ('scratch_peakknown_ang'):

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc 
        
        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        
        np.random.shuffle(dLOS[peak_index])     # shuffle the dLOS in the peak
        shuffle_dlos = dLOS[peak_index]
        
        #print len(shuffle_dlos), np.float(len(shuffle_dlos))/np.float(len(orig_z[now_index]))

        for j, i_now_peak in enumerate(now_peak_index): 

            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + shuffle_dlos[j], **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0

    elif correction['name'].lower() in ('scratch_peakknown_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc  # dLOS values for the pairs 

        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now_peak in enumerate(now_peak_index): 

            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -15.0 + 30.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -15.0 + 30.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0
        
    elif correction['name'].lower() in ('scratch_allpeak_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        now_indices = (now_index[0])

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now in enumerate(now_indices): 

            orig_ra[i_now] = orig_ra[upw_index[i_now]]
            orig_dec[i_now] = orig_dec[upw_index[i_now]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -12.0 + 24.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -12.0 + 24.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now], **cosmo)*cosmo['h']
            
            if z_upw[i_now] < 0.: 
                print z_upw[i_now], comdis_upw, rand2, comdis_upw + rand2
                continue

            orig_z[i_now] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now] = 1.0
            orig_wfc[upw_index[i_now]] -= 1.0
    
    elif correction['name'].lower() in ('scratch_peakknown_gauss_divide'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc  # dLOS values for the pairs 

        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 
    
        append_ra, append_dec, append_z, append_wfc, append_wcomp = [], [], [], [], [] 
        for j, i_now_peak in enumerate(now_peak_index): 

            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            for jj in range(1,5): 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -15.0 + 30.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
                
                while peakpofr <= rand1: 
                    rand1 = np.random.random(1) 
                    rand2 = np.random.random(1) 

                    rand2 = -15.0 + 30.0 * rand1
                    peakpofr = fit_func(rand2, 3.8) 
                
                comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']
    
                if jj == 1: 
                    orig_z[i_now_peak] = comdis2z(comdis_upw + rand2, **cosmo)
                    orig_wfc[i_now_peak] = 0.25
                else: 
                    append_ra.append(orig_ra[upw_index[i_now_peak]]) 
                    append_dec.append(orig_dec[upw_index[i_now_peak]]) 
                    append_z.append(comdis2z(comdis_upw + rand2, **cosmo))
                    append_wfc.append(0.25)
                    append_wcomp.append(orig_wcomp[i_now_peak])
                    
            orig_wfc[upw_index[i_now_peak]] -= 1.0
        
        orig_ra = np.append( orig_ra, np.array(append_ra) ) 
        orig_dec = np.append( orig_dec, np.array(append_dec) ) 
        orig_z = np.append( orig_z, np.array(append_z) ) 
        orig_wfc = np.append( orig_wfc, np.array(append_wfc) ) 
        orig_wcomp = np.append( orig_wcomp, np.array(append_wcomp) ) 

    elif correction['name'].lower() in ('scratch_dividedw_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        now_indices = (now_index[0])

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now in enumerate(now_indices): 

            orig_ra[i_now] = orig_ra[upw_index[i_now]]
            orig_dec[i_now] = orig_dec[upw_index[i_now]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -12.0 + 24.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -12.0 + 24.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now], **cosmo)*cosmo['h']
            
            if z_upw[i_now] < 0.: 
                print z_upw[i_now], comdis_upw, rand2, comdis_upw + rand2
                continue

            orig_z[i_now] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now] = 0.7
            orig_wfc[upw_index[i_now]] -= 0.7
        
    # save to file 
    scratch_file = get_galaxy_data_file('data', **cat_corr) 
    np.savetxt(scratch_file, 
            np.c_[orig_ra, orig_dec, orig_z, orig_wfc, orig_wcomp], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

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
    ''' Given comoving distance and cosmology, determine z 
    Comoving distance *has* to be in Mpc/h
    '''
    z_arr = np.array([0.0+0.05*np.float(i) for i in range(21)]) 
    dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo)*cosmo['h']

    dmz_spline = sp.interpolate.interp1d(dm_arr, z_arr, kind='cubic') 
    
    z = dmz_spline(comdis)
    return z 
