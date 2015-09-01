'''

Correction galaxy/random catalog 

'''

import true as corr_true
import fibcollided as corr_fibcol

from util import catalog as cata 

class correction: 

    def __init__(self, cat_corr): 
        """ A class describing the correction to galaxy catalog of simulations or BOSS data
        """
    
        # correction dictionary 
        self.corr_dict = {
                'true': corr_true, 
                'upweight': corr_fibcol 
                }

        self.cat_corr = cat_corr
    
        corr = cat_corr['correction']
        if corr['name'] not in (self.corr_dict).keys(): 
            raise ValueError()
        else: 
            self.corr_mod = self.corr_dict[corr['name']]

    def build(self): 
        """ Build corrected galaxy catalog
        """

        # run correction 
        (self.corr_mod).build(self.cat_corr)
    
        return None 

    def file(self):
        """ Name of corrected galaxy catalog file 
        """
        cat = cata.catalog(self.cat_corr)
        file_list = cat.file()      # list of file name elements
    
        corr_str = (self.corr_mod).file(self.cat_corr)
        file_list.insert( -1, corr_str )

        return ''.join(file_list)

if __name__=="__main__":
    cat_corr = {'catalog': {'name': 'nseries', 'n_mock': 1}, 'correction': {'name': 'upweight'}} 
    corr_class = correction(cat_corr)   
    print corr_class.file()
    print corr_class.build()

"""
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
                    elif correction['name'].lower() == 'photozenvpeakshot': 
                        # Peak Shot correction using photometric redshift information
                        photozenv_corr.build_photoz_env_dlospeak_fibcol(cat_corr, doublecheck=True)
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
                file_data = np.loadtxt(file_name, skiprows=1, unpack=True, usecols=[0,1,2,3,4,5,6,7])   

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
                file_data = np.loadtxt(file_name, skiprows=1, unpack=True, usecols=[0,1,2,3,4])  

                # assign to data columns class
                for i_col, catalog_column in enumerate(catalog_columns): 
                    setattr(self, catalog_column, file_data[i_col]) 
    
        else: 
            raise NameError('not yet coded') 




        if catalog['name'].lower() == 'nseries': # N-series 

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
                elif correction['name'].lower() == 'photozpeakshot':    
                    # peak shot correction utilizing photoz
                    
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
                elif correction['name'].lower() == 'photozenvpeakshot':     
                    # Photometric Redshift + Galaxy Environment + Peak + Shot Noise Correction 
                    
                    if correction['fit'].lower() in ('gauss'): 
                        # specify best fit (expon or gauss) function to peak 
                        # correction specifier string 
                        corr_str = ''.join(['.photoz.', 'env_d', str(correction['n_NN']), 'NN.',
                            correction['fit'].lower(), '.peakshot.sigma', str(correction['sigma']), 
                            '.fpeak', str(correction['fpeak'])]) 
                    else: 
                        raise NotImplementedError("Only Gaussian Best-Fit coded") 

                    file_name = ''.join([data_dir, 
                        'CutskyN', str(catalog['n_mock']), '.fibcoll', corr_str, cosmo_str, '.dat' ]) 
                else: 
                    raise NameError('not yet coded') 

            elif DorR == 'random': 
                # random catalog 
                file_name = ''.join([data_dir, 'Nseries_cutsky_randoms_50x_redshifts_comp.dat']) 
        

        catalog = (self.cat_corr)['catalog']        # catalog 
        correction = (self.cat_corr)['correction']      # correction 
        
        if catalog['name'].lower() == 'lasdamasgeo':                # LasDamasGeo --------------

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
                data_dir = ''              

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
            if catalog['name'].lower != 'cmass':
                data_dir += 'dr12v5/'
            
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
                elif 'cmasslowz' in catalog['name'].lower(): 
                    # CMASS LOWZ DR12v5 combined sample (high redshift bin) 

                    cmasslowz_str = ''  # to account for the three different types of combined sample
                    if 'e2' in catalog['name'].lower(): 
                        cmasslowz_str = 'e2'
                    elif 'e3' in catalog['name'].lower(): 
                        cmasslowz_str = 'e3'
                    elif 'tot' in catalog['name'].lower(): 
                        cmasslowz_str = 'tot'
                    
                    if '_low' in catalog['name'].lower(): 
                        zbin_str = '-low' 
                    elif '_high' in catalog['name'].lower(): 
                        zbin_str = '-high'
                    else: 
                        raise NameError("Must specify redshift bin of CMASS LOWZ sample") 

                    # CMASS LOWZ combined sample high redshift bin DR12v5 
                    if correction['name'].lower() in ('upweight'):      
                        # upweighted
                        file_name = ''.join([data_dir, 
                            'cmasslowz', cmasslowz_str, '-dr12v5-N', zbin_str, '.dat']) # hardcoded
                    elif correction['name'].lower() in ('peakshot'):    
                        # peakshot
                        # correction string in file name 
                        corr_str = ''.join([
                            '.', correction['fit'].lower(), '.', correction['name'].lower(), 
                            '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                        file_name = ''.join([data_dir, 
                            'cmasslowz', cmasslowz_str, '-dr12v5-N', zbin_str, '.dat', corr_str]) 

            elif DorR == 'random':                  # random catalog 
                if catalog['name'].lower() == 'cmass': 
                    file_name = ''.join([data_dir, 
                        'cmass-dr12v4-N-Reid.ran.dat'])
                elif 'cmasslowz' in catalog['name'].lower(): 
                    cmasslowz_str = ''
                    if 'e2' in catalog['name'].lower(): 
                        cmasslowz_str = 'e2'
                    elif 'e3' in catalog['name'].lower(): 
                        cmasslowz_str = 'e3'
                    elif 'tot' in catalog['name'].lower(): 
                        raise NotImplementedError("Not Implemented Error")
                    
                    if '_low' in catalog['name'].lower(): 
                        zbin_str = '-low' 
                    elif 'high' in catalog['name'].lower(): 
                        zbin_str = '-high'
                    else: 
                        raise NameError("Must specify redshift bin of CMASS LOWZ sample") 

                    file_name = ''.join([data_dir, 
                        'cmasslowz', cmasslowz_str, '-dr12v5-N', zbin_str, '.ran.dat'])
                else: 
                    raise NotImplementedError('lskdfjaklsdfj')

        return file_name 
"""
