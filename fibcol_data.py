import numpy as np
import time 
import os.path
import subprocess
import cosmolopy as cosmos
import warnings 
import matplotlib.pyplot as plt

# Classes ------------------------------------------------------------
class galaxy_data: 
    def __init__(self, DorR, readdata=True, **cat_corr): 
        '''
        Read in the data or random file and store all the appropriate values  
        '''
        self.cat_corr = cat_corr
        catalog = cat_corr['catalog'] 
        correction = cat_corr['correction'] 

        # LasDamas Geo ---------------------------------------------------------------------------------------------------------
        if catalog['name'].lower() == 'lasdamasgeo': 
            omega_m = 0.25
            # Data ---------------------------------------------------------------------------
            if DorR == 'data': 
                catalog_columns = ['ra', 'dec', 'z', 'weight']          # columns that this catalog data will have  
                self.columns = catalog_columns

                # True LasDamasGeo (No FiberCollisions) ---------------------------------------------------------------------- 
                if correction['name'].lower() == 'true': 
                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog['n_mock'])+catalog['letter']+'_no.rdcz.dat'

                    if readdata == True: 
                        if os.path.isfile(file_name) == False: 
                            # if file does not exists, make file  
                            build_true(**cat_corr)

                        # read data 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights (1.0)

                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 
                    '''
                    #file_name = '/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/sdssmock_gamma_lrgFull_zm_oriana'+\
                    #        str("%02d" % catalog['n_mock'])+catalog['letter']+'_no.rdcz.dat'

                    if readdata == True: 
                        # Read data 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])         # ra, dec, CZ

                        for i_col, catalog_column in enumerate(catalog_columns):            # assign to class
                            if catalog_column == 'z': 
                                column_data = file_data[i_col]/299800.0     # divide by speed of light
                            elif catalog_column == 'weight':                    
                                column_data = np.array([1.0 for j in range(len(file_data[0,:]))]) # no weights for true (all 1) 
                            else: 
                                column_data = file_data[i_col]
                            # assign to class
                            setattr(self, catalog_column, column_data) 
                    '''

                # Fibercollisions Corrected by Upweight correction ------------------------------------------------------------ 
                elif correction['name'].lower() in ('upweight', 'fibcol', 'shotnoise', 'hectorsn', 'floriansn'): 
                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog['n_mock'])+catalog['letter']+'_no.rdcz.fibcoll.dat'

                    if readdata == True: 
                        if os.path.isfile(file_name) == False: 
                            # if file does not exists, make file  
                            build_fibercollided(**cat_corr)

                        # read data 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

                # Fibercollisions Corrected by peak correction ---------------------------------------------------------------------
                elif correction['name'].lower() in ('peak', 'peaknbar', 'peakshot', 'peaktest'): 
                    if correction['name'].lower()  == 'peak': 
                        # to correct for poor naming convention 
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
                    
                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog['n_mock'])+catalog['letter']+'_no.rdcz.fibcoll.dat'+corr_str
                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_peakcorrected_fibcol(sanitycheck=True, **cat_corr)
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

                elif correction['name'].lower() in ('vlospeaknbar', 'vlospeakshot'): 
                    # specify peak correction fit (expon or gauss) 
                    if correction['fit'].lower() in ('gauss', 'expon'): 
                        fit_str = correction['fit'].lower()
                        corr_str = '.'+fit_str+'.'+correction['name'].lower()+\
                                '.sigma'+str(correction['sigma'])+'.fpeak'+str(correction['fpeak'])
                    elif correction['fit'].lower() == 'real': 
                        fit_str = 'real'
                        corr_str = '.'+fit_str+'.'+correction['name'].lower()+'.fpeak'+str(correction['fpeak'])
                    else: 
                        raise NameError('peak fit has to be specified: gauss or expon') 

                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog['n_mock'])+catalog['letter']+'_no.rdcz.fibcoll.dat'+corr_str
                    
                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_vlospeak_fibcol(**cat_corr)
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

                # Fibercollisions Corrected by all peak correction ---------------------------------------------------------------------
                elif correction['name'].lower() in ('allpeak', 'allpeakshot'):

                    # specify peak correction fit (expon or gauss) 
                    if correction['fit'].lower() == 'gauss': 
                        fit_str = 'gauss' 
                    elif correction['fit'].lower() == 'expon': 
                        fit_str = 'expon'
                    else: 
                        raise NameError('peak fit has to be specified: gauss or expon') 

                    corr_str = '.'+fit_str+'.'+correction['name'].lower()+\
                            '.sigma'+str(correction['sigma'])

                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog['n_mock'])+catalog['letter']+'_no.rdcz.fibcoll.dat'+corr_str
                    
                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_peakcorrected_fibcol(sanitycheck=True, **cat_corr)
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

                # Test adjustments ------------------------------------------------------------------------------------------------ 
                elif (correction['name'].lower() == 'randrm'):  
                    # get file name 
                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog['n_mock'])+catalog['letter']+'_no.rdcz.fibcoll.dat.'+correction['name'].lower() 
                    
                    # read data
                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_test_adjustments(sanitycheck=True, **cat_corr)
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign data columns to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

                else: 
                    raise NameError('Correction Name Unknown') 

            # Random Catalogs ------------------------------------------------------------------------------------------------------
            elif DorR.lower() == 'random': 
                catalog_columns = ['ra', 'dec', 'z']        # columns that this catalog random will have (NOTE CZ IS CONVERTED TO Z) 
                self.columns = catalog_columns
                
                # True ------------------------------------------------------------------------------------------------------
                if correction['name'].lower() in ('true', 'upweight', 'peaknbar', 'peakshot', 'shotnoise', 'vlospeakshot', 'hectorsn', 'floriansn'): 
                    # uncorrected LasDamasGeo random catalogue
                    # Upweight and peakshot have been included because they only change nbar(z) by < 1% 
                    
                    # hardcoded
                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat'

                    # read data
                    if readdata == True: 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])

                        for i_col, catalog_column in enumerate(catalog_columns): 
                            if catalog_column == 'z': 
                                column_data = file_data[i_col]/299800.0
                            else: 
                                column_data = file_data[i_col]
                            # assign to class
                            setattr(self, catalog_column, column_data)

                # Corrections ------------------------------------------------------------------------------------------------
                else:  
                    pass
                    '''
                    # corrections for random catalog
                    #if correction['name'].lower() == 'upweight': 
                    #    # upweight correction 
                    #    corr_str = '.upweight'

                    if correction['name'].lower() in ('peak', 'peaknbar', 'peaktest'): 
                        # Correction methods that have both peak and tail contributions   
                        # peak/peaknbar = peak + tail correction 
                        # peaktest = fpeak peak correction + remove rest 
                        # peakshot = fpeak peak correction + shot noise for rest

                        # specify peak correction fit (expon or gauss) 
                        if correction['fit'].lower() == 'gauss': 
                            fit_str = 'gauss' 
                        elif correction['fit'].lower() == 'expon': 
                            fit_str = 'expon'
                        else: 
                            raise NameError('peak fit has to be specified: gauss or expon') 
                        
                        if correction['name'].lower() == 'peak': 
                            # just for peak+nbar correction, for consistency 
                            correction['name'] = 'peaknbar' 

                        corr_str = ''.join(['.', fit_str, '.', correction['name'].lower(), 
                            '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])]) 

                    elif correction['name'].lower() in ('allpeak', 'allpeakshot'):
                        # all in the peak
                        # allpeak = all fc galaxies in peak with weight fpeak
                        # allpeakshot = all fc galaxies in peak with weight fpeak, (1-fpeak) weights are then shot noise corrected

                        # specify peak correction fit (expon or gauss) 
                        if correction['fit'].lower() == 'gauss': 
                            fit_str = 'gauss' 
                        elif correction['fit'].lower() == 'expon': 
                            fit_str = 'expon'
                        else: 
                            raise NameError('peak fit has to be specified: gauss or expon') 

                        corr_str = '.'+fit_str+'.'+correction['name'].lower()+'.sigma'+str(correction['sigma'])

                    # Test adjustments --------------------------------------------------------------------------------------------
                    elif (correction['name'].lower() == 'randrm'): 
                        corr_str = '.'+correction['name'].lower()

                    else: 
                        raise NameError("specify correction method") 
                    
                    # file name 
                    file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/', 
                        'sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat.allmocks', corr_str]) 
                    
                    if readdata == True: 
                        # if read data is true, then read data 

                        if os.path.isfile(file_name) == False:                            
                            # if corrected random file does *not* exist
                            print 'Constructing ', file_name
                            build_corrected_randoms(sanitycheck=False, **cat_corr) 

                        # read corrected random file
                        print 'Reading ', file_name                     # just because it takes forever
                        t0 = time.time() 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2]) 
                        print 'took ', (time.time()-t0)/60.0, ' mins'

                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            # assign to class
                            setattr(self, catalog_column, column_data)
                    '''

        # TILING MOCKS --------------------------------------------------------------------------------------------------- 
        elif catalog['name'].lower() == 'tilingmock':
            omega_m = 0.274
            # MOCK --------------------------------------------------------------------------------------------------- 
            if DorR == 'data': 
                # columns that this catalog data will have  
                catalog_columns = ['ra', 'dec', 'z', 'weight']
                self.columns = catalog_columns

                data_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'      # data directory

                # True Original Tiling mocks ---------------------------------------------------------------------------
                if correction['name'].lower() == 'true': 
                    # all weights = 1 (fibercollisions *not* imposed) 
                    file_name = data_dir+'cmass-boss5003sector-icoll012.zlim.dat'
                        
                    if readdata == True: 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_true(**cat_corr) 

                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weight

                        # assign to data columns class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            # assign to class
                            setattr(self, catalog_column, column_data) 

                # Fibercollision Upweight Correction --------------------------------------------------------------------  
                elif correction['name'].lower() in ('upweight', 'fibcol', 'shotnoise', 'floriansn', 'hectorsn'): 
                    file_name = data_dir+'cmass-boss5003sector-icoll012.zlim.fibcoll.dat'

                    if readdata == True: 
                        # if fibercollided file does not exist 
                        if os.path.isfile(file_name) == False: 
                            print 'Contructing ', file_name 
                            build_fibercollided(**cat_corr)                 # build-in fibercollisions using spherematch idl code

                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weight

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            # assign to data column class
                            setattr(self, catalog_column, column_data) 
    
                # Fibercollision peak+tail correction ---------------------------------------------------------------------- 
                elif correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot'): 
                    # Correction methods that have both peak and tail contributions   
                    # peak/peaknbar = peak + tail correction 
                    # peaktest = fpeak peak correction + remove rest 
                    # peakshot = fpeak peak correction + shot noise for rest

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
                    
                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_peakcorrected_fibcol(sanitycheck=True, **cat_corr)
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to data columns class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 
                
                # Fibercollision VLOS peak+tail correction ---------------------------------------------------------------------- 
                elif correction['name'].lower() in ('vlospeaknbar', 'vlospeakshot'): 
                    # specify fit (expon or gauss) to peak 
                    if correction['fit'].lower() in ('gauss', 'expon'): 
                       pass 
                    else: 
                        raise NameError('peak fit has to be specified: gauss or expon') 

                    # correction string in file name 
                    corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])

                    file_name = data_dir+'cmass-boss5003sector-icoll012.zlim.dat'+corr_str 
                    
                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_vlospeak_fibcol(**cat_corr)
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to data columns class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

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
                    
                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_peakcorrected_fibcol(sanitycheck=True, **cat_corr)
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

                else: 
                    raise NameError('Correction Name Unknown') 

            elif DorR.lower() == 'random': 
                # For Randoms ------------------------------------------------------------------------------------------------------
                # columns that this catalog random will have
                catalog_columns = ['ra', 'dec', 'z'] 
                self.columns = catalog_columns

                random_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'
                
                # True Original Random Tiling Mock ----------------------------------------------------------------------  
                if correction['name'].lower() in ('true', 'upweight', 'peaknbar', 'peakshot', 'vlospeakshot', 'shotnoise', 'floriansn', 'hectorsn'): 
                    # uncorrected random catalog for Tiling Mocks 
                    file_name = random_dir + 'randoms-boss5003-icoll012-vetoed.zlim.dat'

                    if readdata == True: 
                        if os.path.isfile(file_name) == False:                            
                            print 'Constructing ', file_name
                            build_corrected_randoms(sanitycheck=False, **cat_corr)       # impose redshift limit

                        # read random file
                        print 'Reading ', file_name                     # just because it takes forever
                        t0 = time.time() 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])
                        print 'took ', (time.time()-t0)/60.0, ' mins'

                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            # assign data column to class
                            setattr(self, catalog_column, column_data)

                # Corrected Random Tiling Mock ------------------------------------------------------------------------
                else: 
                    print "WHat are you doing here?"
                    pass
                    # corrected string for file name 
                    #if correction['name'].lower() == 'upweight': 
                    #    # random corrected for fiber collision upweight correction 
                    #    corr_str = '.upweight'
                    '''
                    if correction['name'].lower() in ('peak', 'peaknbar', 'peaktest'): 
                        # random corrected for fiber collision peak correction 
                
                        if correction['name'].lower() == 'peak': 
                            # correct for poor naming convention 
                            correction['name'] = 'peaknbar'

                        # peak fit method 
                        if correction['fit'].lower() in ('gauss', 'expon'): 
                            pass
                        else: 
                            raise NameError('peak fitting has to be gaussian or exponential') 

                        corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                            '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])]) 

                    elif correction['name'].lower() in ('allpeak', 'allpeakshot'):
                        # random corrected for fiber collision peak correction 

                        # check peak fit method 
                        if correction['fit'].lower() in ('gauss', 'expon'): 
                            pass
                        else: 
                            raise NameError('peak fitting has to be gaussian or exponential') 

                        corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                            '.sigma', str(correction['sigma'])])
                    else: 
                        raise NameError('correction method not coded') 

                    file_name = '/mount/riachuelo1/hahn/data/tiling_mocks/'+\
                            'randoms-boss5003-icoll012-vetoed.zlim.dat'+corr_str

                    if readdata == True: 
                        # if corrected random file does *not* exist
                        if os.path.isfile(file_name) == False:                            
                            print 'Constructing ', file_name
                            build_corrected_randoms(sanitycheck=False, **cat_corr) 

                        # read random file
                        print 'Reading ', file_name                     # just because it takes forever
                        t0 = time.time() 
                        # read corrected random file
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2]) 
                        print 'took ', (time.time()-t0)/60.0, ' mins'

                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            # assign to class
                            setattr(self, catalog_column, column_data)
                    '''
                
        # QPM ------------------------------------------------------------------------------------------------------------ 
        elif catalog['name'].lower() == 'qpm':  
            omega_m = 0.31
            # Mock ------------------------------------------------------------------------------------------------------
            if DorR == 'data': 
                # columns that this catalog data will have  
                catalog_columns = ['ra', 'dec', 'z', 'wfkp', 'wfc', 'comp']
                self.columns = catalog_columns

                data_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'              # data directory

                # True  -----------------------------------------------------------------------------------------------
                if correction['name'].lower() == 'true': 
                    # all weights = 1 (fibercollisions *not* imposed) 
                    file_name = ''.join([data_dir, 'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.dat']) 
                    if readdata == True: 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_true(**cat_corr) 

                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5])         # ra, dec, z, wfkp, wfc, comp

                        # assign to data columns class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            setattr(self, catalog_column, file_data[i_col]) 
                
                # Upweight ------------------------------------------------------------------------------------------------
                elif correction['name'].lower() in ('upweight', 'shotnoise', 'floriansn', 'hectorsn'): 
                    
                    file_name = ''.join([data_dir, 'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.fibcoll.dat']) 
                    
                    if readdata == True: 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_fibercollided(**cat_corr) 
                        
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5])         # ra, dec, z, wfkp, wfc, comp

                        # assign to object data columns 
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            setattr(self, catalog_column, file_data[i_col]) 

                # peak correction methods ------------------------------------------------------------------------------------------
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

                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_peakcorrected_fibcol(**cat_corr)  # build peak corrected file 
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5])         # ra, dec, z, wfkp, wfc, comp

                        # assign to data columns class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

                elif correction['name'].lower() in ('tailupw'):     
                    # only upweight uncorrelated chance alignment fc pairs 
                    corr_str = 'tailupw' 
                    file_name = ''.join([data_dir, 
                        'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.fibcoll.', corr_str, '.dat'
                        ]) 

                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_tailupweight_fibcol(**cat_corr)  # build peak corrected file 
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5])   # ra, dec, z, wfkp, wfc, comp
                        # assign to data columns class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

                elif correction['name'].lower() in ('vlospeaknbar', 'vlospeakshot'): 

                    if correction['fit'].lower() in ('gauss', 'expon'): 
                        # correction specifier string 
                        corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), 
                            '.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])]) 
                    elif correction['fit'].lower() in ('real'): 
                        corr_str = ''.join(['.', correction['fit'].lower(), '.', correction['name'].lower(), '.fpeak', str(correction['fpeak'])]) 
                    else: 
                        raise NameError('peak fit has to be specified: gauss or expon') 


                    file_name = ''.join([data_dir, 'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.fibcoll', corr_str, '.dat']) 

                    if readdata == True: 
                        # if file doesn't exist create file 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_vlospeak_fibcol(**cat_corr)
                        
                        # read file 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5])         # ra, dec, z, wfkp, weights, comp

                        # assign to data columns class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 
                else: 
                    raise NameError('not yet coded') 
                
            # Random ----------------------------------------------------------------------------------------------------------------
            elif DorR == 'random': 
                # data columns
                catalog_columns = ['ra', 'dec', 'z', 'wfkp']

                data_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/' 
                
                #if correction['name'].lower() in ('true'): 
                file_name = ''.join([data_dir, 'a0.6452_rand50x.dr12d_cmass_ngc.vetoed.dat'])             # hardcoded to 50x so it does'nt take forever

                if readdata == True: 
                    if os.path.isfile(file_name) == False: 
                        print 'Constructing ', file_name 
                        build_qpm_true_random() 

                    file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])             # ra, dec, z, wfkp

                    # assign to object data columns
                    for i_col, catalog_column in enumerate(catalog_columns): 
                        setattr(self, catalog_column, file_data[i_col])
                '''
                else: 
                    if correction['name'].lower() in ('upweight', 'shotnoise'): 
                        # upweight
                        corr_str = '.upweight'

                    elif correction['name'].lower() == 'peakshot': 
                        # peakshot
                        corr_str = ''.join(['.peakshot.sigma', str(correction['sigma']), '.fpeak', str(correction['fpeak'])])
                    else: 
                        raise NameError('Not Yet Codeded') 

                    file_name = ''.join([data_dir, 'a0.6452_rand50x.dr12d_cmass_ngc', corr_str, '.dat'])
                    
                    if readdata == True: 
                        if os.path.isfile(file_name) == False: 
                            print 'Constructing ', file_name 
                            build_corrected_randoms(sanitycheck=False, **cat_corr) 

                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4])             # ra, dec, z, nbar, wveto

                        # assign to object data columns
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            setattr(self, catalog_column, file_data[i_col])
                '''
        
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

        else: 
            raise NameError('not yet coded') 

        self.file_name = file_name          # save file name 

        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.7 
        cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        self.cosmo = cosmo 


# ------------------------------------------------------------------------
# Build functions

def build_true(**cat_corr): 
    '''
    Adjust original data for convenience purposes: 
    For Tiling Mock: 
        Reads in true mock catalogs and imposes predetermined redshift limits on them and attaches a flag 
        mostly hardcoded since it's a simple procedure 
    For QPM:
        Handles messy weights 
    '''
    catalog = cat_corr['catalog']
    
    # Tiling Mock --------------------------------------------------------------------------------
    if catalog['name'].lower() == 'tilingmock': 
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

    # QPM ----------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'qpm': 
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
        vetomask = (orig_true_veto == 0)            # only keep galaxies with w_veto = 1
        
        true_file = galaxy_data('data', readdata=False, **cat_corr)
        np.savetxt(true_file.file_name, np.c_[
            true_ra[vetomask], true_dec[vetomask], true_z[vetomask], 
            true_wfkp[vetomask], true_wfc[vetomask], true_comp[vetomask]], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    
    # Las Damas Geo ---------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'lasdamasgeo': 
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

    else: 
        raise NameError('not yet coded') 

def build_qpm_true_random(): 
    '''
    Hacked code to build the true random catalogs for convenience. 
    '''
    # read original random catalog  
    orig_true_random = np.loadtxt('/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/randoms/a0.6452_rand50x.dr12d_cmass_ngc.rdz')             # ra, dec, z, wfkp
    orig_true_random_info = np.loadtxt('/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/randoms/a0.6452_rand50x.dr12d_cmass_ngc.rdz.info')   # galid, comp?
    orig_true_random_veto = np.loadtxt('/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/randoms/a0.6452_rand50x.dr12d_cmass_ngc.veto')       # veto  

    vetomask = (orig_true_random_veto == 0)
    true_random = galaxy_data('random', readdata=False, **{'catalog':{'name':'qpm'}, 'correction':{'name':'true'}})
    
    np.savetxt(true_random.file_name, np.c_[(orig_true_random[:,0])[vetomask], (orig_true_random[:,1])[vetomask], 
        (orig_true_random[:,2])[vetomask], (orig_true_random[:,3])[vetomask], (orig_true_random_info[:,1])[vetomask]], 
        fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def build_fibercollided(**cat_corr): 
    '''
    build Fibercollided mock catalogs using specific idl routines
    '''
    catalog = cat_corr['catalog']
    
    # Las Damas Geo ---------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        # command is a currently bit hardcoded
        fibcollided_cmd = ' '.join(['idl', '-e', '/home/users/hahn/powercode/FiberCollision/LasDamas/Geo/ldg_fibcollmock_wcp_assign,', 
            str(catalog['n_mock']), ",'"+catalog['letter']+"'"])

        subprocess.call(fibcollided_cmd.split()) 

    # Tiling Mock --------------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'tilingmock': 
        fibcollided_cmd = ' '.join(['idl', '-e', '"', "build_wcp_assign, 'tilingmock'", '"'])
        os.system(fibcollided_cmd) 

    # QPM --------------------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'qpm': 
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
    else: 
        raise NameError('not yet coded') 

    return fibcollided_cmd 

def build_peakcorrected_fibcol(sanitycheck=False, **cat_corr): 
    '''
    Build peak corrected fibercollided mock catalogs (using cosmolopy) 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # cosmology *DEPENDS ON MOCK CATALOG!* ---------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        omega_m = 0.25
    elif catalog['name'].lower() == 'qpm': 
        omega_m = 0.31
    elif catalog['name'].lower() == 'tilingmock': 
        omega_m = 0.274
    else: 
        raise NameError('not yet coded!')
    print 'Omega_m = ', omega_m, 'Omega_L = ', 1.0-omega_m      # assumign flatness
    cosmo = {}
    cosmo['omega_M_0'] = omega_m 
    cosmo['omega_lambda_0'] = 1.0 - omega_m 
    cosmo['h'] = 0.7 
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 
    # ---------------------------------------------------------------------------------------
    # lambda function 
    if correction['fit'].lower() == 'gauss': 
        fit_func = lambda x, sig: np.exp(-0.5*x**2/sig**2)
    elif correction['fit'].lower() == 'expon': 
        fit_func = lambda x, sig: np.exp(-1.0*np.abs(x)/sig)
    elif correction['fit'].lower() == 'true': 
        pass 
    else: 
        raise NameError('correction fit has to be specified as gauss or expon') 

    # Las Damas Geo ---------------------------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        # LOS comoving distance for redshift limit
        comdis_lo = cosmos.distance.comoving_distance(0.16, **cosmo)*cosmo['h']
        comdis_hi = cosmos.distance.comoving_distance(0.44, **cosmo)*cosmo['h'] 
        n_mocks = 160
        
        # in case more than one catalog parameter is defined
        if (isinstance(catalog['n_mock'], int) == False) or (isinstance(catalog['letter'], str) == False): 
            raise NameError('only one mock can be corrected at a time') 
        
        # read in fiber collided galaxies
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr) 

        # set peak fraction         
        if correction['name'].lower() in ('allpeak', 'allpeakshot'): 
            f_peak = 1.0
        else: 
            f_peak = correction['fpeak'] 

        appended_ra = [] 
        appended_dec = [] 
        appended_z = [] 
        appended_weight = [] 
        upweight_again = []
            
        # if we want to check that the peak p(r) is generated properly
        if sanitycheck == True: 
            pr_test = [] 
    
        for i_mock in range(len(fibcoll_mock.weight)): 
            # go through every galaxy in fibercollided mock catalog
            while fibcoll_mock.weight[i_mock] > 1: 
                # for galaxies with wcp > 1
                fibcoll_mock.weight[i_mock] = fibcoll_mock.weight[i_mock]-1.0

                # LOS comoving distance of the galaxy 
                comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], **cosmo)*cosmo['h']
                
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
                        # compute the displacement within peak using best-fit ----------------------------------
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
                        # compute the displacement within peak using actual distribution ------------------------  
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
                    if (comdis_imock+rand2 > comdis_hi) or (comdis_imock+rand2 < comdis_lo): 
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

    # Tiling mock -----------------------------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() in ('tilingmock', 'qpm'): 
        # redshift limits comoving distance  
        comdis_lo = cosmos.distance.comoving_distance(0.43, **cosmo)*cosmo['h']         # in units of Mpc/h
        comdis_hi = cosmos.distance.comoving_distance(0.7, **cosmo)*cosmo['h'] 

        if catalog['name'].lower() == 'qpm':
            n_mocks = 100

        print 'D_c(z_min) = ', comdis_lo, ' Mpc/h'
        print 'D_c(z_max) = ', comdis_hi, ' Mpc/h'
        
        # read in mock with fibercollisions imposed
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr) 

        if catalog['name'].lower() == 'qpm': 
            fibcoll_mock.weight = fibcoll_mock.wfc            # fibercollisions weights are ultimately the weights I will be using here
                        
        # read in the true galaxies for tail portion  
        # the tail portion of the peak corrections will be generated similar to the mksample procedures for 
        # assigning the random catalog redshift 
        true_cat_corr = {'catalog':catalog, 'correction': {'name': 'true'}}
        true_data = galaxy_data('data', **true_cat_corr)
        
        if catalog['name'].lower() == 'qpm': 
            true_weight = true_data.wfc
        else: 
            true_weight = true_data.weight
        true_z = true_data.z

        true_weight_cum = true_weight.cumsum()
        true_weight_max = true_weight_cum.max()
        i_true = np.arange(0,len(true_weight_cum))

        # set peak fraction (f_peak in correction direction reflects the correct fpeak, 
        # but for all peak correction methods this is simply used as weight 
        if (correction['name'].lower() == 'allpeak') or (correction['name'].lower() == 'allpeakshot'): 
            f_peak = 1.0
        else: 
            f_peak = correction['fpeak'] 
    
        # list of data columns to be appended 
        appended_ra = [] 
        appended_dec = [] 
        appended_z = []
        appended_dlos = []
        appended_weight = [] 
        if catalog['name'].lower() == 'qpm': 
            appended_wfkp = [] 
            appended_comp = [] 

        upweight_again = []

        if sanitycheck == True: 
            # check that the peak p(r) is generated properly
            pr_test = [] 
    
        # go through each galaxy 
        for i_mock in range(len(fibcoll_mock.weight)): 
            # for galaxies with wcp > 1
            while fibcoll_mock.weight[i_mock] > 1: 
                fibcoll_mock.weight[i_mock] = fibcoll_mock.weight[i_mock]-1.0

                # LOS comoving distance of the galaxy 
                comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], **cosmo)*cosmo['h']
                
                rand_num = np.random.random(1)              # to determine whether the galaxy is in the peak or not 
               
                # peak ------------------------------------------------------------------------------------------------
                if rand_num < f_peak:          # if in the peak 
                    # keep ra and dec
                    appended_ra.append(fibcoll_mock.ra[i_mock])
                    appended_dec.append(fibcoll_mock.dec[i_mock])
                    if catalog['name'].lower() == 'qpm': 
                        appended_wfkp.append(fibcoll_mock.wfkp[i_mock]) 
                        appended_comp.append(fibcoll_mock.comp[i_mock])

                    # appended galaxy has weight of 1.0 
                    appended_weight.append(1.0)
                    
                    if correction['fit'].lower() in ('gauss', 'expon'): 
                        # compute the displacement within peak ----------------------------------
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
                        # compute the displacement within peak using actual distribution ------------------------  
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

                    # in case the displacement falls out of bound 
                    # NOTE: THIS IS NOT CORRECT, BUT IMPLEMENTED FOR SIMPLICITY 
                    if (comdis_imock+rand2 > comdis_hi) or (comdis_imock+rand2 < comdis_lo): 
                        collided_z = comdis2z(comdis_imock-rand2, **cosmo)
                    else: 
                        collided_z = comdis2z(comdis_imock+rand2, **cosmo)

                    appended_z.append(collided_z[0]) 
                    appended_dlos.append(rand2)
                    
                    if correction['name'] == 'allpeakshot': 
                        # for allpeak + shot noise correction 
                        # need to re-upweight the fibercollided galaxies by 1-fpeak 
                        upweight_again.append(i_mock)       # save index to be re-upweighted 

                    if sanitycheck == True:             
                        # append calculated LOS displacement 
                        pr_test.append(rand2) 

                # Tail -----------------------------------------------------------------------------------------------
                else:
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

                    elif correction['name'] == 'peakshot': 
                        # peak + shot noise correction 
                        # record galaxy index that need to be upweighted again 
                        upweight_again.append(i_mock)
                    
                    elif (correction['name'] == 'allpeak') or (correction['name'] == 'allpeakshot'): 
                        # for all peak correction methods, all galaxies should be in the peak!
                        raise NameError('should not happen')

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
        if catalog['name'].lower() == 'qpm': 
            fibcoll_mock.wfkp = np.concatenate([fibcoll_mock.wfkp, appended_wfkp])
            fibcoll_mock.comp = np.concatenate([fibcoll_mock.comp, appended_comp])
        fibcoll_mock.z = np.concatenate([fibcoll_mock.z, appended_z])

        np.savetxt('peak_corrected_dlos_values_'+correction['name'].lower()+'_'+correction['fit'].lower()+'.dat',
                np.c_[appended_dlos], fmt=['%10.5f']) 
    else: 
        raise NameError('Error here') 

    peakcorr = galaxy_data('data', readdata=False, **cat_corr)
    peakcorr_file = peakcorr.file_name 
    
    if catalog['name'].lower() == 'qpm': 
        # QPM has specific formatting
        np.savetxt(peakcorr_file, np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.wfkp, fibcoll_mock.weight, fibcoll_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    else: 
        np.savetxt(peakcorr_file, np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    if sanitycheck == True: 
        np.savetxt(peakcorr_file+'.sanitycheck', np.c_[pr_test], fmt=['%10.5f'], delimiter='\t') 

def build_vlospeak_fibcol(**cat_corr): 
    '''
    Build peak corrected fibercollided mock catalogs using vLOS rather than dLOS 
    '''
    c_speed = 299792.0 
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    if correction['fit'].lower() == 'gauss': 
        fit_func = lambda x, sig: np.exp(-0.5*x**2/sig**2)
    elif correction['fit'].lower() == 'expon': 
        fit_func = lambda x, sig: np.exp(-1.0*np.abs(x)/sig)
    elif correction['fit'].lower() == 'real': 
        pass
    else: 
        raise NameError('correction fit has to be specified as gauss or expon') 

    # Las Damas Geo ---------------------------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        zmin = 0.16
        zmax = 0.44
        n_mocks = 160 
        
        # in case more than one catalog parameter is defined
        if (isinstance(catalog['n_mock'], int) == False) or (isinstance(catalog['letter'], str) == False): 
            raise NameError('only one mock can be corrected at a time') 
        
        # read in fiber collided galaxies
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr) 

        # set peak fraction         
        f_peak = correction['fpeak'] 

        appended_ra = [] 
        appended_dec = [] 
        appended_z = [] 
        appended_weight = [] 
        upweight_again = []
            
        for i_mock in range(len(fibcoll_mock.weight)): 
            # go through every galaxy in fibercollided mock catalog
            while fibcoll_mock.weight[i_mock] > 1: 
                # for galaxies with wcp > 1
                fibcoll_mock.weight[i_mock] = fibcoll_mock.weight[i_mock]-1.0
                z_imock = fibcoll_mock.z[i_mock]

                rand_num = np.random.random(1) 
                if rand_num <= f_peak:          # if in the peak 
                    # keep ra and dec
                    appended_ra.append(fibcoll_mock.ra[i_mock])
                    appended_dec.append(fibcoll_mock.dec[i_mock])

                    # appended galaxy has weight of 1.0 
                    appended_weight.append(1.0)
                    
                    if correction['fit'].lower() in ('gauss', 'expon'): 
                        # compute the displacement within peak using 
                        # parameterized fits for the peak ----------------------------------
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
                    elif correction['fit'].lower() == 'real': 
                        vlos_comb_peak_file = ((fibcoll_mock.file_name).rsplit('/', 1))[0]+'/vLOS_norm_peak_dist_'+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined.dat'
                        vlos_mid, vlos_dist = np.loadtxt(vlos_comb_peak_file, unpack=True, usecols=[0,1])

                        vlos_cdf = vlos_dist.cumsum() 

                        rand1 = np.random.random(1) 
                        
                        cdf_closest_index = min(range(len(vlos_cdf)), key = lambda i: abs(vlos_cdf[i]-rand1[0])) 
                        closest_vlos = vlos_mid[cdf_closest_index] 
                       
                        try: 
                            closest_vloses
                        except NameError:
                            closest_vloses = [closest_vlos]
                        else: 
                            closest_vloses.append(closest_vlos)

                        rand2 = np.array([closest_vlos])

                    # in case the displacement falls out of bound 
                    if (z_imock+rand2/c_speed > zmax) or (z_imock+rand2/c_speed < zmin): 
                        collided_z = z_imock-rand2/c_speed
                    else: 
                        collided_z = z_imock+rand2/c_speed

                    appended_z.append(collided_z[0]) 

                else:                           # if not in the peak 
                    if correction['name'] in ('vlospeak', 'vlospeaknbar'): 
                        # randomly displace 
                        appended_ra.append(fibcoll_mock.ra[i_mock])
                        appended_dec.append(fibcoll_mock.dec[i_mock])
                        appended_weight.append(1.0) 
                        
                        rand1 = np.random.random(1) 
                        appended_z.append(0.16+rand1[0]*0.28)

                    elif correction['name'] == 'vlospeakshot': 
                        upweight_again.append(i_mock)

        if correction['name'] == 'vlospeakshot': 
            # re-upweighting for peak+shotnoise correction 
            for i_upweightagain in upweight_again: 
                fibcoll_mock.weight[i_upweightagain] = fibcoll_mock.weight[i_upweightagain]+1.

        print len(appended_ra), ' galaxies were peak corrected'
        fibcoll_mock.ra = np.concatenate([fibcoll_mock.ra, appended_ra])
        fibcoll_mock.dec = np.concatenate([fibcoll_mock.dec, appended_dec])
        fibcoll_mock.weight = np.concatenate([fibcoll_mock.weight, appended_weight])
        fibcoll_mock.z = np.concatenate([fibcoll_mock.z, appended_z])

    # Tiling mock -----------------------------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() in ('tilingmock', 'qpm'): 
        # redshift limits
        zmin = 0.43
        zmax = 0.7
        if catalog['name'].lower() == 'qpm': n_mocks = 100
        
        # read in mock with fibercollisions imposed
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr) 

        if catalog['name'].lower() == 'qpm': 
            fibcoll_mock.weight = fibcoll_mock.wfc            # fibercollisions weights are ultimately the weights I will be using here
                        
        # read in the true galaxies for tail portion  
        # the tail portion of the peak corrections will be generated similar to the mksample procedures for 
        # assigning the random catalog redshift 
        true_cat_corr = {'catalog':catalog, 'correction': {'name': 'true'}}
        true_data = galaxy_data('data', **true_cat_corr)
        
        if catalog['name'].lower() == 'qpm': 
            true_weight = true_data.wfc
        else: 
            true_weight = true_data.weight
        true_z = true_data.z

        true_weight_cum = true_weight.cumsum()
        true_weight_max = true_weight_cum.max()
        i_true = np.arange(0,len(true_weight_cum))

        # set peak fraction (f_peak in correction direction reflects the correct fpeak, 
        # but for all peak correction methods this is simply used as weight 
        f_peak = correction['fpeak'] 
    
        # list of data columns to be appended 
        appended_ra = [] 
        appended_dec = [] 
        appended_z = []
        appended_weight = [] 
        if catalog['name'].lower() == 'qpm': 
            appended_wfkp = [] 
            appended_comp = [] 

        upweight_again = []
    
        # go through each galaxy 
        for i_mock in range(len(fibcoll_mock.weight)): 
            # for galaxies with wcp > 1
            while fibcoll_mock.weight[i_mock] > 1: 
                fibcoll_mock.weight[i_mock] = fibcoll_mock.weight[i_mock]-1.0

                # redshift of the galaxy 
                z_imock = fibcoll_mock.z[i_mock] 
                
                rand_num = np.random.random(1)              # to determine whether the galaxy is in the peak or not 
                
                # peak ------------------------------------------------------------------------------------------------
                if rand_num < f_peak:          # if in the peak 
                    # keep ra and dec
                    appended_ra.append(fibcoll_mock.ra[i_mock])
                    appended_dec.append(fibcoll_mock.dec[i_mock])
                    if catalog['name'].lower() == 'qpm': 
                        appended_wfkp.append(fibcoll_mock.wfkp[i_mock]) 
                        appended_comp.append(fibcoll_mock.comp[i_mock]) 

                    # appended galaxy has weight of 1.0 
                    appended_weight.append(1.0)
                
                    if correction['fit'].lower() in ('gauss', 'expon'): 
                        # compute the displacement within peak ----------------------------------
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*correction['sigma']
                            peakpofr = fit_func(rand2, correction['sigma']) 
                    elif correction['fit'].lower() in ('real'): 
                        vlos_comb_peak_file = ((fibcoll_mock.file_name).rsplit('/', 1))[0]+'/vLOS_norm_peak_dist_'+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined.dat'
                        vlos_mid, vlos_dist = np.loadtxt(vlos_comb_peak_file, unpack=True, usecols=[0,1])

                        vlos_cdf = vlos_dist.cumsum() 

                        rand1 = np.random.random(1) 
                        
                        cdf_closest_index = min(range(len(vlos_cdf)), key = lambda i: abs(vlos_cdf[i]-rand1[0])) 
                        closest_vlos = vlos_mid[cdf_closest_index] 
                        
                        try: 
                            closest_vloses
                        except NameError:
                            closest_vloses = [closest_vlos]
                        else: 
                            closest_vloses.append(closest_vlos)
                        
                        rand2 = np.array([closest_vlos])

                    #----------------------------------- --------------------------------- 
                    # in case the displacement falls out of bound 
                    # NOTE: THIS IS NOT CORRECT, BUT IMPLEMENTED FOR SIMPLICITY 
                    if (z_imock+rand2/c_speed > zmax) or (z_imock+rand2/c_speed < zmin): 
                        collided_z = z_imock-rand2/c_speed
                    else: 
                        collided_z = z_imock+rand2/c_speed

                    appended_z.append(collided_z[0]) 

                # Tail -----------------------------------------------------------------------------------------------
                else:
                    if correction['name'] == 'vlospeaknbar': 
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

                    elif correction['name'] == 'vlospeakshot': 
                        # peak + shot noise correction 
                        # record galaxy index that need to be upweighted again 
                        upweight_again.append(i_mock)

        if correction['name'] == 'vlospeakshot': 
            # re-upweighting for peak+shotnoise correction 
            for i_upweightagain in upweight_again: 
                fibcoll_mock.weight[i_upweightagain] = fibcoll_mock.weight[i_upweightagain]+1.

        print len(appended_ra), ' galaxies were peak corrected'
        fibcoll_mock.ra = np.concatenate([fibcoll_mock.ra, appended_ra])
        fibcoll_mock.dec = np.concatenate([fibcoll_mock.dec, appended_dec])
        fibcoll_mock.weight = np.concatenate([fibcoll_mock.weight, appended_weight])
        if catalog['name'].lower() == 'qpm': 
            fibcoll_mock.wfkp = np.concatenate([fibcoll_mock.wfkp, appended_wfkp])
            fibcoll_mock.comp = np.concatenate([fibcoll_mock.comp, appended_comp])
        fibcoll_mock.z = np.concatenate([fibcoll_mock.z, appended_z])
    else: 
        raise NameError('Error here') 

    '''
    fig = plt.figure(1) 
    sub = fig.add_subplot(111)
    closest_vloses_hist, vlos_bins = np.histogram(closest_vloses, bins=300, range=[-5000.0, 5000.0], normed=True)
    vlos_low = vlos_bins[:-1]
    vlos_high = vlos_bins[1:] 
    vlos_mid = np.array([0.5*(vlos_low[i]+vlos_high[i]) for i in range(len(vlos_low))])
    sub.plot(vlos_mid, closest_vloses_hist) 
    plt.show() 
    '''

    peakcorr = galaxy_data('data', readdata=False, **cat_corr)
    peakcorr_file = peakcorr.file_name 
    
    if catalog['name'].lower() == 'qpm': 
        # QPM has specific formatting
        np.savetxt(peakcorr_file, np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.wfkp, fibcoll_mock.weight, fibcoll_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    else: 
        np.savetxt(peakcorr_file, np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def build_test_adjustments(sanitycheck=False, **cat_corr): 
    '''
    For any tests to galaxy mock catalogs 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    # LasDamasGeo-----------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        
        # in case more than one catalog parameter is defined
        if (isinstance(catalog['n_mock'], int) == False) or (isinstance(catalog['letter'], str) == False): 
            raise NameError('only one mock can be corrected at a time') 
    
        corrected_ra = [] 
        corrected_dec = [] 
        corrected_z = []
        corrected_weight = [] 

        # randomly remove galaxies
        if correction['name'].lower() == 'randrm': 

            # read in original mock catalog
            orig_cat_corr = {'catalog':catalog, 'correction': {'name': 'true'}}
            orig_mock = galaxy_data('data', **orig_cat_corr) 
    
            for i_mock in range(len(orig_mock.weight)):
                # Loop through galaxies 
                rand_num = np.random.random(1) 

                if rand_num > 0.03:             # hardcoded number of galaxies
                    corrected_ra.append(orig_mock.ra[i_mock])
                    corrected_dec.append(orig_mock.dec[i_mock])
                    corrected_z.append(orig_mock.z[i_mock])
                    corrected_weight.append(orig_mock.weight[i_mock])
                else: 
                    pass 
        else: 
            raise NameError('Yet to be coded') 

        corrected_ra = np.array(corrected_ra)
        corrected_dec = np.array(corrected_dec)
        corrected_z = np.array(corrected_z)
        corrected_weight = np.array(corrected_weight)

        peakcorr = galaxy_data('data', readdata=False, **cat_corr)
        peakcorr_file = peakcorr.file_name 
        np.savetxt(peakcorr_file, np.c_[corrected_ra, corrected_dec, corrected_z, corrected_weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

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

# --------------------------------------------------------------------
# Utility Functions
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
