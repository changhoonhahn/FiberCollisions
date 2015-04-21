import numpy as np
import os.path
import subprocess
import cosmolopy as cosmos
import fibcol_data as fc_data
import fibcol_utility as fc_util

# Classes ------------------------------------------------------------
class nbar: 
    def __init__(self, **cat_corr):
        '''
        read/construct nbar(z) file corresponding to catalog_correction dictionary 
        '''
        catalog = cat_corr['catalog'] 
        correction  = cat_corr['correction'] 

        # set zlow, zmid, zhigh (consistent with CMASS nbar(z) files) 
        self.zlow = np.array([0.005*np.float(i) for i in range(201)])
        self.zhigh = self.zlow+0.005 
        self.zmid = self.zlow+0.0025 
    
        # store catalog and correction info 
        self.cat_corr = cat_corr
        
        # LasDamasGeo ----------------------------------------------------------------------------------------
        if catalog['name'].lower() == 'lasdamasgeo': 
            # get nbar file name 
            self.file_name = get_nbar_file(**cat_corr) 

            ldg_nbar = 9.44451*10.0**-5         # constant nbar(z) value for true

            # for no correction ----------------------------------------------------------------------------
            if correction['name'] .lower() in ('true', 'upweight', 'peakshot', 'vlospeakshot'): 
                # upweight and peakshot change the nbar(z) only by a negligible amount so a correction is unnecesary   
                # constant nbar(z)=9.44451e-5 
                self.nbar = np.array([ldg_nbar for i in range(len(self.zlow))]) 

            # for upweight, peak, all peak, and peak test correction ----------------------------------------
            else:
                # if file exists
                if os.path.isfile(self.file_name) == True: 
                    # read file 
                    self.nbar = np.loadtxt(self.file_name, unpack=True, usecols=[3])
                else: 
                    # if nbar does not exist then make it 
                    print 'Constructing ', self.file_name 

                    # check that the corrected nbar-ngal files exist
                    corr_ngal_file = get_nbar_ngal_file('allmocks', **cat_corr)

                    if os.path.isfile(corr_ngal_file) == False: 
                        print 'Constructing ', corr_ngal_file
                        write_nbar_ngal('allmocks', **cat_corr) 
                    else: 
                        print 'Using ', corr_ngal_file

                    # read corrected nbar_ngal
                    corr_ngal = np.loadtxt(corr_ngal_file, unpack=True, usecols=[3]) 
       
                    # check that the true nbar_ngal fiel exists
                    true_cat_corr = {'catalog':catalog, 'correction':{'name':'true'}}
                    true_ngal_file = get_nbar_ngal_file('allmocks', **true_cat_corr) 

                    if os.path.isfile(true_ngal_file) == False: 
                        print 'Constructing ', true_ngal_file
                        write_nbar_ngal('allmocks', **true_cat_corr) 
                    
                    # read true nbar_ngal
                    true_ngal = np.loadtxt(true_ngal_file, unpack=True, usecols=[3])

                    # use corrected nbar_ngal file to determine corrected nbar(z) 
                    self.nbar = np.zeros(len(self.zlow))
                    for i in range(len(true_ngal)): 
                        if true_ngal[i] != 0: 
                            self.nbar[i] = ldg_nbar*(corr_ngal[i]/true_ngal[i])
                            #print (corr_ngal[i]/true_ngal[i])
                        else: 
                            self.nbar[i] = 0.0
                    self.writenbar()
        
        # Tiling MOck ------------------------------------------------------------------------------------ 
        elif catalog['name'].lower() == 'tilingmock': 
            # get file name 
            if correction['name'].lower() in ('true', 'upweight', 'peakshot', 'vlospeakshot'): 
                self.file_name = '/mount/riachuelo1/hahn/data/tiling_mocks/nbar-cmass-boss5003sector-icoll012.dat'  #hard coded
            else: 
                self.file_name = get_nbar_file(**cat_corr)          

            if os.path.isfile(self.file_name) == True: 
                self.nbar = np.loadtxt(self.file_name, unpack=True, usecols=[3])

            # if nbar does not exist then make it 
            else:
                print 'Constructing ', self.file_name 
                # check that the corrected nbar-ngal files exist
                corr_ngal_file = get_nbar_ngal_file('random', **cat_corr)

                if os.path.isfile(corr_ngal_file) == False: 
                    print 'Constructing ', corr_ngal_file
                    write_nbar_ngal('random', **cat_corr) 
                else: 
                    print 'Using ', corr_ngal_file

                # read corrected nbar_ngal
                corr_ngal = np.loadtxt(corr_ngal_file, unpack=True, usecols=[3]) 
       
                # check that the true nbar_ngal fiel exists
                true_cat_corr = {'catalog':catalog, 'correction':{'name':'true'}}
                true_ngal_file = get_nbar_ngal_file('random', **true_cat_corr) 

                if os.path.isfile(true_ngal_file) == False: 
                    print 'Constructing ', true_ngal_file
                    write_nbar_ngal('random', **true_cat_corr) 
                
                true_nbar_file = get_nbar_file(**true_cat_corr)         # true nbar file name 

                # read true nbar_ngal and nbar files 
                true_ngal = np.loadtxt(true_ngal_file, unpack=True, usecols=[3])
                true_nbar = np.loadtxt(true_nbar_file, unpack=True, usecols=[3])
                
                # use corrected nbar_ngal file to determine corrected nbar(z) 
                self.nbar = np.zeros(len(self.zlow))
                for i in range(len(true_ngal)): 
                    if true_ngal[i] != 0: 
                        self.nbar[i] = true_nbar[i]*(corr_ngal[i]/true_ngal[i])
                    else: 
                        self.nbar[i] = 0.0
                self.writenbar()

        else: 
            raise NameError("Not yet Coded!") 

    def writenbar(self): 
        '''
        Write class object nbar to ASCII file 
        '''
        np.savetxt(self.file_name, 
                np.c_[self.zmid, self.zlow, self.zhigh, self.nbar], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e'], delimiter='\t')

# Functions -----------------------------------------------------------------
def get_nbar_file(**cat_corr): 
    '''
    get nbar file name
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    # las damas geo ------------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        # construct nbar file 
        nbar_file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'     # nbar directory 
        
        # correction string 
        if correction['name'].lower() in ('true', 'upweight', 'shotnoise', 'floriansn', 'hectorsn'): 
            # ture or upweight
            corr_param_str = ''
        else: 
            try: 
                correction['fit']
            except KeyError: 
                pass 
            else: 
                # specify peak correction fit (expon or gauss) 
                if correction['fit'].lower() in ('gauss', 'expon', 'real'): fit_str = correction['fit'].lower()  
                else: raise NameError('fit type not specified') 
            
            if correction['name'].lower() in ('peak', 'peaknbar','peaktest', 'peakshot', 'vlospeakshot'): 
                # for anything that involves the peak and tail  

                if correction['name'].lower() == 'peak':   
                    # correct for poor namign convention 
                    correction['name'] = 'peaknbar'
                
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    corr_param_str = ''.join(['.', fit_str, '.', correction['name'].lower(), 
                        '.sigma', str(correction['sigma']), 'fpeak', str(correction['fpeak'])])
                elif correction['fit'].lower() in ('real'): 
                    corr_param_str = ''.join(['.', fit_str, '.', correction['name'].lower(), '.fpeak', str(correction['fpeak'])])

            # for anything that involves only the peak
            elif (correction['name'].lower() == 'allpeak') or (correction['name'].lower() == 'allpeakshot'): 
                corr_param_str = ''.join(['.', fit_str, '.', correction['name'].lower(), 
                    '.sigma', str(correction['sigma'])])

            # test adjustments  
            elif (correction['name'].lower() == 'randrm'): 
                corr_param_str = ''

            else: 
                raise NameError('Correction method not supported') 

        if correction['name'].lower() == 'shotnoise': 
            corr_str = 'upweight' 
        else: 
            corr_str = correction['name'].lower() 
        file_name = nbar_file_dir+'-'.join(['nbar', catalog['name'].lower(), corr_str])+corr_param_str+'.dat'

    # Tiling Mock ------------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'tilingmock': 
        # construct nbar file 
        nbar_file_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'     # nbar directory 
        
        # correction string 
        if correction['name'].lower() in ('true', 'upweight', 'shotnoise', 'floriansn', 'hectorsn'): 
            # true or upweight
            corr_param_str = ''

        # Any correction emthod involving the peak  --------------------------------------------------
        else: 
            try: 
                correction['fit']
            except KeyError: 
                pass 
            else: 
                # specify peak correction fit (expon or gauss) 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    fit_str = correction['fit'].lower()  
                else: 
                    raise NameError('fit type not specified') 
            
            # correction specifying string: 
            if correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'vlospeakshot'): 
                # for anything that involves the peak and tail  

                if correction['name'].lower() == 'peak': 
                    # correct for bad naming convention 
                    correction['name'] = 'peaknbar'
                
                corr_param_str = ''.join(['.', fit_str, '.', correction['name'].lower(), 
                    '.sigma', str(correction['sigma']), 'fpeak', str(correction['fpeak'])])

            elif correction['name'].lower() in ('allpeak', 'allpeakshot'): 
                # for anything that involves only the peak

                corr_param_str = ''.join(['.', fit_str, '.', correction['name'].lower(), 
                    '.sigma', str(correction['sigma'])])

            elif (correction['name'].lower() == 'randrm'): 
                # test adjustments  
                corr_param_str = ''

            else: 
                print correction['name'].lower()
                raise NameError('Correction method not supported') 


        file_name = nbar_file_dir+'-'.join(['nbar', catalog['name'].lower(), correction['name']])+corr_param_str+'.dat'
        
        if correction['name'].lower() == 'true': 
            # crude hardcoded for true correction method
            file_name = '/mount/riachuelo1/hahn/data/tiling_mocks/nbar-cmass-boss5003sector-icoll012.dat'

    else: 
        raise NameError('catalog not yet coded') 

    return file_name 

def get_nbar_ngal_file(DorR, **cat_corr): 
    '''
    get nbar_ngal file name
    '''
    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction'] 
    
    file_prefix = 'nbar-ngal-'

    catalog_str = catalog['name'].lower() 
    if DorR.lower() in ('data', 'random', 'allmocks'): 
        DorR_str = DorR.lower()             # can be 'data', 'random', and 'allmocks'
    else: 
        raise NameError('DorR error') 

    # LasDamasGeo ------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'

        # if data, specify specific catalog #  
        if DorR.lower() == 'data':
            catalog_str = catalog_str+'-'+str(catalog['n_mock'])+catalog['letter']
    
        # specify correction method 
        # No extra correction parameters ---------------------------------------------------------------------
        if correction['name'].lower() in ('true', 'upweight', 'shotnoise', 'floriansn', 'hectorsn'): 
            if correction['name'].lower() in ('shotnoise', 'floriansn', 'hectorsn'): 
                corr_str = 'upweight'
            else:
                corr_str = correction['name'].lower()  

        # correction methods that involve peak and tail ---------------------------------------------
        elif correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'vlospeakshot'):

            if correction['name'].lower() == 'peak': 
                # correct for poor naming convention
                correction['name'] = 'peaknbar'

            # specify peak correction fit (expon or gauss) 
            if correction['fit'].lower() in ('gauss', 'expon'): 
                corr_str = ''.join([correction['fit'].lower(), '.', correction['name'].lower(), 
                    '.sigma', str(correction['sigma']), 'fpeak', str(correction['fpeak'])])
            elif correction['fit'].lower() in ('real'): 
                corr_str = ''.join([correction['fit'].lower(), '.', correction['name'].lower(), '.fpeak', str(correction['fpeak'])]) 
            else: 
                raise NameError('peak fit has to be specified: gauss or expon') 

        # correction methods that involve only peak -------------------------------------------------
        elif correction['name'].lower() in ('allpeak', 'allpeakshot'): 
            
            # specify peak correction fit (expon or gauss) 
            if correction['fit'].lower() in ('gauss', 'expon'): 
                pass
            else: 
                raise NameError('peak fit has to be specified: gauss or expon') 

            corr_str = ''.join([correction['fit'].lower(), '.', correction['name'].lower(), '.sigma', str(correction['sigma'])])
        
        # test adjustments ----------------------------------------------------------------------------
        elif (correction['name'].lower() == 'randrm'): 
            corr_str = correction['name'].lower() 

        else: 
            raise NameError('correction method unknown')
    
    # Tiling mock --------------------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'tilingmock': 
        file_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'
        
        # specify correction method 
        # No correction parameters ----------------------------------------------------------------------------------------
        if correction['name'].lower() in ('true', 'upweight', 'shotnoise', 'floriansn', 'hectorsn', 'randrm'): 
            if correction['name'].lower() in ('shotnoise', 'floriansn', 'hectorsn'):
                corr_str = 'upweight'
            else:
                corr_str = correction['name'].lower()  

        # correction methods that involve peak and tail ---------------------------------------------
        elif correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'vlospeakshot'):

            if correction['name'].lower() == 'peak': 
                # correct for poor naming convention
                correction['name'] = 'peaknbar'

            # specify peak correction fit (expon or gauss) 
            if correction['fit'].lower() in ('gauss', 'expon'): 
                pass
            else: 
                raise NameError('peak fit has to be specified: gauss or expon') 

            corr_str = ''.join([correction['fit'].lower(), '.', correction['name'].lower(), 
                '.sigma', str(correction['sigma']), 'fpeak', str(correction['fpeak'])])

        # correction methods that involve only peak -------------------------------------------------
        elif correction['name'].lower() in ('allpeak', 'allpeakshot'): 
            
            # specify peak correction fit (expon or gauss) 
            if (correction['fit'].lower() == 'gauss') or (correction['fit'].lower() == 'expon'): 
                pass
            else: 
                raise NameError('peak fit has to be specified: gauss or expon') 

            corr_str = ''.join([correction['fit'].lower(), '.', correction['name'].lower(), '.sigma', str(correction['sigma'])])
        else: 
            raise NameError('correction method not available')

    elif catalog['name'].lower() == 'qpm': 
        file_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'

        if correction['name'].lower() in ('true'): 
            corr_str = correction['name'].lower()  

    elif catalog['name'].lower() == 'cmass': 
        file_dir = '/mount/riachuelo1/hahn/data/'
        #####  
        #####  
        #####  FIX LATER
        #####  
        #####  
        corr_str = correction['name'].lower() 
    else: 
        raise NameError('not yet coded!')

    # combine to form filename  
    file_name = ''.join([file_dir, file_prefix, catalog_str, '-', DorR_str, '-', corr_str, '.dat'])
    return file_name

def write_nbar_ngal(DorR, **cat_corr): 
    '''
    write ngal values for nbar(z) redshift bins to a file so it doesn't have to be repeated
    '''
    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction'] 
   
    # import z values and shell volumes from arbitrary nbar file  
    zcen, zlow, zhigh, shell_vol = np.loadtxt('/mount/riachuelo1/hahn/data/nbar-junk.dat', unpack=True, usecols=[0,1,2,5])
    z_values =[zcen, zlow, zhigh, shell_vol]

    # LasDamasGeo --------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        # import mock/random data  
        if DorR.lower() == 'allmocks': 
            Ngal = 0 
            for i_mock in range(1,41): 
                for letter in ['a', 'b', 'c','d']:
                    i_catalog = catalog 
                    i_catalog['n_mock'] = i_mock 
                    i_catalog['letter'] = letter
                    i_cat_corr = {'catalog': i_catalog, 'correction':correction}
                    i_data = fc_data.galaxy_data('data', **i_cat_corr) 
                    print 'Reading ', i_data.file_name 
                   
                    Ngal = Ngal + len(i_data.z) 
                    try: 
                        z_dist
                    except NameError: 
                        z_dist = i_data.z
                        z_weights = i_data.weight
                    else: 
                        z_dist = np.concatenate([z_dist, i_data.z]) 
                        z_weights = np.concatenate([z_weights, i_data.weight]) 
            print Ngal 
            print len(z_dist) 
            print len(z_weights) 
        else:
            data = fc_data.galaxy_data(DorR, **cat_corr) 
    
            z_dist = data.z 
            try:  
                data.weight                         # see if weights exists (does not exist for certain randoms and trues) 
            except AttributeError: 
                Ngal = np.float(len(data.z))
                z_weights = np.array([1.0 for i in range(len(z_dist))]) 
            else: 
                Ngal = np.sum(data.weight)          # number of galaxies account for weights 
                z_weights = data.weight

        nbar_ngal = np.zeros(len(z_values[0]))
        for i_z, zmid in enumerate(z_values[0]):
            zlim = (z_dist >= (z_values[1])[i_z]) & (z_dist < (z_values[2])[i_z])
            nbar_ngal[i_z] = np.sum(z_weights[zlim])

    # Tiling Mock/CMASS ------------------------------------------------------------------------------------------------- 
    elif catalog['name'].lower() in ('cmass', 'tilingmock'): 
        # import mock/random data  
        if DorR.lower() == 'allmocks': 
            print catalog['name'].lower(), ' only has one mock'
            data = fc_data.galaxy_data('data', **cat_corr) 
        else:
            data = fc_data.galaxy_data(DorR, **cat_corr) 
    
        z_dist = data.z 
        try:  
            data.weight                         # see if weights exists (does not exist for certain randoms and trues) 
        except AttributeError: 
            Ngal = np.float(len(data.z))
            z_weights = np.array([1.0 for i in range(len(z_dist))]) 
        else: 
            Ngal = np.sum(data.weight)          # number of galaxies account for weights 
            z_weights = data.weight

        nbar_ngal = np.zeros(len(z_values[0]))
        for i_z, zmid in enumerate(z_values[0]):
            zlim = (z_dist >= (z_values[1])[i_z]) & (z_dist < (z_values[2])[i_z])
            nbar_ngal[i_z] = np.sum(z_weights[zlim])

    # QPM ----------------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'qpm': 
        # import mock/random data  
        if DorR.lower() == 'allmocks': 
            Ngal = 0 
            for i_mock in range(1,45): 
                i_catalog = catalog 
                i_catalog['n_mock'] = i_mock 
                i_cat_corr = {'catalog': i_catalog, 'correction':correction}
                i_data = fc_data.galaxy_data('data', **i_cat_corr) 
                print 'Reading ', i_data.file_name 
                   
                Ngal = Ngal + len(i_data.z) 
                try: 
                    z_dist
                except NameError: 
                    z_dist = i_data.z
                    z_weights = i_data.weight
                else: 
                    z_dist = np.concatenate([z_dist, i_data.z]) 
                    z_weights = np.concatenate([z_weights, i_data.weight]) 
            print Ngal 
            print len(z_dist) 
            print len(z_weights) 
        else:
            data = fc_data.galaxy_data(DorR, **cat_corr) 
    
            z_dist = data.z 
            try:  
                data.weight                         # see if weights exists (does not exist for certain randoms and trues) 
            except AttributeError: 
                Ngal = np.float(len(data.z))
                z_weights = np.array([1.0 for i in range(len(z_dist))]) 
            else: 
                Ngal = np.sum(data.weight)          # number of galaxies account for weights 
                z_weights = data.weight

        nbar_ngal = np.zeros(len(z_values[0]))
        for i_z, zmid in enumerate(z_values[0]):
            zlim = (z_dist >= (z_values[1])[i_z]) & (z_dist < (z_values[2])[i_z])
            nbar_ngal[i_z] = np.sum(z_weights[zlim])
    
    else: 
        raise NameError('Im an error') 

    # write nbar_ngal data to ask ascii file
    nbar_ngal_file = get_nbar_ngal_file(DorR, **cat_corr)
    np.savetxt(nbar_ngal_file,            
            np.c_[z_values[0], z_values[1], z_values[2], nbar_ngal],
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t')

def append_corr_nbar(DorR, sanitycheck=False, **cat_corr):
    '''
    append corrected interpolated nbar(z) to corrected data or random file
    '''
    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction'] 

    # read in data/random galaxies
    gal_data = fc_data.galaxy_data(DorR, **cat_corr)

    if DorR.lower() == 'random':
        print 'random ', gal_data.file_name, ' read' 
        print 'Nran = ', len(gal_data.ra)

    # read in corrected nbar file
    corr_nbar = nbar(**cat_corr)
    print corr_nbar.file_name

    # interpolate within redshift limits
    if catalog['name'].lower() == 'lasdamasgeo':
        zlim = (corr_nbar.zmid > 0.16) & (corr_nbar.zmid < 0.44)        # for las damas geo
    elif catalog['name'].lower() == 'tilingmock': 
        zlim = (corr_nbar.zmid > 0.43) & (corr_nbar.zmid < 0.7)        # for tilingmock  
    else:
        raise NameError("not yet coded!")

    # numpy interpolate
    nbar_arr = np.interp(gal_data.z, corr_nbar.zmid[zlim], corr_nbar.nbar[zlim])

    if DorR.lower() == 'random':
        print 'Nran = ', len(nbar_arr)

    # if santiy check is true then plot the interpolated nbar values
    if sanitycheck == True:
        fig = plt.figure(1)
        sub = fig.add_subplot(111)
        sub.scatter(gal_data.z, nbar_arr, color='red', s=5, label=DorR)
        sub.plot(corr_nbar.zmid, corr_nbar.nbar, color='black', lw=2, label="Corrected nbar(z)")
        sub.set_xlim([0.16, 0.44])
        sub.set_ylim([9.3e-05, 10.0e-05])
        sub.legend(loc='upper right', scatterpoints=1)
        fig.savefig(''.join(['/home/users/hahn/research/figures/boss/',
        'corrected-', DorR, '-nbarz-sanitycheck.png']), bbox_inches="tight")
        fig.clear()

    # write corr nbar appended data/random
    # corr nbar appended file name
    gal_corr_nbar_file = ''.join([gal_data.file_name, '.corrnbar'])
    if DorR.lower() == 'data':
        np.savetxt(gal_corr_nbar_file,
                np.c_[gal_data.ra, gal_data.dec, gal_data.z, nbar_arr, gal_data.weight],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e', '%10.5f'], delimiter='\t')
    elif DorR.lower() == 'random':
        np.savetxt(gal_corr_nbar_file,
                np.c_[gal_data.ra, gal_data.dec, gal_data.z, nbar_arr],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e'], delimiter='\t')
