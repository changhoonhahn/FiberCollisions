'''


Code to analyze nbar(z) of mock catalogs


Author(s) : ChangHoon Hahn 


'''

import numpy as np
import os.path
import subprocess
import cosmolopy as cosmos

# --- Local --- 
import fibcol_data as fc_data
import fibcol_utility as fc_util

# Classes ------------------------------------------------------------
class nbar: 
    def __init__(self, clobber=False, **cat_corr):
        ''' read/construct nbar(z) file corresponding to catalog_correction dictionary 

        Parameters
        ----------
        clobber : True/False
        cat_corr : catalog and correction dictionary 

        Notes
        -----
        * nbar(z) is built by first using a nbar-ngal file, then scaling the nbar(z) file 
        * Upweight, Peak+shot correciton methods change nbar(z) by a negligible amount so a correction is unnecssary
        * only coded for lasdamasgeo so far

        '''
        catalog = cat_corr['catalog'] 
        correction  = cat_corr['correction'] 

        # set zlow, zmid, zhigh (consistent with CMASS nbar(z) files) 
        self.zlow = np.arange(0.0, 1.005, 0.005)
        self.zhigh = self.zlow+0.005 
        self.zmid = self.zlow+0.0025 
        n_zbin = len(self.zlow) 
    
        if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo -----------------------------------------------

            self.file_name = get_nbar_file(**cat_corr)  # nbar file name 
            
            ldg_nbar = 0.0000944233     # constant nbar(z) value for true

            if correction['name'] .lower() in ('true', 'upweight', 'peakshot'): 
                # no correction 
                self.nbar = np.array([ldg_nbar for i in range(n_zbin)]) 

            else:
                # nbar(z) needs to be corrected 

                if (os.path.isfile(self.file_name) == True) and (clobber == False): 

                    self.nbar = np.loadtxt(self.file_name, unpack=True, usecols=[3])

                else: 

                    print 'Constructing ', self.file_name 

                    # check that the corrected nbar-ngal files exist
                    corr_ngal_file = get_nbar_ngal_file('allmocks', **cat_corr)

                    if (os.path.isfile(corr_ngal_file) == False) or clobber: 
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

                    # determine corrected nbar(z) using corrected nbar_ngal file
                    self.nbar = np.zeros(len(self.zlow))
                    for i in range(len(true_ngal)): 
                        if true_ngal[i] != 0: 
                            self.nbar[i] = ldg_nbar*(corr_ngal[i]/true_ngal[i])
                            #print (corr_ngal[i]/true_ngal[i])
                        else: 
                            self.nbar[i] = 0.0

                    self.writenbar()

        else: 
            raise NameError("Not yet Coded!") 

    def writenbar(self): 
        ''' Write class object nbar to ASCII file 
        '''
        np.savetxt(self.file_name, 
                np.c_[self.zmid, self.zlow, self.zhigh, self.nbar], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e'], delimiter='\t')

def build_averaged_nbar(n_mock, **cat_corr): 
    ''' Build averaged nbar(z) file for survey
    
    Parameters
    ----------
    n_mock : number of mock catalogs
    cat_corr : catalog correction dictionary

    Notes
    -----

    '''
    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction'] 
    
    # redshift bins 
    zlow = np.arange(0.0, 1.0, 0.005) 
    zhigh = zlow + 0.005
    zmid = zlow + 0.0025
    
    comp_area = 7341./41252.96
    
    mock_flags = [] 
    if catalog['name'].lower() == 'lasdamasgeo': 
        for i_mock in range(1, n_mock+1): 
            for letter in ['a', 'b', 'c', 'd']: 
                mock_flags.append({'n_mock': i_mock, 'letter': letter})
    else: 
        for i_mock in range(1, n_mock+1): 
            mock_flags.append({'n_mock': i_mock})

    
    # loop through mocks
    for i_m, mock_flag in enumerate(mock_flags): 

        i_catalog = catalog.copy()  
        i_catalog['n_mock'] = mock_flag['n_mock'] 
        if catalog['name'].lower() == 'lasdamasgeo': 
            i_catalog['letter'] = mock_flag['letter'] 

        i_cat_corr = {'catalog': i_catalog, 'correction': correction}
        # import galaxy data 
        data = fc_data.galaxy_data('data', **i_cat_corr)
        cosmo = data.cosmo  # survey cosmology 

        nbar = [] 
        for i_z in range(len(zmid)): 

            z_range = (data.z >= zlow[i_z]) & (data.z < zhigh[i_z])
            
            if catalog['name'].lower() == 'lasdamasgeo': 
                shell_weight =  np.sum(data.weight[z_range]) 
            else: 
                shell_weight = np.sum(data.wfc[z_range])

            shell_volume = comp_area * (cosmos.distance.comoving_volume(zhigh[i_z], **cosmo) - 
                cosmos.distance.comoving_volume(zlow[i_z], **cosmo))*cosmo['h']**3
            shell_nbar = shell_weight/shell_volume
            nbar.append(shell_nbar)

        nbar = np.array(nbar)
            
        if i_m == 0: 
            sum_nbar = nbar
        else: 
            sum_nbar += nbar 

    avg_nbar = sum_nbar/np.float(n_mock)

    # zcen,zlow,zhigh,nbar,wfkp,shell_vol,total weighted gals
    nbar_file = get_nbar_file(**cat_corr) 
    print nbar_file
    np.savetxt(nbar_file,
            np.c_[zmid, zlow, zhigh, avg_nbar],
            fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e'], delimiter='\t')

# Functions -----------------------------------------------------------------
def get_nbar_file(**cat_corr): 
    ''' Return nbar file name

    Parameters
    ----------
    DorR : 'data', 'random', 'allmocks' (You want to use 'allmocks') 
    cat_corr : catalog correction dictionary 
    
    Returns
    -------
    file_name 
    '''
    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction'] 
    
    if catalog['name'].lower() == 'lasdamasgeo': 
        catalog['n_mock'] = 1
        catalog['letter'] = 'a' 
    elif catalog['name'].lower() in ('qpm', 'nseries'): 
        catalog['n_mock'] = 1
    else: 
        raise NameError('not yet coded!')
    
    data_file = fc_data.get_galaxy_data_file('data', **cat_corr)
    data_dir = '/'.join(data_file.split('/')[:-1])+'/'      # directory
    data_file_name = data_file.split('/')[-1]

    # correction specifier 
    if catalog['name'].lower() == 'lasdamasgeo':
        corr_str = '.'.join(data_file_name.split('.')[2:])
    elif catalog['name'].lower() == 'qpm': 
        corr_str = '.'.join(data_file_name.split('.')[4:])
    elif catalog['name'].lower() == 'nseries':
        corr_str = '.'.join(data_file_name.split('.')[1:])
    else: 
        raise NameError('not yet coded!')

    # combine to form filename  
    file_name = ''.join([data_dir, 'nbar-', catalog['name'].lower(), '-', corr_str])
    return file_name

def get_nbar_ngal_file(DorR, **cat_corr): 
    ''' Return nbar_ngal file name

    Parameters
    ----------
    DorR : 'data', 'random', 'allmocks' (You want to use 'allmocks') 
    cat_corr : catalog correction dictionary 
    
    Returns
    -------
    file_name 
    '''
    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction'] 
    
    if DorR.lower() == 'allmocks': 
        data_file = fc_data.get_galaxy_data_file('data', **cat_corr)
    else: 
        data_file = fc_data.get_galaxy_data_file(DorR, **cat_corr)
    data_dir = '/'.join(data_file.split('/')[:-1])+'/'      # directory
    data_file_name = data_file.split('/')[-1]

    if DorR.lower() in ('data', 'random', 'allmocks'): 
        DorR_str = DorR.lower() 
    else: 
        raise NameError('DorR error') 

    if catalog['name'].lower() == 'lasdamasgeo':

        if DorR.lower() == 'data':
            catalog_str = ''.join([catalog['name'].lower(), '-', str(catalog['n_mock']), catalog['letter']]) 
        elif DorR.lower() == 'allmocks': 
            catalog_str = ''.join([catalog['name'].lower(), '-', 'allmocks']) 
        else:
            catalog_str = ''.join([catalog['name'].lower(), '-', 'random']) 
        
        # correction specifier 
        corr_str = '.'.join(data_file_name.split('.')[2:])

    elif catalog['name'].lower() == 'qpm': 
        
        if DorR.lower() == 'data':
            catalog_str = ''.join([catalog['name'].lower(), '-', str(catalog['n_mock'])]) 
        elif DorR.lower() == 'allmocks': 
            catalog_str = ''.join([catalog['name'].lower(), '-', 'allmocks']) 
        else: 
            catalog_str = ''.join([catalog['name'].lower(), '-', 'random']) 
        
        # correction specifier 
        corr_str = '.'.join(data_file_name.split('.')[4:])

    else: 
        raise NameError('not yet coded!')

    # combine to form filename  
    file_name = ''.join([data_dir, 
        'nbar-ngal-', catalog_str, '-', corr_str])
    return file_name

def write_nbar_ngal(DorR, **cat_corr): 
    ''' Construct ngal values for nbar(z) redshift bins for specified mock catalog
    so that it doesn't have to be repeated
    
    Parameters
    ----------
    DorR : data, random, allmocks 
    cat_corr : catalog correction specifier 
    
    Notes
    -----
    * If DorR == 'allmock', then nbar_ngal is calculated for the combination of all the mocks

    '''
    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction'] 
   
    # import z values and shell volumes from nbar file  
    zcen, zlow, zhigh, shell_vol = np.loadtxt('/mount/riachuelo1/hahn/data/nbar-junk.dat', 
            unpack=True, usecols=[0,1,2,5])
    z_values =[zcen, zlow, zhigh, shell_vol]

    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ----------------------
        
        if DorR.lower() == 'allmocks':  # combine mocks 

            Ngal = 0 
            for i_mock in range(1,41): 
                for letter in ['a', 'b', 'c','d']:
                    i_catalog = catalog 
                    i_catalog['n_mock'] = i_mock 
                    i_catalog['letter'] = letter
                    i_cat_corr = {'catalog': i_catalog, 'correction':correction}
                    
                    i_data = fc_data.galaxy_data('data', **i_cat_corr) 
                    print 'Reading ', i_data.file_name 
                   
                    Ngal += len(i_data.z) 

                    # save redshifts and weights
                    try: 
                        z_dist = np.concatenate([z_dist, i_data.z]) 
                        z_weights = np.concatenate([z_weights, i_data.weight]) 
                    except NameError: 
                        z_dist = i_data.z
                        z_weights = i_data.weight

            print Ngal 
            print len(z_dist) 
            print len(z_weights) 

        else:       # data or random 

            data = fc_data.galaxy_data(DorR, **cat_corr) 
    
            z_dist = data.z 
            try:  
                data.weight                         # see if weights exists (does not exist for certain randoms and trues) 
            except AttributeError: 
                z_weights = np.array([1.0 for i in range(len(z_dist))]) 
            else: 
                z_weights = data.weight

        # calculate ngal(z) 
        nbar_ngal = np.zeros(len(z_values[0]))
        for i_z, zmid in enumerate(z_values[0]):
            zlim = (z_dist >= (z_values[1])[i_z]) & (z_dist < (z_values[2])[i_z])
            nbar_ngal[i_z] = np.sum(z_weights[zlim])

    elif catalog['name'].lower() in ('cmass', 'tilingmock'):        # Tiling Mock/CMASS --------------------------------

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

    elif catalog['name'].lower() == 'qpm':                          # QPM ----------------------------------------
        
        if DorR.lower() == 'allmocks':  # combine mocks 

            Ngal = 0 
            for i_mock in range(1,51):  # (only 50 mocks for now) 
                i_catalog = catalog 
                i_catalog['n_mock'] = i_mock 
                i_cat_corr = {'catalog': i_catalog, 'correction':correction}
                i_data = fc_data.galaxy_data('data', **i_cat_corr) 
                print 'Reading ', i_data.file_name 
                   
                Ngal += len(i_data.z) 
                try: 
                    z_dist
                except NameError: 
                    z_dist = i_data.z
                    z_weights = i_data.weight
                else: 
                    z_dist = np.concatenate([z_dist, i_data.z]) 
                    z_weights = np.concatenate([z_weights, i_data.weight]) 

            if (Ngal != len(z_dist)) or (Ngal != len(z_weights)): 
                raise NameError("mismatch between Ngal and z_weights") 

        else:
            data = fc_data.galaxy_data(DorR, **cat_corr) 
    
            z_dist = data.z 
            
            # check whether weights exists
            try:  
                data.weight                         
            except AttributeError: 
                z_weights = np.array([1.0 for i in range(len(z_dist))]) 
            else: 
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

def append_corr_nbar(DorR, **cat_corr):
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

if __name__=='__main__':
    cat_corr = {'catalog': {'name': 'lasdamasgeo'}, 
            'correction': {'name': 'true'}} 
    build_averaged_nbar(10, **cat_corr)
    #nbar(**cat_corr) 
