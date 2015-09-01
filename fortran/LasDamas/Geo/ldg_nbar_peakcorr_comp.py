import os.path 

def readnbarz(file="/mount/riachuelo1/hahn/data/nbar-junk.dat"): 
    '''
    Read in nbar(z) file in order to get redshift arrays and comoving volume measurements (for consistency)  
    outputs [zcen, zlow, zhigh, shell_vol]
    '''
    zcen, zlow, zhigh, shell_vol = np.loadtxt(file, unpack=True, usecols=[0,1,2,5])         # read columns
    return [zcen, zlow, zhigh, shell_vol] 

def write_nbar_ngal(DorR='data', catalog='lasdamasgeo', catalog_param=[1,'a'], 
        corr='peaknbar', corr_param=[5.3, 0.1]): 
    '''
    write ngal values for nbar(z) redshift bins to a file so it doesn't have to be repeated 
    caution: hacked together
    '''
    z_values = readnbarz(file="/mount/riachuelo1/hahn/data/nbar-junk.dat")                  # import z-values

    if DorR.lower() == 'random': 
        # True Random catalog
        if corr.lower() == 'true': 
            # Hardcoded for LasDamasGeo
            if catalog.lower() == 'lasdamasgeo':
                true_rand_file = '/mount/chichipio2/rs123/MOCKS/randoms/sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat'
            else: 
                raise NameError('Only LasDamasGeo coded so far') 

            # read in files 
            true_rand_ra, true_rand_dec, true_rand_cz = np.loadtxt(true_rand_file, unpack=True, usecols=[0,1,2])
            
            # get z values and the weights
            z_dist = true_rand_cz/299800.0
            z_weights = np.array([1 for i in range(len(z_dist))])

        # Peak corrected Random catalog 
        elif corr.lower() == 'peaknbar': 
            if catalog.lower() == 'lasdamasgeo':
                corr_str = '.peak.sigma'+str(corr_param[0])+'.fpeak'+str(corr_param[1])
                corr_rand_file = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'+\
                        'sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat.allmocks'+corr_str 
            else: 
                raise NameError('Only LasDamasGeo coded so far') 
            # read corrected ranomd file 
            corr_rand_ra, corr_rand_dec, corr_rand_z = np.loadtxt(corr_rand_file, unpack=True, usecols=[0,1,2])
            
            # get z values and weights
            z_dist = corr_rand_z
            z_weights = np.array([1 for i in range(len(z_dist))])

    # data catalog         
    elif DorR.lower() == 'data': 
        # For true data
        if corr.lower() == 'true': 
            if catalog.lower() == 'lasdamasgeo': 
                true_file = ''.join(['/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/', 
                    'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog_param[0]), catalog_param[1], 
                    '_no.rdcz.dat'])
            # read true file 
            true_ra, true_dec, true_cz = np.loadtxt(true_file, unpack=True, usecols=[0,1,2]) 
            
            # get z values and weights
            z_dist = true_cz/299800.0
            z_weights = np.array([1 for i in range(len(z_dist))])

        # For peak corrected data 
        elif corr.lower() == 'peaknbar': 
            if catalog.lower() == 'lasdamasgeo': 
                if len(catalog_param) == 2: 
                    if corr.lower() == 'peaknbar':            # peak nbar corrected
                        fibcoll_str = 'fibcoll.'
                        corr_str = '.peak.sigma'+str(corr_param[0])+'.fpeak'+str(corr_param[1])
                    corr_file = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/', 
                        'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog_param[0]), catalog_param[1], 
                        '_no.rdcz.', fibcoll_str, 'dat', corr_str])
                else:  
                    raise NameError("Specify LasDamasGeo Catalog Parameters (e.g. [1,'a'])") 
            else: 
                raise NameError('Only LasDamasGeo coded so far')

            # read-in corrected file
            corr_ra, corr_dec, corr_z, corr_w = np.loadtxt(corr_file, unpack=True, usecols=[0,1,2,3]) 
             
            # get z values and weights
            z_dist = corr_z
            z_weights = corr_w 

        # upweighted corrected data
        elif corr.lower() == 'upweight': 
            if catalog.lower() == 'lasdamasgeo': 
                corr_file = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'+\
                        'sdssmock_gamma_lrgFull_zm_oriana'+str("%02d" % catalog_param[0])+catalog_param[1]+\
                        '_no.rdcz.fibcoll.dat'
            else: 
                raise NameError('Only LasDamasGeo coded so far') 
            
            # read-in corrected file
            corr_ra, corr_dec, corr_z, corr_w = np.loadtxt(corr_file, unpack=True, usecols=[0,1,2,3]) 
             
            # get z values and weights
            z_dist = corr_z
            z_weights = corr_w 
   
    nbar_ngal = np.zeros(len(z_values[0]))
    for i_z, zmid in enumerate(z_values[0]):
        zlim = (z_dist >= (z_values[1])[i_z]) & (z_dist < (z_values[2])[i_z])
        nbar_ngal[i_z] = np.sum(z_weights[zlim]) 
    
    # write nbar_ngal data to ask ascii file 
    nbar_ngal_file = get_nbar_ngal_file(DorR=DorR, catalog=catalog, catalog_param=catalog_param, 
            corr=corr, corr_param=corr_param)
    np.savetxt(nbar_ngal_file,
            np.c_[z_values[0], z_values[1], z_values[2], z_values[3], nbar_ngal], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def get_nbar_ngal_file(DorR='data', catalog='lasdamasgeo', catalog_param=[1,'a'], 
        corr='peaknbar', corr_param=[5.3, 0.1]): 
    '''
    get nbar_ngal file name 
    '''
    file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
    file_prefix = 'nbar-ngal-'

    # specify catalog
    catalog_str = catalog.lower()
    if DorR.lower() == 'data':
        catalog_str = catalog_str+'-'+str(catalog_param[0])+catalog_param[1]

    # specify data or random
    DorR_str = DorR.lower() 

    # specify correction 
    if corr.lower() == 'true': 
        corr_str = 'true'
    elif corr.lower() == 'peaknbar': 
        corr_str = 'peak.sigma'+str(corr_param[0])+'fpeak'+str(corr_param[1])
    elif corr.lower() == 'upweight': 
        corr_str = 'fibcoll.upweight' 
    
    file_name = file_dir+file_prefix+catalog_str+DorR_str+corr_str+'.dat'
    return file_name 
            
def nbar_comparison(catalog='lasdamasgeo', catalog_param_list=[[1,'a']], 
        corr='peaknbar', corr_param_list=[[5.3, 0.1]]): 
    '''
    nbar comparison (currently hardcoded for LasDamasGeo Peak+Nbar Fibercollision Correction 
    '''
    z_values = readnbarz(file="/mount/riachuelo1/hahn/data/nbar-junk.dat")                  # import z-values
   
    # configure plot
    prettyplot() 
    fig = plt.figure(1, figsize=(10,8)) 
    sub = fig.add_subplot(111) 

    # import true random (hardcoded for lasdamasgeo)  
    true_random_file = get_nbar_ngal_file(DorR='random', catalog=catalog, catalog_param=catalog_param_list[0], 
        corr='true')
    true_random_ngal = np.loadtxt(true_random_file, unpack=True, usecols=[4]) 
        
    # only redshift bins within 0.16 < z < 0.44  
    within_zlim = (z_values[0] >= 0.16) & (z_values[0] < 0.44)
    
    # Loop through true and corrected data catalog parameters
    for i_catalog, catalog_param in enumerate(catalog_param_list):  
        # import true data
        true_file = get_nbar_ngal_file(DorR='data', catalog=catalog, catalog_param=catalog_param, corr='true')

        if os.path.isfile(true_file) == False:  # if file doesn't exist
            print true_file, ' does not exist... will compute'
            write_nbar_ngal(DorR='data', catalog=catalog, catalog_param=catalog_param, 
                    corr='true')
        true_ngal = np.loadtxt(true_file, unpack=True, usecols=[4])
     
        # for the list of correction parameters get correction file
        for i_corr, corr_param in enumerate(corr_param_list): 
            # import corrected data
            corr_file = get_nbar_ngal_file(DorR='data', catalog=catalog, catalog_param=catalog_param, corr='peaknbar',
                    corr_param=corr_param)
            if os.path.isfile(corr_file) == False: 
                print corr_file, ' does not exist... will compute'
                write_nbar_ngal(DorR='data', catalog=catalog, catalog_param=catalog_param, corr='peaknbar', 
                        corr_param=corr_param) 
            corr_ngal = np.loadtxt(corr_file, unpack=True, usecols=[4])
            
            if (i_catalog == 0) and (i_corr == 0):
                sub.plot((z_values[0])[within_zlim], (corr_ngal[within_zlim])/(true_ngal[within_zlim]), 
                    color=(214/255., 39/255., 40/255.), lw=2, 
                    label=r'$\bar{n}_{corrected}(z)/\bar{n}_{true}(z)$') 
            else: 
                sub.plot((z_values[0])[within_zlim], (corr_ngal[within_zlim])/(true_ngal[within_zlim]), 
                    color=(214/255., 39/255., 40/255.), lw=2) 

            # import corrected random nbar
            corr_random_file = get_nbar_ngal_file(DorR='random', catalog=catalog, corr='peaknbar', corr_param=corr_param)
            if os.path.isfile(corr_random_file) == False: 
                print corr_random_file, ' does not exist... will compute'
                write_nbar_ngal(DorR='random', catalog=catalog, corr='peaknbar', corr_param=corr_param) 
            corr_random_ngal = np.loadtxt(corr_random_file, unpack=True, usecols=[4])
    
    # Loop through corrected random catalogs
    for i_corr, corr_param in enumerate(corr_param_list): 
        sub.plot((z_values[0])[within_zlim], (corr_random_ngal[within_zlim])/(true_random_ngal[within_zlim]), 
                color='black', lw=2, 
                label=r'$\bar{n}_{corrected\;random}(z)/\bar{n}_{random}(z)$' ) 

    # configure axis
    sub.set_ylabel(r'$\bar{n}(z)$ ratio') 
    sub.set_xlim([0.16, 0.44]) 
    sub.set_ylim([0.9, 1.3])
    sub.set_xlabel('z') 
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/', 
        'nbar_comp_', catalog, '_data_random_comparison.png']), bbox_inches="tight")
    fig.clear() 

def nbar_comparison_data_over_random(catalog='lasdamasgeo', catalog_param_list=[[1,'a']], 
        corr='peaknbar', corr_param=[5.3, 0.1]): 
    '''
    comparison of n_g/n_r for true and corrected data/random combination  
    '''
    z_values = readnbarz(file="/mount/riachuelo1/hahn/data/nbar-junk.dat")                  # import z-values
    
    # only redshift bins within 0.16 < z < 0.44  
    within_zlim = (z_values[0] >= 0.16) & (z_values[0] < 0.44)
   
    # configure plot
    prettyplot() 
    fig = plt.figure(1, figsize=(10,8)) 
    sub = fig.add_subplot(111) 

    # import true random (hardcoded for lasdamasgeo)  
    true_random_file = get_nbar_ngal_file(DorR='random', catalog=catalog, catalog_param=catalog_param_list[0], 
        corr='true')
    true_random_ngal = np.loadtxt(true_random_file, unpack=True, usecols=[4]) 
    N_rand_true = np.sum(true_random_ngal) 
    
    # import corrected random nbar
    corr_random_file = get_nbar_ngal_file(DorR='random', catalog=catalog, corr='peaknbar', corr_param=corr_param)
    if os.path.isfile(corr_random_file) == False: 
        print corr_random_file, ' does not exist... will compute'
        write_nbar_ngal(DorR='random', catalog=catalog, corr='peaknbar', corr_param=corr_param) 
    corr_random_ngal = np.loadtxt(corr_random_file, unpack=True, usecols=[4])
    N_rand_corr = np.sum(corr_random_ngal) 

    print 'Nrandom True = ', N_rand_true, ' Nrandom Corrected = ', N_rand_corr           # should be the same 

    # Loop through true and corrected data catalog parameters
    sum_true_ngal = np.zeros(len(z_values[0]))
    sum_corr_ngal = np.zeros(len(z_values[0]))
    for i_catalog, catalog_param in enumerate(catalog_param_list):  
        # import true data
        true_file = get_nbar_ngal_file(DorR='data', catalog=catalog, catalog_param=catalog_param, corr='true')

        if os.path.isfile(true_file) == False:  # if file doesn't exist
            print true_file, ' does not exist... will compute'
            write_nbar_ngal(DorR='data', catalog=catalog, catalog_param=catalog_param, 
                    corr='true')
        true_ngal = np.loadtxt(true_file, unpack=True, usecols=[4])

        sum_true_ngal = sum_true_ngal+true_ngal
     
        # import corrected data
        corr_file = get_nbar_ngal_file(DorR='data', catalog=catalog, catalog_param=catalog_param, corr='peaknbar',
                corr_param=corr_param)
        if os.path.isfile(corr_file) == False: 
            print corr_file, ' does not exist... will compute'
            write_nbar_ngal(DorR='data', catalog=catalog, catalog_param=catalog_param, corr='peaknbar', 
                    corr_param=corr_param) 
        corr_ngal = np.loadtxt(corr_file, unpack=True, usecols=[4])
        
        sum_corr_ngal = sum_corr_ngal+corr_ngal
           
    # need to scale the data ngal values 
    sum_sum_true = np.sum(sum_true_ngal) 
    sum_sum_corr = np.sum(sum_corr_ngal) 
    alpha_true = sum_sum_true/N_rand_true
    alpha_corr = sum_sum_corr/N_rand_corr 
            
    sub.plot((z_values[0])[within_zlim], (alpha_corr*corr_random_ngal[within_zlim])/(sum_corr_ngal[within_zlim]), 
            color=(214/255., 39/255., 40/255.), lw=2, 
            label=r'$\bar{n}_{g,corr}(z)/\bar{n}_{r,corr}(z)$') 

    sub.plot((z_values[0])[within_zlim], (alpha_true*true_random_ngal[within_zlim])/(sum_true_ngal[within_zlim]), 
            color='black', lw=2, 
            label=r'$\bar{n}_{g,true}(z)/\bar{n}_{r,true}(z)$' ) 

    # configure axis
    sub.set_ylabel(r'$\bar{n}(z)$ ratio') 
    sub.set_xlim([0.16, 0.44]) 
    sub.set_ylim([0.9, 1.3])
    sub.set_xlabel('z') 
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/', 
        'nbar_comp_', catalog, '_data_over_random_comparison.png']), bbox_inches="tight")
    fig.clear() 

def nbar_comparison_data_corr(catalog='lasdamasgeo', catalog_param_list=[[1,'a']], 
        corr='peaknbar', corr_param=[5.3, 0.1]): 
    '''
    comparison of n_g/n_r for true and corrected data/random combination  
    '''
    z_values = readnbarz(file="/mount/riachuelo1/hahn/data/nbar-junk.dat")                  # import z-values
    
    # only redshift bins within 0.16 < z < 0.44  
    within_zlim = (z_values[0] >= 0.16) & (z_values[0] < 0.44)
   
    # configure plot
    prettyplot() 
    fig = plt.figure(1, figsize=(10,8)) 
    sub = fig.add_subplot(111) 

    # Loop through true and corrected data catalog parameters
    sum_true_ngal = np.zeros(len(z_values[0]))
    sum_corr_ngal = np.zeros(len(z_values[0]))
    for i_catalog, catalog_param in enumerate(catalog_param_list):  
        # import true data
        true_file = get_nbar_ngal_file(DorR='data', catalog=catalog, catalog_param=catalog_param, corr='true')

        if os.path.isfile(true_file) == False:  # if file doesn't exist
            print true_file, ' does not exist... will compute'
            write_nbar_ngal(DorR='data', catalog=catalog, catalog_param=catalog_param, 
                    corr='true')
        true_ngal = np.loadtxt(true_file, unpack=True, usecols=[4])

        sum_true_ngal = sum_true_ngal+true_ngal
     
        # import corrected data
        corr_file = get_nbar_ngal_file(DorR='data', catalog=catalog, catalog_param=catalog_param, corr=corr,
                corr_param=corr_param)
        print corr_file

        if os.path.isfile(corr_file) == False: 
            print corr_file, ' does not exist... will compute'
            write_nbar_ngal(DorR='data', catalog=catalog, catalog_param=catalog_param, corr=corr, 
                    corr_param=corr_param) 
        # get corr ngal 
        corr_ngal = np.loadtxt(corr_file, unpack=True, usecols=[4])
        
        sum_corr_ngal = sum_corr_ngal+corr_ngal
   
    if corr == 'peaknbar': 
        corr_param_str = ''.join([str(i) for i in corr_param])
    elif corr == 'upweight': 
        corr_param_str = ''

    sub.plot((z_values[0])[within_zlim], (sum_corr_ngal[within_zlim])/(sum_true_ngal[within_zlim]), 
            color=(214/255., 39/255., 40/255.), lw=2, 
            label=r'$\bar{n}_{g,'+corr+corr_param_str+r'}(z)/\bar{n}_{g,true}(z)$') 

    # configure axis
    sub.set_ylabel(r'$\bar{n}(z)$ ratio') 
    sub.set_xlim([0.16, 0.44]) 
    sub.set_ylim([0.9, 1.3])
    sub.set_xlabel('z') 
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/', 
        'nbar_comp_', catalog, '_data_', corr, corr_param_str, '_corrected_comparison.png']), bbox_inches="tight")
    fig.clear() 

def plot_upweighted_gal_redshift():
    '''
    plot the redshifts of the upweighted galaxies 
    everything hardcoded for now 
    '''
    ldg_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
    ldg_fibcoll_file = 'sdssmock_gamma_lrgFull_zm_oriana01a_no.rdcz.fibcoll.dat'

    ra, dec, z, wcp = np.loadtxt(ldg_dir+ldg_fibcoll_file, unpack=True, usecols=[0,1,2,3])

    fig = plt.figure(1) 
    sub = fig.add_subplot(111)
    sub.hist(z[wcp > 1], 20, normed=1, histtype='bar') 
    sub.set_xlim([0.16, 0.44])

    plt.show() 


if __name__=='__main__': 
    '''
    for DorR in ['data']: #, 'random']: 
        for param in [[i,j] for i in range(1,6) for j in ['a', 'b', 'c', 'd']]: 
            write_nbar_ngal(DorR=DorR, catalog='lasdamasgeo', catalog_param=param, 
                    corr='upweight')#, corr_param=[5.3,0.1])
    '''
    ldg_catalog_list = [[i, j] for i in range(1,6) for j in ['a', 'b', 'c', 'd']]
    peakcorr_param_list = [[5.3, i*0.1] for i in range(11)]
    #nbar_comparison(catalog='lasdamasgeo', catalog_param_list=ldg_catalog_list) 

    #nbar_comparison_data_over_random(catalog='lasdamasgeo', catalog_param_list=ldg_catalog_list, 
    #        corr='peaknbar', corr_param=[5.3, 0.1])
    #nbar_comparison_data_corr(catalog='lasdamasgeo', catalog_param_list=ldg_catalog_list, 
    #    corr='upweight')
    #nbar_comparison_data_corr(catalog='lasdamasgeo', catalog_param_list=ldg_catalog_list, 
    #    corr='peaknbar', corr_param=[5.3, 0.0]) 
    plot_upweighted_gal_redshift()
