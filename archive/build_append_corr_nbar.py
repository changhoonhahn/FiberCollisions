import fibercollisions as fc
import sys 

def append_corr_nbar(DorR='data', catalog='LasDamasGeo', catalog_param={'n_mock':[1,'a']}, 
        corr='peaknbar', corr_param={'sigma':5.3, 'fpeak':0.5}, sanitycheck=False): 
    '''
    append corrected interpolated nbar(z) to corrected data or random file 
    '''
    # read in data/random galaxies 
    gal_data = fc.galaxy_data(DorR=DorR, catalog=catalog, catalog_param=catalog_param, corr=corr, corr_param=corr_param)
    if DorR == 'random': 
        print 'Nran = ', len(gal_data.ra) 
    # read in corrected nbar file 
    corr_nbar = fc.nbar(catalog=catalog, corr=corr, corr_param=corr_param)
    
    # interpolate within redshift limits 
    if catalog.lower() == 'lasdamasgeo': 
        zlim = (corr_nbar.zmid > 0.16) & (corr_nbar.zmid < 0.44)        # for las damas geo 
    else: 
        raise NameError("not yet coded!")
    # numpy interpolate
    nbar_arr = np.interp(gal_data.z, corr_nbar.zmid[zlim], corr_nbar.nbar[zlim])
    
    if DorR == 'random': 
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
    if DorR == 'data': 
        np.savetxt(gal_corr_nbar_file, 
                np.c_[gal_data.ra, gal_data.dec, gal_data.z, nbar_arr, gal_data.weight],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e', '%10.5f'], delimiter='\t')
    elif DorR == 'random': 
        np.savetxt(gal_corr_nbar_file, 
                np.c_[gal_data.ra, gal_data.dec, gal_data.z, nbar_arr],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e'], delimiter='\t')

if __name__=='__main__': 
    print "RUNNING!"
    # hardcoded for las damas geo 
    dorr = sys.argv[1] 
    catalog = sys.argv[2]

    if catalog.lower() == 'lasdamasgeo': 
        catalog_param_0 = sys.argv[3]
        catalog_param_1 = sys.argv[4] 
        
        corr = sys.argv[5]
        
        if (corr.lower() == 'peaknbar') or (corr.lower() == 'peak'): 
            corr_sigma = sys.argv[6]
            corr_fpeak = sys.argv[7]

            print 'runninnnnnggg'
            append_corr_nbar(DorR=dorr, catalog=catalog, catalog_param={'n_mock':[np.int(catalog_param_0), catalog_param_1]}, 
                    corr=corr, corr_param={'sigma':np.float(corr_sigma), 'fpeak':np.float(corr_fpeak)}, sanitycheck=True)
