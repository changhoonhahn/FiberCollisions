

def build_random(**cat_corr): 
    ''' Build the random catalogs from original data 

    Paramter
    --------
    cat_corr : Catalog and Correction dictionary
    '''
    catalog = cat_corr['catalog']
    if 'cmass' in catalog['name'].lower():          # CMASS -------------------------------- 
        data_dir = '/mount/riachuelo1/hahn/data/CMASS/'
        if 'cmasslowz' in catalog['name'].lower(): 
            data_dir += 'dr12v5/'

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

        elif 'cmasslowz' in catalog['name'].lower():   
            # CMASS LOWZ combined data
            
            # three different CMASS LOWZ  
            cmasslowz_str = ''
            if 'e2' in catalog['name'].lower(): 
                cmasslowz_str = 'E2' 
            elif 'e3' in catalog['name'].lower(): 
                cmasslowz_str = 'E3'
            elif 'tot' in catalog['name'].lower(): 
                cmasslowz_str = 'TOT'

            if 'high' in catalog['name'].lower(): 
                zmin, zmax = 0.5, 0.75
            elif '_low' in catalog['name'].lower():
                zmin, zmax = 0.2, 0.5
            else: 
                raise NameError("CMASSLOWZ Catalog must specify high or lowr edshift bin") 
            
            if 'tot' not in catalog['name'].lower(): 
                # random data fits file
                data_file = ''.join([data_dir, 'random0_DR12v5_CMASSLOWZ', cmasslowz_str, '_North.fits.gz'])
                # old version 'cmasslowz-dr12v4-N-Reid.ran.fits'
                cmass = mrdfits(data_file) 
            
                # mask file 
                mask_file = ''.join([data_dir, 'mask_DR12v5_CMASSLOWZ', cmasslowz_str, '_North.fits.gz'])
                mask = mrdfits(mask_file) 
                ipoly = cmass.ipoly # polygon index
                comp = mask.weight[ipoly]
            
                # redshift limit 
                zlimit = np.where((cmass.z >= zmin) & (cmass.z < zmax))
            else: 
                if 'high' in catalog['name'].lower(): 
                    hl_str = '_high'
                elif '_low' in catalog['name'].lower():
                    hl_str = '_low'

                cc_cmd = 'cat ' 
                for comb_str in ['', 'e2', 'e3']: 
                    cc = {'catalog': {'name': 'cmasslowz'+hl_str+comb_str}, 
                            'correction': {'name': 'upweight'}}
                    cc_file = get_galaxy_data_file('random', **cc)
                    cc_cmd += cc_file+' ' 

                random_file = get_galaxy_data_file('random', **cat_corr) 
                cc_cmd += '> '+random_file 
                return 

        else: 
            raise NotImplementedError("Only CMASS and CMASS+LOWZ combined sample implemented") 
    
        head_str = 'columns : ra, dec, z, nbar, comp'
        #ra, dec, z, nz, comp 
        random_file = get_galaxy_data_file('random', **cat_corr) 
        np.savetxt(random_file, 
                np.c_[
                    (cmass.ra)[zlimit], (cmass.dec)[zlimit], (cmass.z)[zlimit], 
                    (cmass.nz)[zlimit], comp[zlimit]], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f'], delimiter='\t', header=head_str) 

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
