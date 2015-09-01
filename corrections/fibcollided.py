'''

Galaxy catalogs with fiber collisions 

Everything is hardcoded. Code can be improved but too lazy. 


'''
import correction as corrclass


def file(cat_corr): 
    """ Specify correction string
    """
    cat = cat_corr['catalog']

    if 'cmass' in cat['name'].lower(): 
        return ''
    
    corr_str = '.fibcoll'

    return corr_str 

def build(cat_corr): 
    ''' Build Fibercollided mock catalogs using specific idl routines or by using the given fiber collision weights

    Parameters
    ----------
    cat_corr : catalog correction dictionary 

    Notes
    -----
    '''
    catalog = cat_corr['catalog']

    if catalog['name'].lower() == 'nseries':          # N-series 

        # original file 
        data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
        orig_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.rdzwc']) 
        orig_ra, orig_dec, orig_z, orig_wfc = np.loadtxt(orig_file, unpack=True, usecols=[0,1,2,4])
    
        # file with mask completeness
        mask_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.mask_info']) 
        orig_wcomp = np.loadtxt(mask_file, unpack=True, usecols=[0]) 

        header_str = 'Columns : ra, dec, z, nbar, w_cp, comp' 
        data_list = [orig_ra, orig_dec, orig_z, orig_wfc, orig_wcomp]   # data column list 
        data_fmt = ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f']   # data format

    else: 
        raise NotImplementedError('not yet coded') 
    
    # corrected file 
    corr_class = corrclass.correction(cat_corr) 
    output_file = corr_class.file()
    
    np.savetxt(output_file, (np.vstack(np.array(data_list))).T, fmt=data_fmt, delimiter='\t', header=header_str) 

    return None

    """
    if 'cmass' in catalog['name'].lower():                  # CMASS ---------------------

        if catalog['name'].lower() == 'cmass': 
            data_file = ''.join([data_dir, 'cmass-dr12v4-N-Reid.dat.fits']) # fits file 
            data = mrdfits(data_file)

            zlimit = np.where((data.z >= 0.43) & (data.z <= 0.7))
        '''
        elif 'cmasslowz' in catalog['name'].lower(): 
            # CMASS LOWZ combined sample
            cmasslowz_str = ''
            if 'e2' in catalog['name'].lower(): 
                cmasslowz_str = 'E2'
            elif 'e3' in catalog['name'].lower(): 
                cmasslowz_str = 'E3'
            elif 'tot' in catalog['name'].lower(): 
                cmasslowz_str = 'TOT'

            if '_low' in catalog['name'].lower(): 
                zmin, zmax = 0.2, 0.5
            elif 'high' in catalog['name'].lower(): 
                zmin, zmax = 0.5, 0.75
            else: 
                raise NameError("redshift bin must be specified") 
        
            if 'tot' not in catalog['name'].lower(): 
                # original combined data sample
                data_file = ''.join([data_dir, 'galaxy_DR12v5_CMASSLOWZ', cmasslowz_str, '_North.fits.gz'])
                data = mrdfits(data_file) 

                zlimit = np.where((data.z >= zmin) & (data.z < zmax))  # redshift limit
            else: 
                # Concatenate the other CMASSLOWZ 
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
        ''' 
        head_str = 'columns : ra, dec, z, nbar, w_systot, w_noz, w_cp, comp'

        # save to file 
        np.savetxt(output_file, 
                np.c_[
                    (data.ra)[zlimit], (data.dec)[zlimit], (data.z)[zlimit], (data.nz)[zlimit],
                    (data.weight_systot)[zlimit], (data.weight_noz)[zlimit], 
                    (data.weight_cp)[zlimit], (data.comp)[zlimit]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', 
                    '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t', 
                header=head_str) 

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
        
        np.savetxt(output_file, 
                np.c_[
                    orig_ra[vetomask], orig_dec[vetomask], orig_z[vetomask], 
                    orig_nbar[vetomask], orig_wfc[vetomask]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f'], 
                delimiter='\t') 
    
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
    """
