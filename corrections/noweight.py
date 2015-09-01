

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
