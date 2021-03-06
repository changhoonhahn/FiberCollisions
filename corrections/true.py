'''

"True" galaxy catalogs 

Everything is hardcoded. Code can be improved but too lazy. 


'''
import numpy as np
# --- Local ---
from corrections import Corrections

class TrueCorr(Corrections): 

    def __init__(self, cat_corr, **kwargs): 
        ''' Child class of the Corrections class in corrections.py
        '''
        super(TrueCorr, self).__init__(cat_corr, **kwargs)
        
        self.corrstr() 
    
    def corrstr(self): 
        ''' Correction string for true is nothing
        '''
        self.corr_str = ''
        return self.corr_str

    def build(self): 
        """ Construct "true" galaxy catalogs with no systematic effects 

        Notes
        -----
        * Tiling Mock : 
            Reads in true mock catalogs and imposes predetermined redshift 
            limits on them and attaches a flag mostly hardcoded since it's a 
            simple procedure 
        * QPM: 
            Handles messy weights 
        * PATCHY: Returns mocks with only necessary columns and w_fc = 1
        * Everything is very hardcoded
        """
        catalog = (self.cat_corr)['catalog']

        data_cols = self.datacolumns()
        data_fmts = self.datacols_fmt()
        data_hdrs = self.datacols_header()

        output_file = self.file()
        
        if catalog['name'].lower() == 'nseries':       # N Series
            # read in original files from Jeremy and adjust them to make them
            # easier to use for fiber collisions
            # read rdzw file 
            data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
            orig_file = ''.join([
                data_dir, 
                'CutskyN', str(catalog['n_mock']), '.rdzwc'
                ]) 
            orig_ra, orig_dec, orig_z, orig_wfc, orig_zupw, orig_upw_index = np.loadtxt(
                    orig_file, 
                    unpack = True, 
                    usecols=[0,1,2,4,5,6]
                    )
        
            # file with completeness
            mask_file = ''.join([
                data_dir, 
                'CutskyN', str(catalog['n_mock']), '.mask_info'
                ]) 

            orig_wcomp = np.loadtxt(
                    mask_file, 
                    unpack = True, 
                    usecols = [0]
                    ) 

            # true wfc 
            true_wfc = np.array([ 1.0 for i in range(len(orig_wfc)) ]) 
            
            # write to file 
            data_list = [orig_ra, orig_dec, orig_z, true_wfc, orig_wcomp, orig_zupw, orig_upw_index]

        elif catalog['name'].lower() == 'qpm': 
            # import original true data 
            orig_file = ''.join([
                '/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
                'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 
            orig_data = np.loadtxt(orig_file) 

            orig_info_file = orig_file+'.info'
            orig_info = np.loadtxt(orig_info_file)    # gal_id, comp, z_real, z_red, mass_halo, flag_sta, id_halo
        
            # consistency issue with #46
            if catalog['n_mock'] in (46, 52, 53, 54, 56, 61, 707, 756, 794, 819, 831, 835, 838): 
                orig_veto_file = ''.join(['/mount/riachuelo1/hahn/data/QPM/dr12d/', 
                    'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
            else:
                orig_veto_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
                    'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
            orig_veto = np.loadtxt(orig_veto_file) 
            n_gal = len(orig_veto)
       
            # assign RA, Dec, z, and w_veto 
            true_ra     = orig_data[:,0]
            true_dec    = orig_data[:,1]
            true_z      = orig_data[:,2]
            true_wfkp   = orig_data[:,3]
            true_wfc = np.repeat(1.0, n_gal)    # fiber collisions weights are all 1 for true

            # check to make sure that the redshifts correspond btw rdz file and rdz.info file 
            if np.max(np.abs(orig_info[:,3]/true_z-1.0)) > 10**-5: 
                raise ValueError('redshifts between the data file and info file dont match') 

            true_comp = orig_info[:,1]         # completness weights

            # remove veto mask 
            # Only keep galaxies with veto = 0 (for veto values in .veto file) 
            vetomask = np.where(orig_veto == 0)            
            data_list = [
                    true_ra[vetomask], true_dec[vetomask], true_z[vetomask], 
                    true_wfc[vetomask], true_comp[vetomask]
                    ]

        elif catalog['name'].lower() == 'bigmd':                
            # Big MultiDark  
            P0 = 20000.0
            # read original file
            data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
            if 'version' in catalog.keys():
                if catalog['version'] == 'nowiggles':
                    # simulation with no wiggles 
                    orig_file = ''.join([data_dir, 'nowiggles/BigMD-cmass-dr12v4-nowiggle-veto.dat']) 
                else: 
                    raise NotImplementedError
                    #orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-wcp-veto.dat'])  # hardcoded
                    #orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-standardHAM-veto.dat'])
                    #orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-quadru-veto.dat'])
            else:       # default 
                orig_file = ''.join([data_dir, 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak-veto.dat']) 

            orig_ra, orig_dec, orig_z, orig_wfkp, orig_veto, orig_wfc = np.loadtxt(
                    orig_file, 
                    unpack=True, 
                    usecols=[0,1,2,3,4,5])

            # true wfc = 1 for all galaxies 
            true_wfc = np.repeat(1., len(orig_wfc))
            #nbar = (1.0/P0) * (1.0/orig_wfkp - 1.0) 

            vetomask = np.where(orig_veto == 1)     # if veto = 1 then keep; otherwise discard 
            data_list = [
                    orig_ra[vetomask], 
                    orig_dec[vetomask], 
                    orig_z[vetomask], 
                    true_wfc[vetomask]
                    ]

        elif catalog['name'].lower() == 'qso_bigmd':
            # QSO Big MultiDark  
            P0 = 20000.0
            # read original file
            data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
            if catalog['version'] == 'evo':
                # simulation with no wiggles 
                orig_file = ''.join([data_dir, 
                    'QSO-BigMD-1481.43deg2-bin0_v1.0_evo.dat'])
            elif catalog['version'] == 'noevo':
                orig_file = ''.join([data_dir, 
                    'QSO-BigMD-1481.43deg2-bin0_v1.0_noevo.dat'])
            elif catalog['version'] == 'eboss': 
                orig_file = ''.join([data_dir, 
                    'eboss_v1.0-QSO-NS-eboss_v1.0_bin0.dat'])
            elif catalog['version'] == 'ebossv1.5': # eBOSS v1.5 
                P0 = 6000.  
                orig_file = ''.join([data_dir, 
                     'eboss_v1.5-QSO-eboss_v1.5.dat'])
            elif catalog['version'] == 'ebossnew':
                orig_file = ''.join([data_dir, 
                    'QSOs_y1_comp_cut_0.5_double_weight_col_Z.dat']) 
            elif 'jackknife' in catalog['version']: 
                n_jack = catalog['version'].split('jackknife')[-1]
                orig_file = ''.join([data_dir, 
                    'QSO-bin0_', str(n_jack), '.dat']) 
            elif 'v2' in catalog['version']: 
                P0 = 6000.  
                if 'z' in catalog['version']: 
                    orig_file = ''.join([data_dir, 'BigMDPL-QSOZ.dat']) 
                elif 'nsat' in catalog['version'] : 
                    orig_file = ''.join([data_dir, 'BigMDPL-QSO-NSAT.dat']) 
                else: 
                    orig_file = ''.join([data_dir, 'BigMDPL-QSO.dat']) 

            elif catalog['version'] == 'other_mock':    
                # other mocks that Sergio sent me Nov 1, 2016
                orig_file = ''.join([data_dir, 
                    'mockPM_QSO_', str("%04d" % catalog['n_mock']), '.dat'])
            else: 
                raise NotImplementedError

            orig_ra, orig_dec, orig_z, orig_wfkp = np.loadtxt(
                    orig_file, 
                    unpack=True, 
                    usecols=[0,1,2,3])

            orig_nbar = (1.0/P0) * (1.0/orig_wfkp - 1.0) 

            # true wfc = 1 for all galaxies 
            true_wfc = np.repeat(1., len(orig_z))
            
            if catalog['version'] == 'ebossnew': 
                wfkpzero = np.where(orig_wfkp != 0.)  
                data_list = [
                        orig_ra[wfkpzero], 
                        orig_dec[wfkpzero], 
                        orig_z[wfkpzero], 
                        orig_nbar[wfkpzero], 
                        true_wfc[wfkpzero]
                        ]
            else: 
                data_list = [orig_ra, orig_dec, orig_z, orig_nbar, true_wfc]


        elif catalog['name'].lower() == 'tilingmock':   # tiling mock 
            # read in original file and impose redshift limits 
            orig_ra, orig_dec, orig_z, orig_w = np.loadtxt(
                    ''.join(['/mount/riachuelo1/hahn/data/tiling_mocks/', 
                        'cmass-boss5003sector-icoll012.dat']),
                    unpack=True, 
                    usecols=[0,1,2,3]) 
            
            zlow, zhigh = self.survey_zlimits()
            zlim = np.where((orig_z > zlow) & (orig_z < zhigh))
            
            true_wfc = np.repeat(1.0, len(orig_ra))
            data_list = [orig_ra[zlim], orig_dec[zlim], orig_z[zlim], true_wfc[zlim]]

        else: 
            raise NameError('not yet coded') 
            
        # write to file 
        np.savetxt(
                self.file(), 
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmts, 
                delimiter='\t', 
                header=data_hdrs
                ) 

        return None 

"""
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
        '''

        if catalog['name'].lower() == 'qpm':          # QPM 

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
    
    #DEPRECIATED
            elif catalog['name'].lower() == 'patchy':       # PATCHY mocks ------------------------

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
               
                header_str = "Columns : ra, dec, z, nz, wfc" 
                data_list = [orig_ra[vetomask], orig_dec[vetomask], orig_z[vetomask], 
                            orig_nbar[vetomask], new_wfc[vetomask]] 
                data_fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f']
        
            elif catalog['name'].lower() == 'lasdamasgeo':          
                
                orig_true_file = ''.join(['/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/', 
                    'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.dat']) 
                orig_true_data = np.loadtxt(orig_true_file, unpack=True, usecols=[0,1,2])         # ra, dec, ***CZ***

                true_ra = orig_true_data[0]
                true_dec = orig_true_data[1]
                true_z = orig_true_data[2]/299800.0         # convert cz to z

                true_weight = np.array([1.0 for j in range(len(true_z))])   # no weights for true (all 1) 

                header_str = "Columns : ra, dec, z, weight"
                data_list = [true_ra, true_dec, true_z, true_weight]
                data_fmt = ['%10.5f', '%10.5f', '%10.5f', '%10.5f']

"""
