'''

"True" random data for simulated/observed catalogs 

'''

import numpy as np
# --- Local ---
from defutility.fitstables import mrdfits
from corrections import Corrections
from util.direc import direc

class Rand(Corrections): 

    def __init__(self, cat_corr, **kwargs): 
        """ Child clas of Correction for random data of 
        simulation/data catalogs
        """

        super(Rand, self).__init__(cat_corr, **kwargs)

        self.corrstr() 

    def file(self): 
        """ Override parent class file attribute 
        Random catalog file 
        """

        cat_name = ((self.cat_corr)['catalog'])['name'].lower()
    
        data_dir = direc('data', self.cat_corr)

        if 'cmass' in cat_name: # CMASS

            if cat_name == 'cmass': 
                # CMASS random catalog 
                file_name = 'cmass-dr12v4-N-Reid.ran.dat'

            elif 'cmasslowz' in cat_name:  

                # CMASS LOWZ combined random catalog
                if 'e2' in cat_name: 
                    cmasslowz_str = 'e2'
                elif 'e3' in cat_name: 
                    cmasslowz_str = 'e3'
                else: 
                    cmasslowz_str = ''
                        
                if '_low' in cat_name: 
                    zbin_str = '_LOW' 
                elif 'high' in cat_name: 
                    zbin_str = '_HIGH'
                else: 
                    raise NameError("Must specify redshift bin of CMASS LOWZ sample") 

                file_name = ''.join([
                    'random0_DR12v5_CMASSLOWZ', 
                    cmasslowz_str.upper(), 
                    zbin_str, 
                    '_North.ran.dat'
                    ])
            else: 
                raise NotImplementedError()

        elif cat_name == 'nseries':   # Nseries

            file_name = 'Nseries_cutsky_randoms_50x_redshifts_comp.dat'

        else: 
            raise NotImplementedError()
        
        return  ''.join([data_dir, file_name])

    def build(self): 
        ''' Build the random catalogs from original data 
        Override parent class build attribute 

        '''
        cat_name = ((self.cat_corr)['catalog'])['name'].lower()

        if cat_name == 'nseries':     # Nseries ----------------------------

            # original random catalog 
            data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
            orig_rand_file = ''.join([data_dir, 'Nseries_cutsky_randoms_50x_redshifts.dat']) 
            ra, dec, z = np.loadtxt(orig_rand_file, unpack=True, usecols=[0,1,2]) # RA, Decl, Redhsift
        
            # sector completeness catalog
            orig_comp_file = ''.join([data_dir, 'Nseries_cutsky_randoms_50x_maskinfo.dat'])  
            comp = np.loadtxt(orig_comp_file, unpack=True, usecols=[0])
            
            header_str = 'Columns : ra, dec, z, comp'
            data_list = [ra, dec, z, comp]
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f']

        elif 'cmass' in cat_name:          # CMASS -------------------------------- 

            data_dir = direc('data', self.cat_corr) 

            if cat_name == 'cmass': 
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

            elif 'cmasslowz' in cat_name:   
                # CMASS LOWZ combined data
                
                # three different CMASS LOWZ  
                if 'e2' in cat_name: 
                    cmasslowz_str = 'E2' 
                elif 'e3' in cat_name: 
                    cmasslowz_str = 'E3'
                else: 
                    cmasslowz_str = ''

                if 'high' in cat_name: 
                    zmin, zmax = 0.5, 0.75
                elif '_low' in cat_name:
                    zmin, zmax = 0.2, 0.5
                else: 
                    raise NameError("CMASSLOWZ Catalog must specify high or lowr edshift bin") 
                
                # random data fits file
                data_file = ''.join([
                    data_dir, 
                    'random0_DR12v5_CMASSLOWZ', 
                    cmasslowz_str, 
                    '_North.fits.gz'
                    ])
                cmass = mrdfits(data_file) 
            
                # mask file 
                mask_file = ''.join([
                    data_dir, 
                    'mask_DR12v5_CMASSLOWZ', 
                    cmasslowz_str, 
                    '_North.fits.gz'
                    ])
                mask = mrdfits(mask_file) 
                ipoly = cmass.ipoly # polygon index
                comp = mask.weight[ipoly]
            
                # redshift limit 
                zlimit = np.where((cmass.z >= zmin) & (cmass.z < zmax))

            else: 
                raise NotImplementedError("Only CMASS and CMASS+LOWZ combined sample implemented") 
        
            header_str = 'columns : ra, dec, z, nbar, comp'  #ra, dec, z, nz, comp 
            data_list = [(cmass.ra)[zlimit], (cmass.dec)[zlimit], (cmass.z)[zlimit], (cmass.nz)[zlimit], comp[zlimit]]
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f']
        
        else:
            raise NotImplementedError()

        # write to corrected file 
        output_file = self.file(cat_corr, **kwargs)
        np.savetxt(
                output_file, 
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmt, 
                delimiter='\t', 
                header=header_str
                ) 

        return None 

"""
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
"""
