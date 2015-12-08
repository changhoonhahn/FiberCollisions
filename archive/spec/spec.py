'''

Spectrum class of FiberCollisions project

'''

import numpy as np
import os.path
import subprocess
import cosmolopy as cosmos

from data import Data 
from fft import Fft
from util.direc import direc
from util.fortran import Fcode
from fourier_corr.pk_corr import fourier_tophat_Pk

# Classes ------------------------------------------------------------
class Spec(object): 

    def __init__(self, spectype, cat_corr, ell=None, **kwargs):
        """ Class that describes power/bispectrum measurements 
        
        specify catalog, version, mock file number, file specifications (e.g. Nrandom), 
        fiber collision correction method, correction specifications (e.g. sigma, fpeak)

        Parameters 
        ----------
        spectype : 'power' or 'bispec'
        cat_corr : catalog and correction dictionary 


        """

        if spectype not in ['pk', 'bk']: 
            raise ValueError()
        else: 
            self.type = spectype

        if 'spec' not in cat_corr.keys(): 

            if ell is None: 
                raise ValueError

            # default spectrum parameters
            cat_corr['spec'] = {
                    'P0': 20000, #P0 
                    'Lbox': 3600, 
                    'Ngrid':360, 
                    'ell': ell 
                    }
            self.ell = ell

            if 'Ngrid' in kwargs.keys(): 
                cat_corr['spec']['Ngrid'] = kwargs.pop('Ngrid')
        else: 

            if 'ell' not in cat_corr['spec'].keys():
                raise ValueError

            if ell is not None: 
                if ell != cat_corr['spec']['ell']: 
                    raise ValueError

            self.ell = cat_corr['spec']['ell']
        
        self.cat_corr = cat_corr.copy()
        self.kwargs = kwargs

        self.file_name = self.file()
    
    def read(self): 
        """ Read power/bispectrum of simulated/observed data catalog
        """
    
        if self.ell == 0:   # monopole
                
            col_index = [0, 1]
            data_cols = ['k', 'p0k']

        elif self.ell == 2:     # quadrupoel
            
            col_index = [0, 2, 1, 3]
            data_cols = ['k', 'p2k', 'p0k', 'p4k']

        elif self.ell == 4:     # hexadecapole
            
            col_index = [0, 3, 1, 2]
            data_cols = ['k', 'p4k', 'p0k', 'p2k']

        else: 
            raise NotImplementedError()

        spec_data = np.loadtxt(
                    self.file_name, 
                    unpack = True, 
                    usecols = col_index
                    )

        for i_col, col in enumerate(data_cols): 
            setattr(self, col, spec_data[i_col])

        return None 

    def file(self):
        """ power/bispectrum file 
        """

        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']
        specdict = (self.cat_corr)['spec']

        # powerspectrum or bispectrum 
        if self.type in ('pk'): 
            spec_str = 'POWER_'
        elif self.type == 'bk':
            spec_str = 'BISP_'

        #if 'quad' not in specdict.keys(): 
        #    specdict['quad'] = False
        
        #if specdict['quad']:          
        #    spec_str += 'Q_'
        if self.ell != 0: 
            spec_str += 'Q_'

        gal_data = Data('data', self.cat_corr, **self.kwargs)
        self.data_file = gal_data.file_name
        gal_file = (gal_data.file_name).split('/')[-1]

        rand_data = Data('random', self.cat_corr, **self.kwargs)
        self.random_file = rand_data.file_name

        spec_dir = direc('spec', self.cat_corr)

        if self.type == 'pk': 
            specparam_str = ''.join([
                '.grid', str(specdict['Ngrid']), 
                '.P0', str(specdict['P0']), 
                '.box', str(specdict['Lbox'])
                ])

        elif self.type == 'bk': 
            # (hardcoded)
            spectrum_str = ''.join([
                '.grid', str(specdict['Ngrid']), 
                '.nmax40.ncut3.s3', 
                '.P0', str(specdict['P0']), 
                '.box', str(specdict['Lbox'])
                ])

        else: 
            raise NotImplementedError()
    
        file_name = ''.join([
            spec_dir, 
            spec_str,
            gal_file, 
            specparam_str
            ])

        return file_name

    def build(self): 
        """ Calculate power/bispectrum of simulated/observed data catalog 
        """
        
        if 'clobber' not in (self.kwargs).keys(): 
            bool_clobber = False
        else: 
            bool_clobber = self.kwargs['clobber']
        
        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']
        specdict = (self.cat_corr)['spec'] 

        if corrdict['name'] == 'fourier_tophat':
            if self.ell != 2: 
                raise ValueError

            true_cat_corr = {
                    'catalog': catdict, 
                    'correction': {'name': 'true'}
                    }
            tr_gal = Data('data', true_cat_corr)

            fourier_tophat_Pk(self.cat_corr, self.file_name, tr_gal.file_name)
            return None

        spec_type = self.type

        codeclass = Fcode(spec_type, self.cat_corr) 
        spec_code = codeclass.code
        spec_exe = codeclass.fexe()
        
        # code and exe modification time 
        speccode_t_mod, specexe_t_mod = codeclass.mod_time()

        if specexe_t_mod < speccode_t_mod: 
            codeclass.compile()

        # fft files 
        datafft = Fft('data', self.cat_corr, **self.kwargs)
        if not os.path.isfile(datafft.file_name+'_0') or bool_clobber:
            datafft.build()

        randfft = Fft('random', self.cat_corr, **self.kwargs)
        if not os.path.isfile(randfft.file_name+'_0'): 
            randfft.build()
        
        spec_cmd = codeclass.commandline_call(
                datafft = datafft.file_name, 
                randfft = randfft.file_name, 
                powerfile = self.file_name
                )
        print spec_cmd

        if any([not os.path.isfile(self.file_name), bool_clobber]):
            print ''
            print '-----------------------'
            print 'Constructing '
            print self.file_name  
            print '-----------------------'
            print ''
            print spec_cmd
            print '-----------------------'

            subprocess.call(spec_cmd.split())
        else: 
            print ''
            print '-----------------------'
            print self.file_name  
            print 'Already Exists'
            print '-----------------------'
            print ''

        return None

if __name__=='__main__':
    cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'fourier_tophat', 'fs': 1.0, 'rc': 0.43, 'k_fit': 0.7, 'k_fixed': 0.84}
            }
            #'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 3.9, 'fpeak': 0.68} 
            #}
    spectrum = Spec('pk', cat_corr, ell=2, Ngrid=960)
    print spectrum.file()
    print spectrum.build()

"""
def build_fibcol_bispec(**cat_corr): 
    '''
    Given catalog_correction dictionary, construct bispec file 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    spec = cat_corr['spec'] 
    
    # data file name 
    data = fc_data.galaxy_data('data', readdata=False, **cat_corr) 

    # FFt file names 
    fft_file = get_fibcol_fft_file('data', **cat_corr) 
    fft_rand_file = get_fibcol_fft_file('random', **cat_corr) 

    bispec = fc_spec.spec('bispec', **cat_corr) 
    bispec_file = bispec.file_name 

    bispec_code = fc_util.fortran_code('bispec', **cat_corr) 
    bispec_exe = fc_util.fortran_code2exe(power_code)
    
    # code and exe modification time 
    power_code_mod_time = time.ctime(os.path.getmtime(power_code))
    power_exe_mod_time = time.ctime(os.path.getmtime(power_exe))

    # if code was changed since exe file was last compiled then 
    # compile power code 
    if (power_exe_mod_time < power_code_mod_time) or (os.path.isfile(bispec_exe) == False): 
        fc_util.compile_fortran_code(bispec_code) 

    if catalog['name'].lower() == 'lasdamasgeo': 
        bispec_cmd = ' '.join([bispec_exe, '2', fft_rand_file, fft_file, bispec_file]) 
        print power_cmd
        subprocess.call(power_cmd.split()) 
            
    return power_file  
"""
