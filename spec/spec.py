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

# Classes ------------------------------------------------------------
class Spec(object): 

    def __init__(self, spectype, cat_corr, **kwargs):
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
            # default spectrum parameters
            cat_corr['spec'] = {
                    'P0': 20000, #P0 
                    'Lbox': 3600, 
                    'Ngrid':360, 
                    'quad': False
                    }
        
        self.cat_corr = cat_corr.copy()
        self.kwargs = kwargs

        self.file_name = self.file()
    
    def read(self): 
        """ Read power/bispectrum of simulated/observed data catalog
        """

        spec_dict = self.cat_corr['spec']
    
        if self.type == 'pk': 
                
            if not spec_dict['quad']: 
                col_index = [0, 1]
                data_cols = ['k', 'p0k']
            else: 
                raise NotImplementedError()

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

        spec_dir = ''.join([
            '/mount/riachuelo1/hahn/', 
            self.type, 
            '/'
            ]) 
        
        # powerspectrum or bispectrum 
        if self.type == 'pk': 
            spec_str = 'POWER_'
        elif self.type == 'bk':
            spec_str = 'BISP_'

        if 'quad' not in specdict.keys(): 
            specdict['quad'] = False
        
        if specdict['quad']:          
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
        
        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']
        specdict = (self.cat_corr)['spec'] 
       
        if 'quad' not in specdict.keys():
            specdict['quad'] = False

        if not specdict['quad']:       # quadrupole/monopole
            spec_type = self.type
        else:  
            # not sure what to do here yet 
            NotImplementedError()

        codeclass = Fcode(spec_type, self.cat_corr) 
        spec_code = codeclass.code
        spec_exe = codeclass.fexe()
        
        # code and exe modification time 
        speccode_t_mod, specexe_t_mod = codeclass.mod_time()

        if specexe_t_mod < speccode_t_mod: 
            codeclass.compile()

        # fft files 
        datafft = Fft('data', self.cat_corr, **self.kwargs)
        if not os.path.isfile(datafft.file_name): 
            datafft.build()

        randfft = Fft('random', self.cat_corr, **self.kwargs)
        if not os.path.isfile(randfft.file_name): 
            randfft.build()
        
        spec_cmd = codeclass.commandline_call(
                datafft = datafft.file_name, 
                randfft = randfft.file_name, 
                powerfile = self.file_name
                )

        if 'clobber' not in (self.kwargs).keys(): 
            bool_clobber = False

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
            'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 3.9, 'fpeak': 0.68} 
            }
    spectrum = Spec('pk', cat_corr)
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
