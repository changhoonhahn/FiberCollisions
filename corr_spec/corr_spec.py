
import numpy as np
from Spectrum.spec import Spec

from corr_fft import CorrFft
from corr_corrdata import CorrCorrData

class CorrSpec(Spec): 
    def __init__(self, spectype, cat_corr, ell=None, **kwargs): 
        '''
        '''
        super(CorrSpec, self).__init__(spectype, cat_corr, ell=ell, **kwargs)
    
        try: 
            self.gal_data = CorrCorrData(self.cat_corr, **self.kwargs)
        except NotImplementedError:
            pass
        try: 
            self.datafft = CorrFft('data', self.cat_corr, **self.kwargs)
        except NotImplementedError:
            pass
        try: 
            self.randfft = CorrFft('random', self.cat_corr, **self.kwargs)
        except NotImplementedError:
            pass

        self.file_name = self.file()

    def file(self): 
        '''
        Return power/bispectrum file name strings. Spectrum.spec plus exceptions.
        '''

        if self.cat_corr['catalog']['name'].lower() == 'nseriesbox': 
            spec_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'

            file_name = ''.join([
                spec_dir, 
                'power3600z_BoxN', 
                str(self.cat_corr['catalog']['n_mock']), 
                '.dat'])
            return file_name
        else: 
            return super(CorrSpec, self).file()

    def read(self): 
        '''
        Read bispectrum/powerspectrum file
        '''
        if self.cat_corr['catalog']['name'].lower() == 'nseriesbox': 
            col_index = [0, 6, 2, 3, 4]
            data_cols = ['k', 'p0k', 'p2k', 'p4k', 'p6k']

            spec_data = np.loadtxt(
                        self.file_name, 
                        unpack = True, 
                        usecols = col_index
                        )
            for i_col, col in enumerate(data_cols): 
                if col != 'k': 
                    setattr(self, col, (2.*np.pi)**3 * spec_data[i_col])
                else: 
                    setattr(self, col, spec_data[i_col])
            return None
        else: 
            return super(CorrSpec, self).read()

if __name__=='__main__': 
    cat_corr = {'catalog': {'name': 'nseriesbox', 'n_mock': 1}, 'correction': {'name': 'true'}}
    spectrum = CorrSpec('bk', cat_corr, Ngrid=360)
    print spectrum.file()
    spectrum.read()
    print spectrum.p4k
    #print spectrum.build()
