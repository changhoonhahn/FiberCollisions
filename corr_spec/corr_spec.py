
import numpy as np
from Spectrum.spec import Spec

from corr_fft import CorrFft
from corr_corrdata import CorrCorrData
    
from fourier_corr.pk_corr import fourier_tophat_Pk
#from pk_corr import fourier_tophat_Pk

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
            if self.cat_corr['correction']['name'].lower() == 'true': 
                corr_str = ''
            elif self.cat_corr['correction']['name'].lower() == 'fourier_tophat': 
                corr = self.cat_corr['correction'].copy()
                corr_str = ''.join([
                    '.', corr['name'].lower(), 
                    '.fs', str(round(corr['fs'], 1)), 
                    '.rc', str(round(corr['rc'], 2)), 
                    '.kfit', str(round(corr['k_fit'], 2)), 
                    '.kfixed', str(round(corr['k_fixed'], 2))])
            else: 
                raise NotImplementedError

            file_name = ''.join([
                spec_dir, 
                'power3600z_BoxN', 
                str(self.cat_corr['catalog']['n_mock']), 
                corr_str, 
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
    
    def build(self):
        '''
        '''
        if self.cat_corr['correction']['name'].lower() == 'fourier_tophat': 
            fourier_tophat_cat_corr = self.cat_corr.copy()

            self.cat_corr = {
                    'catalog': self.cat_corr['catalog'].copy(), 
                    'correction': {'name': 'true'},
                    'spec': self.cat_corr['spec'].copy()
                    }
            if self.cat_corr['catalog']['name'] != 'nseriesbox': 
                self.gal_data = CorrCorrData(self.cat_corr)
            tr_file_name = self.file()
            self.cat_corr = fourier_tophat_cat_corr.copy()
            fourier_tophat_Pk(self.cat_corr, self.file_name, tr_file_name)

            return None
        else: 
            return super(CorrSpec, self).build()

if __name__=='__main__': 
    for i_mock in range(1,8): 
        cat_corr = {
                'catalog': {'name': 'nseriesbox', 'n_mock': i_mock}, 
                'correction': {'name': 'fourier_tophat', 'fs': 1.0, 'rc': 0.43, 'k_fit': 4.3, 'k_fixed': 4.34}
                }
        spectrum = CorrSpec('pk', cat_corr, ell=0, Ngrid=960)
        print spectrum.file()
        spectrum.build()
        #'correction': {'name': 'fourier_tophat', 'fs': 1.0, 'rc': 0.43, 'k_fit': 4.0, 'k_fixed': 4.34}, 
