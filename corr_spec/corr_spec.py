from Spectrum.spec import Spec

from corr_fft import CorrFft
from corr_corrdata import CorrCorrData

class CorrSpec(Spec): 
    
    def __init__(self, spectype, cat_corr, ell=None, **kwargs): 
        '''
        '''
        super(CorrSpec, self).__init__(spectype, cat_corr, ell=ell, **kwargs)

        self.gal_data = CorrCorrData(self.cat_corr, **self.kwargs)
        self.datafft = CorrFft('data', self.cat_corr, **self.kwargs)
        self.randfft = CorrFft('random', self.cat_corr, **self.kwargs)

        self.file_name = self.file()

if __name__=='__main__': 
    cat_corr = {'catalog': {'name': 'nseries', 'n_mock': 1}, 'correction': {'name': 'upweight'}}
    spectrum = CorrSpec('bk', cat_corr, Ngrid=360)
    print spectrum.file()
    print spectrum.build()
