"""

"""
from Spectrum.fft import Fft
from corr_corrdata import CorrCorrData

class CorrFft(Fft): 
    
    def __init__(self, DorR, cat_corr, **kwargs): 
        '''
        '''
        super(CorrFft, self).__init__(DorR, cat_corr, **kwargs)

        if self.type == 'data': 
            self.galdata = CorrCorrData(self.cat_corr, **self.kwargs)  # data class 
        self.file_name = self.file()

    def file(self): 
        '''
        Return power/bispectrum file name strings. Spectrum.spec plus exceptions.
        '''
        if self.type == 'data': 
            if self.cat_corr['correction']['name'].lower() == 'floriansn': 
                fc_file = super(CorrFft, self).file()
                return '.floriansn.dat'.join(fc_file.split('.dat'))
            elif self.cat_corr['correction']['name'].lower() == 'hectorsn': 
                fc_file = super(CorrFft, self).file()
                return '.hectorsn.dat'.join(fc_file.split('.dat'))
            else: 
                return super(CorrFft, self).file()
        else: 
            return super(CorrFft, self).file()
