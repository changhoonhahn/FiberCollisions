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
