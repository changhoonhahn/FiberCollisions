
from Spectrum.corrdata import CorrData

# --- Corrections ---
from corrections.true import TrueCorr
from corrections.dlospeak import DlospeakCorr
from corrections.fibcollided import UpweightCorr
from corrections.dlospeak_env import DlospeakEnvCorr
from corrections.dlospeak_flex import DlospeakFlexCorr
from corrections.dlospeak_known import DlospeakKnownCorr
from corrections.dlospeak_photoz import DlospeakPhotozCorr
from corrections.dlospeak_shuffle import DlospeakShuffleCorr
from corrections.dlospeak_tailonly import DlospeakTailonlyCorr
from corrections.dlospeak_peakonly import DlospeakPeakonlyCorr
from corrections.photoz_corr import PhotozCorr
from corrections.fourier_tophat import FourierTophatCorr

class CorrCorrData(CorrData): 

    def __init__(self, cat_corr, **kwargs): 

        self.corrclass_dict = { 
                'true': TrueCorr,
                'upweight': UpweightCorr, 
                'photoz': PhotozCorr,
                'dlospeak': DlospeakCorr, 
                'dlospeakenv': DlospeakEnvCorr, 
                'dlospeakphotoz': DlospeakPhotozCorr,
                'dlospeakknown': DlospeakKnownCorr,
                'dlospeak.flex': DlospeakFlexCorr,
                'dlospeak.shuffle': DlospeakShuffleCorr,
                'dlospeak.tailonly': DlospeakTailonlyCorr, 
                'dlospeak.peakonly': DlospeakPeakonlyCorr,
                'fourier_tophat': FourierTophatCorr
                }

        super(CorrCorrData, self).__init__(cat_corr, **kwargs)
