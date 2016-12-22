'''

Tophat convolved true Powerspectrum monopole/quadrupole. 
Convolving the tophat function with the true P_l(k) reproduces
the effects of fiber collisions on the power spectrum. So the 
tophat convolved power spectrum should be compared to the 
upweighted power spectrum. 


'''
import numpy as np

# --- Local --- 
from true import TrueCorr
from fourier_corr import p_k_mu as pkmu
from fourier_corr import fourier_corr as fourcorr

from corr_spec.corr_spec import CorrSpec


class TophatConvCorr(Corrections): 

    def __init__(self, cat_corr, **kwargs): 
        ''' Child class of Correction class in corrections.py. This child class solely 
        exists to have the Tophat Convolution correction within the same class framework
        as the other spec.spec.Spec class objects.

        Notes
        -----
        '''

        super(TophatConvCorr, self).__init__(cat_corr, **kwargs)
        
        self.corrstr() 
    
    def corrstr(self): 
        ''' Specify correction string that identifies the correction. The keys included
        in the corrstr include 
    
        fs : float
            Specifies the fraction of the survey that suffers form fiber collisions. For 
            Nseries and CMASS fs = 0.6. 
        rc : float
            Specifies the comoving physical scale of the fiber collision angular scale. 
            At z = 0.55, rc = 0.43. (This is one of the key assumptions made in the 
            correction method). 
        fold : int  
            The P(k, mu) value is supplemented with a folded P(k, mu) for high k values. 
            "fold" specifies the number of folds. fold = 5 or 10. 10 is stable; 5 is not. 
        rebin : int 
            Due to the noise in P(k, mu) at high k, the k values are rebin at k > 0.1. 
            rebin specifies the number of bins. rebin = 20. 
        '''
        corr = (self.cat_corr)['correction']
        # correction keys 
        if not all(x in corr.keys() for x in ['fs', 'rc', 'fold', 'rebin']): 
            raise KeyError(
                    "Specify fs (collided fraction), rc (fibcollision comoving radius), \n k_fit (fitting range of k), k_fixed (for the power law) in the correction dictionary"
                    )

        self.corr_str = ''.join([
            '.', corr['name'].lower(), 
            '.fs', str(round(corr['fs'], 1)), 
            '.rc', str(round(corr['rc'], 2)), 
            '.fold', str(corr['fold']), 
            '.rebin', str(corr['rebin'])
            ])
        return self.corr_str

    def build(self): 
        ''' Build the tophat convolved power spectrum. The correction   
        '''
        catdict = (self.cat_corr)['catalog'] # catalog dictionary
        corrdict = (self.cat_corr)['correction'] # correction dictionary
        # Only Nseries is supporte
        if catdict['name'].lower() != 'nseries': 
            raise ValueError('Only Nseries catalog is supported for this correction') 
        # correction parameters
        fs = corrdict['fs']
        rc = corrdict['rc']
        fold = corrdict['fold']
        rebin = corrdict['rebin']
        # true power spectrum  
        true_cat_corr = {
                'catalog': cat_dict,
                'correction': {'name': 'true'}
                }
        true_spec = CorrSpec('pk', true_cat_corr, ell= 


    def _read_true_Pk(catdict):
        ''' Convenience function to read in true Powerspectrum given the
        catalog dictionary. Currently only implemented for Nseries. This is only set up 
        this way because
        '''

