'''

FiberCollision Correction using peak of line-of-sight 
displacement distribution. Code could use overall editing


'''
import numpy as np
import scipy as sp 
import cosmolopy as cosmos
from scipy.stats import norm
from scipy.stats import expon
from scipy.stats import binom
import time 

# --- Local ---
from util.direc import direc
from util.catalog import Catalog
from corrections import Corrections

class FourierTophatCorr(Corrections): 

    def __init__(self, cat_corr, **kwargs): 
        '''
        Child class of Correction class in corrections.py. This child class solely 
        exists to have the FourierTophat correction within the same class framework
        as the other spec.spec.Spec class objects.

        Notes
        -----
        '''

        super(FourierTophatCorr, self).__init__(cat_corr, **kwargs)
        
        self.corrstr() 
    
    def corrstr(self): 
        """ Specify correction string for dlos peak correction 
        """
        corr = (self.cat_corr)['correction']

        if not all(x in corr.keys() for x in ['fs', 'rc', 'k_fit', 'k_fixed']): 
            raise KeyError("Specify fs (collided fraction), rc (fibcollision comoving radius), \n k_fit (fitting range of k), k_fixed (for the power law) in the correction dictionary")

        self.corr_str = ''.join([
            '.', corr['name'].lower(), 
            '.fs', str(round(corr['fs'], 1)), 
            '.rc', str(round(corr['rc'], 2)), 
            '.kfit', str(round(corr['k_fit'], 2)), 
            '.kfixed', str(round(corr['k_fixed'], 2))
            ])
        return self.corr_str
