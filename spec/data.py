'''

Code to handle galaxy data for FiberCollisions 

Author(s): ChangHoon Hahn 


'''


import numpy as np
import scipy as sp 
import time 
import random
import os.path
import subprocess
import cosmolopy as cosmos

# --- Local ----
from util.catalog import Catalog 
from corrections.rand import Rand
from corrections.true import TrueCorr
from corrections.dlospeak import DlospeakCorr
from corrections.fibcollided import UpweightCorr
from corrections.dlospeak_env import DlospeakEnvCorr
from corrections.dlospeak_known import DlospeakKnownCorr
from corrections.dlospeak_photoz import DlospeakPhotozCorr
from corrections.photoz_corr import PhotozCorr

# Classes ------------------------------------------------------------
class Data(object): 

    def __init__(self, DorR, cat_corr, **kwargs): 
        """ A class describing galaxy/random catalog of simulations or BOSS data. 
        Corresponds to an ASCII file with galaxy/random catalog 

        Parameters
        ----------

        DorR : 'data' or 'random'
        cat_corr :  Catalog correction Dictionary 

        """ 

        self.cat_corr = cat_corr.copy()    
        self.kwargs = kwargs    

        if DorR == 'data':         

            # correction class dictionary 
            corrclass_dict = { 
                    'true': TrueCorr,
                    'upweight': UpweightCorr, 
                    'photoz': PhotozCorr,
                    'dlospeak': DlospeakCorr, 
                    'dlospeakenv': DlospeakEnvCorr, 
                    'dlospeakphotoz': DlospeakPhotozCorr,
                    'dlospeakknown': DlospeakKnownCorr
                    }

            corr_name = ((self.cat_corr)['correction'])['name']
            if corr_name not in corrclass_dict.keys():
                raise NameError()

            self.corrclass = corrclass_dict[corr_name](cat_corr, **kwargs)

        elif DorR == 'random': 
            self.corrclass = Rand(cat_corr, **kwargs)
        else: 
            raise ValueError("DorR")

        self.file_name = self.file()  # file name 
        self.type = DorR    # type (data or random)

        # galaxy properties
        self.ra = None
        self.dec = None
        self.z = None
        self.wfc = None
        self.weight = None 
        self.comp = None

        self.cosmos = None   # cosmology of catalog 

    def read(self): 
        """ Read galaxy/random catalog data 
        """

        data_cols = self.corrclass.datacolumns()
        self.datacolumns = self.corrclass.datacolumns()
    
        datah = np.loadtxt(
                self.file_name, 
                skiprows=1, 
                unpack=True, 
                usecols=range(len(data_cols))
                )

        for i_col, col in enumerate(data_cols): 
            setattr(self, col, datah[i_col])
        
        return None
    
    def build(self): 
        """ Calculate galaxy/random catalog
        """
        (self.corrclass).build()

        return None 

    def file(self): 
        """ Name of ASCII file of Data/Random catalogy
        """

        self.file_name = (self.corrclass).file()
        
        return self.file_name

    def cosmo(self): 
    
        try: 
            if self.kwargs['cosmology'] == 'survey': 
                # survey cosmology
                cat = Catalog(self.cat_corr)
                self.cosmo = cat.cosmo()

                return self.cosmo
            else: 
                # default fiducial cosmology (hardcoded)
                omega_m = 0.31 

        except KeyError: 
            omega_m = 0.31  # default 

        # survey cosmology 
        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.676
        cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        self.cosmos = cosmo 

        return self.cosmos

if __name__ == '__main__':

    for i_mock in xrange(1,85): 
        for corr in ['true', 'upweight', 'photoz']:
            cat_corr = {
                    'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
                    'correction': {'name': 'photoz'}
                    }
            corrclass = Data('data', cat_corr, clobber=True)
            corrclass.build()
