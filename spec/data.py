'''

Code to handle galaxy data for FiberCollisions 

Author(s): ChangHoon Hahn 


'''
import scipy as sp 
import time 
import random
import os.path
import subprocess
import cosmolopy as cosmos

# --- Local ----
from util import catalog as cata 
from corrections import Corr
from corrections import Rand

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
        if DorR not in ['data', 'random']: 
            raise ValueError("DorR")

        self.cat_corr = cat_corr    # catalog and correction metadata 
        self.kwargs = kwargs    # pass extra input along
        if DorR == 'data':         
            self.corrclass = Corr(self.cat_corr, **self.kwargs)
        else: 
            self.corrclass = Rand(self.cat_corr, **self.kwargs)


        self.file_name = self.File()  # file name 
        self.type = DorR    # type (data or random)

        # galaxy properties
        self.ra = None
        self.dec = None
        self.z = None
        self.weight = None 

        self.cosmo = None   # cosmology of catalog 
    
    def Build(self): 
        """ Calculate galaxy/random catalog
        """
        self.file_name = (self.corrclass).file()

        return None 

    def File(self): 
        """ Name of ASCII file of Data/Random catalogy
        """

        self.file_name = (self.corrclass).file()
        
        return self.file_name

    def Cosmo(self): 
    
        try: 
            if self.kwargs['cosmology'] == 'survey': 
                # survey cosmology
                cat = cata.catalog(self.cat_corr)
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
        self.cosmo = cosmo 

        return self.cosmo 

if __name__ == '__main__':
    cat_corr = {'catalog': {'name': 'cmass', 'n_mock': 1}, 'correction': {'name': 'upweight'}}
    corrdata = Data('random', cat_corr)
    print corrdata.file_name
