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
import warnings 
import matplotlib.pyplot as plt

# --- Local ----
import correction.correction 
import correction.true

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
        self.cat_corr = cat_corr    # catalog and correction metadata 
        self.file_name = self.File(DorR, cat_corr)  # file name 
        self.Type = DorR    # data or random 

        # galaxy properties
        self.ra = None
        self.dec = None
        self.z = None
        self.weight = None 

        self.cosmo = None   # cosmology of catalog 
    
    def Build(self): 
        """ Calculate galaxy/random catalog
        """
        correction(self.cat_corr)

    def File(self): 
        """ Name of ASCII file of Data/Random catalogy
        """

        DorR = self.Type        # data or random 

    def Cosmo(self): 
        # survey cosmology 
        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.676
        cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        self.cosmo = cosmo 

def comdis2z(comdis, **cosmo): 
    ''' Given comoving distance and cosmology, determine z 
    Comoving distance *has* to be in Mpc/h
    '''
    z_arr = np.array([0.0+0.05*np.float(i) for i in range(21)]) 
    dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo)*cosmo['h']

    dmz_spline = sp.interpolate.interp1d(dm_arr, z_arr, kind='cubic') 
    
    z = dmz_spline(comdis)
    return z 
