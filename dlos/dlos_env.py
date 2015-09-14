'''

Code to investigate environment dependence on the 
line-of-sight displacement

Author(s): ChangHoon Hahn


'''
import numpy as np
import scipy as sp 
import os.path
import cosmolopy as cosmos

# --- Local ---
from dlos import Dlos
from spec.data import Data
from galenv.galenv import d_NN

class DlosEnv(Dlos): 

    def __init__(self, cat_corr, n_NN=3, **kwargs):
        """ Child class of Dlos class that describes the nth nearest neighbor distance of 
        upweighted galaxies in fiber collided pair

        Notes
        -----
        * Very clunky because it has to communicate with dLOS parent class 
        """
        self.n_NN = n_NN

        super(DlosEnv, self).__init__(cat_corr, **kwargs)

        self.dlos = None
        self.env = None 

        self.file_name = self.file()
        self.dlos_file = super(DlosEnv, self).file()

    def file(self): 
        """ Name of dLOS + galaxy environment file
        """

        dlos_filename = super(DlosEnv, self).file()

        nNN_str = ''.join([
            'DLOS_d', 
            str(self.n_NN),
            'NN_'
            ]) 
        
        file_name = nNN_str.join( 
                dlos_filename.split('DLOS_')
                ) 

        return file_name 

    def build(self): 
        """ Calculate the nth nearest neighbor distances for fiber collided pairs 
        """
        self.kwargs.pop('clobber', None)
        if self.dlos == None: 
            super(DlosEnv, self).read()

        # dNN for upweighted galaxy in fiber collided 
        # pairs. Each fiber collided pair corresponds to 
        # a value of dLOS 
        dataclass = Data('data', self.cat_corr) 
        dataclass.read() 
        NN_dist = d_NN_dataclass( 
                self.upw_ra, 
                self.upw_dec, 
                self.upw_z, 
                dataclass, 
                n_NN = self.n_NN 
                ) 
    
        # each value of d_NN corresponds to a dLOS value 
        # in dLOS file 
        print self.file_name
        np.savetxt(self.file_name, 
                np.c_[NN_dist], 
                fmt=['%10.5f'],
                header='Columns : d_NN'
                ) 

        return None 

    def read(self, **kwargs): 
        """ Read both dLOS and dNN values 
        """

        if not os.path.isfile(self.file_name):
            self.build()
        elif 'clobber' in self.kwargs.keys():
            if self.kwargs['clobber']: 
                self.build()
        
        if not os.path.isfile(self.dlos_file):
            super(DlosEnv, self).build()
        
        # read dLOS file from parent class
        super(DlosEnv, self).read()

        self.env = np.loadtxt(
                self.file_name, 
                skiprows=1, 
                unpack=True, 
                usecols=[0]
                )

        return None 

def fpeak_dNN(dNN, cat_corr, n_NN=3, **kwargs): 
    ''' fpeak value given nth nearest neighbor distance

    Notes
    -----
    * Upweight correction automatically assigned
    '''
    cat_corr['correction'] = {'name': 'upweight'}

    comb_dlos = DlosEnv(cat_corr, n_NN=n_NN)
    # fpeak file 
    dlos_dir =  '/'.join( (comb_dlos.file_name).split('/')[:-1] )+'/'
    dlos_file = (comb_dlos.file_name).split('/')[-1]
    fpeak_file = ''.join([dlos_dir, 'fpeak_env_', dlos_file]) 
    
    if 'truefpeak' in kwargs.keys():
        if kwargs['truefpeak']: 
            dNN_avg, fpeaks, sigmas = np.loadtxt(fpeak_file, skiprows=2, unpack=True, usecols=[0,1,2])
    else: 
        dNN_avg, fpeaks, sigmas = np.loadtxt(fpeak_file, skiprows=2, unpack=True, usecols=[0,3,2])
    
    if isinstance(dNN, float):
        dNN_list = np.array([dNN])
    else: 
        dNN_list = dNN
    
    #min_indx = [] 
    #for dNN_i in dNN_list:
    #    min_indx.append((np.abs(dNN_avg - dNN_i)).argmin())

    #interp_fpeak = sp.interpolate.interp1d(dNN_avg, fpeaks, kind='linear')
        
    return np.interp( dNN, dNN_avg, fpeaks )

if __name__=="__main__": 
    pass
    #cat_corr = {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}} 
    #for nnn in [4,10]: 
    #    comb_dlos = DlosEnv(cat_corr, n_NN=nnn)
    #    comb_dlos.fpeak_env(stepsize=2)
    #print fpeak_dNN([1.2, 3.4, 5.7], cat_corr, n_NN=3) 
