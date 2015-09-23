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

class DlosPhotoz(Dlos): 

    def __init__(self, cat_corr, **kwargs):
        """ Child class of Dlos class that describes line-of-sight displacement 
        using the photometric redshift of the collided galaxy. 

        dLOS_photoz = Dc(z_upw) - Dc(z_photoz)

        Notes
        -----
        * Very clunky because it has to communicate with dLOS parent class 
        """

        super(DlosPhotoz, self).__init__(cat_corr, **kwargs)

        if self.cat_corr['catalog']['name'] != 'nseries': 
            raise NotImplementedError()

        self.dlos = None
        self.dlos_photoz = None 

        self.file_name = self.file()
        self.dlos_file = super(DlosPhotoz, self).file()

    def file(self): 
        """ Name of dLOS + galaxy environment file
        """

        dlos_filename = super(DlosPhotoz, self).file()

        photoz_str = 'DLOS_photoz_'
        
        file_name = photoz_str.join( 
                dlos_filename.split('DLOS_')
                ) 

        return file_name 

    def build(self): 
        """ Calculate the line-of-sight displacement using assigned 
        photometric redshift
        """
        self.kwargs.pop('clobber', None)

        # Read in mock catalog with assigned photometric redshifts
        # and calculate the line-of-sight displacement between the 
        # upweighted galaxy and the photometric redshift of the 
        # collided galaxy 
        photoz_cat_corr = {
                'catalog': self.cat_corr['catalog'].copy(), 
                'correction': {'name': 'photoz'}
                }
        dataclass = Data('data', photoz_cat_corr) 
        dataclass.read() 

        cosmo = dataclass.cosmo()

        coll = np.where(dataclass.wfc == 0) 
        
        dlos_actual = (cosmos.distance.comoving_distance(dataclass.z[coll], **cosmo) - \
                cosmos.distance.comoving_distance(dataclass.zupw[coll], **cosmo)) * cosmo['h']
        dlos_photoz = (cosmos.distance.comoving_distance(dataclass.photoz[coll], **cosmo) - \
                cosmos.distance.comoving_distance(dataclass.zupw[coll], **cosmo)) * cosmo['h']

        # each value of d_NN corresponds to a dLOS value 
        # in dLOS file 
        print self.file_name
        np.savetxt(self.file_name, 
                np.c_[dlos_actual, dlos_photoz], 
                fmt=['%10.5f', '%10.5f'],
                header='Columns : dLOS, dLOS_photoz'
                ) 

        return None 

    def read(self, **kwargs): 
        """ Read both dLOS and dLOS_photoz values 
        """

        if not os.path.isfile(self.file_name):
            self.build()
        elif 'clobber' in self.kwargs.keys():
            if self.kwargs['clobber']: 
                self.build()

        # read dLOS file from parent class
        super(DlosPhotoz, self).read()

        self.dlos, self.dlos_photoz = np.loadtxt(
                self.file_name, 
                skiprows=1, 
                unpack=True, 
                usecols=[0, 1]
                )
        return None 

if __name__=="__main__": 
    cat_corr = {'catalog': {'name': 'nseries', 'n_mock': 1}, 'correction': {'name': 'photoz'}} 
    dlos_class = DlosPhotoz(cat_corr) 
    dlos_class.build()
