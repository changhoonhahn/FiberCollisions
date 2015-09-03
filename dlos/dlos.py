"""

Line-of-Sight Displacement 

Author(s): ChangHoon Hahn (NYU CCPP)

"""
import numpy as np
import scipy as sp 

from util.direct import direc

class Dlos: 

    def __init__(self, cat_corr, **kwargs):
        """ Class describing line-of-sight displacement 
        """

        if (cat_corr['correction'])['name'] != 'upweight': 
            cat_corr['correction'] = {'name': 'upweight'}
        self.cat_corr = cat_corr    
        self.kwargs = kwargs

        self.dlos = None
        self.targ_z = None
        self.neigh_z = None 

        self.file_name = self.file()   

    def file(self): 
        """ Name of line-of-sight displacement file 
        """
        catdict = self.cat_corr['catalog']
        corrdict = self.cat_corr['correction']

        dlos_dir = direc('data', self.cat_corr)
        dlos_str = 'DLOS_'

        dataclass = Data(self.type, self.cat_corr, **self.kwargs)  # data class 
        data_file = dataclass.file_name # galaxy data file

        # add bells and whistles here later
        # add bells and whistles here later
        # add bells and whistles here later

        file_name = ''.join([
            dlos_dir, 
            dlos_str, 
            data_file
            ])

        return file_name 

    def read(self): 
        """ Read dLOS data 
        """
        
        data_cols = self.datacolumns
        data_list = np.loadtxt(
                self.file_name, 
                unpack=True, 
                usecols=range(len(data_cols))
                )
        
        for i_col, data_col in data_cols:       # store data as attributes
            setattr(self, data_col, data_list[i_col])   

        return None

    def build(self): 
        """ Calculate dLOS of fiber collided pairs in the galaxy catalog.
        Uses IDL code with spherematch.pro. 
        """



    def datacolumns(self):
        """ Data columns of dLOS data
        """
        data_cols = ['dlos', 'targ_ra', 'targ_dec', 'targ_z', 'neigh_ra', 'neigh_dec', 'neigh_z']

        return data_cols 

    def datacols_fmt(self): 
        """ Data column data format 
        """
        data_fmt = ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f']

        return data_fmt
