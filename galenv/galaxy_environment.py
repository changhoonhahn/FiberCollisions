'''

Galaxy environment to be used in the FCs project


Authro(s): ChangHoon Hahn 

'''
import numpy as np 
import scipy as sp
import os.path
import time 
import cosmolopy as cosmos

# -- Local -- 
from spec.data import Data
from util.util import radecz_to_xyz

def d_NN(ra, dec, redshift, cat_corr, n=3): 
    ''' Calculate nth nearest non-fiber collided neighbor distances in catalog 
    for input ra, dec, z. When counting neighbors only non fiber collided galaxies are
    considered. 
    
    (n = 3 by default) 

    Parameters
    ----------
    ra : Target RA array 
    dec : Target Dec array 
    redshift : Target Redshift array 
    n : n th nearest neighbor (default 3) 
    cat_corr : catalog correction dictionary 

    Notes 
    -----

    '''

    corrdict = cat_corr['correction']
    
    if corrdict['name'] != 'upweight': 
        corrdict = {'name': 'upweight'} 


    dataclass = Data('data', cat_corr)
    dataclass.read()
    data_cosmo = dataclass.cosmo()
    
    if 'weight' in dataclass.datacolumns:    
        # accounts for poorly planned column names 
        fc_weights = dataclass.weight
    else: 
        fc_weights = dataclass.wfc
    
    hasfibers = np.where(fc_weights > 0) 
    
    cat_x, cat_y, cat_z = radecz_to_xyz(
            (dataclass.ra)[hasfiber], 
            (dataclass.dec)[hasfiber], 
            (dataclass.z)[hasfiber], 
            **data_cosmo
            )

    targ_x, targ_y, targ_z = radecz_to_xyz(
            ra, 
            dec, 
            redshift, 
            **data_cosmo
            )

    # set up KD Tree
    tree = scipy.spatial.KDTree( zip(cat_x, cat_y, cat_z) )       

    # query KD Tree for n+1 neighbors because it counts itself
    distance, index = tree.query( zip( targ_x, targ_y, targ_z ), k = n+1 ) 

    return distance[:, n]
"""
def n_nearest(cat_corr, n=3): 
    ''' Calculate nth nearest neighbor distances for upweighted galaxies
    in specified galaxy catalog (n = 3 by default) 
    '''

    corrdict = cat_corr['correction']
            
    if corrdict['name'] != 'upweight': 
        corrdict = {'name': 'upweight'} 

    gal_data = fc_data.galaxy_data('data', readdata=True, **cat_corr)       # read data 
    
    if 'weight' in gal_data.__dict__.keys():        # accounts for poorly planned column names 
        fc_weights = gal_data.weight
    else: 
        fc_weights = gal_data.wfc
    
    in_catalog = fc_weights > 0 
    upweighted = fc_weights > 1         # indices of upweighted galaxies
    
    x, y, z = fc_util.radecz_to_xyz(gal_data.ra, gal_data.dec, gal_data.z, **gal_data.cosmo) 

    tree = sp.spatial.KDTree( zip( x[in_catalog], y[in_catalog], z[in_catalog] ) )   # set up KD Tree
    
    # query KD Tree for n+1 neighbors because it counts itself
    distance, index = tree.query( zip( x[upweighted], y[upweighted], z[upweighted] ), k = n+1 ) 
    print np.min(distance[:, n]), np.max(distance[:, n])
    return distance[:, n]
"""

def dlos_d_NN(n=3, clobber=False, **cat_corr): 
    ''' Get nth nearest neighbor distance for dLOS file 
    '''

    NN_dist_file = dlos_d_NN_file(n=n, **cat_corr)      # dNN file  
    
    if not os.path.isfile(NN_dist_file) or clobber: 
        build_dlos_d_NN(n=n, **cat_corr)  
    else: 
        pass

    dNN = np.loadtxt(NN_dist_file, unpack=True, usecols=[0])
    return dNN 
