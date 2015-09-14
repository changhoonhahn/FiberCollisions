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
from util.util import radecz_to_xyz

def d_NN(radecz1, radecz2, n_NN=3, cosmos='default'): 
    ''' Calculate nth nearest neighbor distance of a set of ra, dec, z 
    w.r.t second another set of ra, dec, z provided. 
    
    (n = 3 by default) 

    --------------------------------------------------------------------------
    Parameters
    --------------------------------------------------------------------------
    radecz1 : [ra1, dec1, z1]
    radecz2 : [ra2, dec2, z2]
    n_NN : n th nearest neighbor (default 3) 
    cosmos : Unless otherwise specified, default cosmology is 'default'. 
    This automatically sets up default cosmology with Omega_M = 0.31 
    --------------------------------------------------------------------------
    Notes 
    --------------------------------------------------------------------------
    * Currently assumed that radecz1 is contained within radecz2. In order to 
    deal with overlap, the n+1th distance is actually returned. This should 
    ultimately be corrected, but for current purposes it is overlooked. Will be
    revisited. 

    '''

    if cosmos == 'default': 
        
        omega_m = 0.31 
        data_cosmos = {} 
        data_cosmos['omega_M_0'] = omega_m 
        data_cosmos['omega_lambda_0'] = 1.0 - omega_m 
        data_cosmos['h'] = 0.676
        data_cosmos = cosmos.distance.set_omega_k_0(cosmo) 

    else: 
        data_cosmos = cosmos
   
    targ_ra, targ_dec, targ_redshift = radecz1
    cat_ra, cat_dec, cat_redshift = radecz2

    # convert to cartesian coordinates to set up 
    # KDE tree
    targ_x, targ_y, targ_z = radecz_to_xyz(
            targ_ra, 
            targ_dec, 
            targ_redshift, 
            **data_cosmos
            )

    cat_x, cat_y, cat_z = radecz_to_xyz(
            cat_ra, 
            cat_dec, 
            cat_redshift, 
            **data_cosmos
            )

    # set up KD Tree
    tree = sp.spatial.KDTree( zip(cat_x, cat_y, cat_z) )       

    # query KD Tree for n+1 neighbors because it counts itself
    distance, index = tree.query( 
            zip( targ_x, targ_y, targ_z ), 
            k = n_NN + 1 
            ) 

    return distance[:, n_NN]

def d_NN_dataclass(ra, dec, redshift, dataclass, n_NN = 3):
    """ Calculate nth nearest non-fiber collided neighbor distances in catalog 
    for input ra, dec, z. When counting neighbors only non fiber collided galaxies are
    considered. 
    
    Parameters
    ----------
    ra : Right Ascension  
    dec : Declination
    redshift : Redshift (z)
    dataclass : Data class object. Correction dictionary should be {'name': 'upweight'}  
    n_NN : n th nearest neighbor (default 3) 

    """
    
    # cosmology of data class
    data_cosmos = dataclass.cosmo()
    
    if 'weight' in dataclass.datacolumns:    
        # accounts for poorly planned column names 
        fc_weights = dataclass.weight
    else: 
        fc_weights = dataclass.wfc
    
    hasfibers = np.where(fc_weights > 0) 
    
    radecz1 = [ra, dec, redshift] 
    radecz2 = [(dataclass.ra)[hasfibers], (dataclass.dec)[hasfibers], (dataclass.z)[hasfibers]]

    distance = d_NN(
            radecz1, 
            radecz2, 
            n_NN = n_NN, 
            cosmos = data_cosmos
            )

    return distance 

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
