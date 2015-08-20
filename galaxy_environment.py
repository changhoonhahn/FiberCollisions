'''

Code for getting galaxy environment to be used in the FCs project


Authro(s): ChangHoon Hahn 

'''


import numpy as np 
import os.path
import time 
import subprocess
import cosmolopy as cosmos
import scipy

# -- Local -- 
import fibcol_nbar as fc_nbar
import fibcol_data as fc_data
import fibcol_spec as fc_spec
import fibcol_dlos as fc_dlos
import fibcol_utility as fc_util
import fibcol_fft as fc_fft
import plot_fibcol as fc_plot

def n_nearest(n=3, **cat_corr): 
    ''' Get nth nearest neighbor distances for upweighted galaxies (n = 3 by default) 
    '''

    # make sure that correction is upweight 
    if (cat_corr['correction'])['name'].lower() != 'upweight': 
        cat_corr['correction'] = {'name': 'upweight'} 

    gal_data = fc_data.galaxy_data('data', readdata=True, **cat_corr)       # read data 
    
    if 'weight' in gal_data.__dict__.keys():        # accounts for poorly planned column names 
        fc_weights = gal_data.weight
    else: 
        fc_weights = gal_data.wfc
    
    in_catalog = fc_weights > 0 
    upweighted = fc_weights > 1         # indices of upweighted galaxies
    
    x, y, z = fc_util.radecz_to_xyz(gal_data.ra, gal_data.dec, gal_data.z, **gal_data.cosmo) 

    tree = scipy.spatial.KDTree( zip( x[in_catalog], y[in_catalog], z[in_catalog] ) )   # set up KD Tree
    
    # query KD Tree for n+1 neighbors because it counts itself
    distance, index = tree.query( zip( x[upweighted], y[upweighted], z[upweighted] ), k = n+1 ) 
    print np.min(distance[:, n]), np.max(distance[:, n])
    return distance[:, n]

def d_NN(ra, dec, redshift, n=3, **cat_corr): 
    ''' Get nth nearest neighbor distances for input ra, dec, z (n = 3 by default) 

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

    # make sure that correction is upweight 
    cat_corr['correction'] = {'name': 'upweight'} 

    gal_data = fc_data.galaxy_data('data', **cat_corr)      # read data 
    
    if 'weight' in gal_data.__dict__.keys():    
        # accounts for poorly planned column names 
        fc_weights = gal_data.weight
    else: 
        fc_weights = gal_data.wfc
    
    in_catalog = np.where(fc_weights > 0)     # only galaxy with fibers
    
    # convert RA, Dec, z to x,y,z
    x, y, z = fc_util.radecz_to_xyz(
            (gal_data.ra)[in_catalog], (gal_data.dec)[in_catalog], (gal_data.z)[in_catalog], 
            **gal_data.cosmo)    # mock catalog

    targ_x, targ_y, targ_z = fc_util.radecz_to_xyz(ra, dec, redshift, **gal_data.cosmo)     # target

    tree = scipy.spatial.KDTree( zip(x, y, z) )   # set up KD Tree
    
    # query KD Tree for n+1 neighbors because it counts itself
    distance, index = tree.query( zip( targ_x, targ_y, targ_z ), k = n+1 ) 
    return distance[:, n]

def dlos_d_NN_file(n=3, **cat_corr): 
    ''' File name of dNN file for dLOS file 
    '''
    dLOS = fc_dlos.dlos(readdata=False, **cat_corr) 

    NN_dist_str = 'DLOS_d_'+str(n)+'NN_'
    NN_dist_file = NN_dist_str.join((dLOS.file_name).split('DLOS_'))
    return NN_dist_file 

def build_dlos_d_NN(n=3, **cat_corr): 
    ''' build nth nearested neighbor distance for dLOS file 
    '''
    # dlos file
    dLOS = fc_dlos.dlos(**cat_corr) 

    if 'targ_ra' not in dLOS.__dict__.keys():   # make sure target info is stored
        raise NameError('Target information not included in dLOS') 
    
    # obtain dNN values for dLOS target galaxies
    NN_dist = d_NN( dLOS.targ_ra, dLOS.targ_dec, dLOS.targ_z, n=n, **cat_corr) 
    
    # write dNN values to file
    NN_dist_file = dlos_d_NN_file(n=n, **cat_corr) 
    print 'Writing ', NN_dist_file
    np.savetxt(NN_dist_file, 
            np.c_[NN_dist], 
            fmt=['%10.5f']) 

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
