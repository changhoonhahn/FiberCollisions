'''

Photometric Redshift to supplement fiber collision correction 

Author(s): ChangHoon Hahn 


'''
import numpy as np
import scipy as sp 
import random
import os.path

# --- Local ----
import pyspherematch as pysph
from spec.data import Data
from defutility.fitstables import mrdfits


def match_zspec_zphoto_cmass(): 
    ''' Match photometric redshift to spectroscopic redshift. Spherematch 
    photometric and spectroscopic catalogs to determine redshift errors of 
    photmetric redshifts wrt spectroscopic redshifts. 
    
    Returns [z_spec, z_photo]

    --------------------------------------------------------------------------
    Parameters 
    --------------------------------------------------------------------------
    None

    --------------------------------------------------------------------------
    Notes
    --------------------------------------------------------------------------
    * Using pyspherematch
    * What is photometric redshift error in photometric catalog? 
    * Details on CMASS redshifts in: https://trac.sdss3.org/wiki/BOSS/clustering/WGCatalogCode
    
    '''

    # read in photometric catalog. See notes for 
    # catalog details 
    photo_file = ''.join([
        '/mount/riachuelo1/hahn/photoz/',
        'DR8LRG_CMASS_phtz_VF.cat'
        ]) 
    photo_ra, photo_dec, photo_z, photo_zerr = np.loadtxt(
            photo_file, 
            unpack=True, 
            skiprows=1, 
            usecols=[0,1,2,3]) 

    # read in spectroscopic catalog. DR12v4 Full 
    # catalog. 
    spec_file = ''.join([
        '/mount/riachuelo1/hahn/data/CMASS/',
        'cmass-dr12v4-N-Reid-full.dat.fits'
        ]) 
    spec = mrdfits(spec_file) 
    spec_ra     = spec.ra
    spec_dec    = spec.dec
    spec_z      = spec.z
    spec_imatch = spec.imatch
    
    # spherematch spectroscopic catalog RA and Decl to 
    # photometric catalog RA and Decl. Afterwards, loop
    # through each of the matches and only keep galaxies
    # that have spectroscopic redshifts. This way we can 
    # calculate the errors in photometric redshift. 
    match_spec, match_photo, d = pysph.spherematch(
            spec_ra, spec_dec, photo_ra, photo_dec, tol=None) 
    
    zspec, zphoto = [], [] 

    for i_spec, i_m in enumerate(match_spec):     
        i_match = spec_imatch[i_m]
        # imatch==0:  Not matched to a redshift.
        # imatch==1:  Matched to a new BOSS redshift.
        # imatch==2:  Matched to an old SDSS redshift.
        # imatch==3:  Fiber collision corrected (using fiberc_zfail)
        # imatch==4:  Star (remove from target list when computing completeness).
        # imatch==5:  Redshift failure corrected (using fiberc_zfail)
        # imatch==6:  Fiber collision, but we cannot find a galaxy* 
        #       redshift for it, this is used for flagging the collision, but if it is not 
        #       overwritten by fiberc_zfail as we find a nearest neighbor (imatch=3), 
        #       it will remain to be 6.
        # imatch==7:  Redshift failure, but we cannot find a galaxy* 
        #       redshift for it, this is initially used for flagging the redshift failures, 
        #       but if it is not overwritten by fiberc_zfail as we find a nearest neighbor 
        #       (imatch=5), it will remain to be 7. 

        if i_match not in (1,2):
            continue 

        if d[i_m] > 10**-10: 
            continue 
        
        zspec.append(spec_z[i_spec])
        zphoto.append(photo_z[match_photo[i_m]]) 
    
    zspec = np.array(zspec)
    zphoto = np.array(zphoto) 

    return [zspec, zphoto]
    
def cmass_deltaz_zspec_zphoto():
    """ Summary statistics of delta z/(1+z) as a function of z.  
    Mu and stddev of delta z/(1+z) for bins of redshift. 
    """

    deltaz_file = ''.join([
        '/mount/riachuelo1/hahn/photoz/',
        'cmass_deltaz_zspec_zphotoz_gauss_stat.dat'
        ]) 

    if not os.path.isfile(deltaz_file):  
        build_cmass_deltaz_zspec_zphoto()

    z_mid, z_low, z_high, mu_deltaz, sigma_deltaz = np.loadtxt(
            deltaz_file, 
            skiprows = 1,
            unpack = True, 
            usecols = [0,1,2,3,4]
            ) 

    return [z_mid, z_low, z_high, mu_deltaz, sigma_deltaz]

def build_cmass_deltaz_zspec_zphoto():
    ''' Calculate summary statistics of delta z/(1+z) as a function of z.  
    Mu and stddev of delta z/(1+z) for bins of redshift. 

    --------------------------------------------------------------------------
    Notes
    --------------------------------------------------------------------------
    * Using pyspherematch
    * Details on CMASS redshifts in: https://trac.sdss3.org/wiki/BOSS/clustering/WGCatalogCode
    * Assuming delta_z/z distribution is gaussian
    
    '''
    
    # matching spectroscopic and photometric redshifts 
    z_spec, z_photo = match_zspec_zphoto_cmass()
    
    # calculate delta_z 
    delta_z = (z_spec - z_photo) / (1.0 + z_spec)
    
    # calculate the average and standard deviations of 
    # delta_z in redshift bins assuming delta_z is a Gaussian
    z_low = np.arange(0.43, 0.7, 0.01)
    z_mid = z_low + 0.005
    z_high = z_low + 0.01

    mu_deltaz, sigma_deltaz = [], []  

    for i_z in xrange(len(z_low)): 

        zbin = np.where(
                (z_spec >= z_low[i_z]) & 
                (z_spec < z_high[i_z])
                ) 
        
        mu_deltaz.append(np.mean(delta_z[zbin]))
        sigma_deltaz.append(np.std(delta_z[zbin]))
    
    # save (zmid, zlow, zhigh, mu_deltaz, sigma_deltaz) to file 
    output_name = ''.join([
        '/mount/riachuelo1/hahn/photoz/',
        'cmass_deltaz_zspec_zphotoz_gauss_stat.dat'
        ]) 
    output_header = 'Columns: zmid, zlow, zhigh, mean(delta z/(1+z)), stddev(delta z/(1+z))'
    np.savetxt(
            output_name, 
            np.c_[
                z_mid, z_low, z_high, mu_deltaz, sigma_deltaz 
                ], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%.5e'], 
            header=output_header, 
            delimiter='\t'
            )
    return None 

if __name__=='__main__': 
    pass
