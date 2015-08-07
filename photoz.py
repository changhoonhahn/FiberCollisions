'''

Photometric Redshift to supplement fiber collision correction 

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
import fibcol_data as fc_data 
import fibcol_dlos as fc_dlos
import fibcol_nbar as fc_nbar
import galaxy_environment as genv
import pyspherematch as pysph
from utility.fitstables import mrdfits
import bovy_plot as bovy

def build_fibcol_assign_photoz(qaplot=False, **cat_corr): 
    ''' Assign photometric redshifts to fiber collided galaxies in mock catalog based on their
    actual redshifts in order to reproduce Delta z / z relation
    
    Parameters
    ----------
    * cat_corr : catalog and correction dictionary for mock catalog 


    Notes
    -----
    * Implemented for Nseries 
    * Appends photo zs at the last column  

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    if catalog['name'] == 'cmass': 
        raise NameError("Only available for mock catalogs NOT CMASS") 
    
    # read data of mock catalog 
    if catalog['name'].lower() == 'nseries': 
        correction['name'] = 'original'
        data_file = fc_data.get_galaxy_data_file('data', **cat_corr) 
        ra, dec, redshift, wfc, zupw, upw_index = np.loadtxt(data_file,  
                unpack=True, usecols=[0, 1, 2, 4, 5, 6]) 
        
        correction['name'] = 'wcompfile'
        wcomp_file = fc_data.get_galaxy_data_file('data', **cat_corr) 
        wcomp = np.loadtxt(wcomp_file, unpack=True, usecols=[0])
    else: 
        raise NotImplementedError('Not yet implemented') 

    fced = np.where(wfc == 0)   # fiber collided galaxies
    #fced = [range(len(redshift))]
    
    # read delta z/(1+z) Gaussian values from table
    photo_dir = '/mount/riachuelo1/hahn/photoz/'
    output_name = ''.join([photo_dir, 
        'cmass_photoz_error_gauss_table.dat']) 
    Dz_zmid, Dz_zlow, Dz_zhigh, mean_Dz, sigma_Dz = np.loadtxt(output_name, 
            skiprows=1, unpack=True, usecols=[0,1,2,3,4])
    
    photoz = np.array([-999. for i in range(len(redshift))])
    for i_fc in fced[0]: 
        closest_zbin = min(range(len(Dz_zmid)), key=lambda i: abs(Dz_zmid[i] - redshift[i_fc])) 
        photoz[i_fc] = redshift[i_fc] - \
                np.random.normal(mean_Dz[closest_zbin], sigma_Dz[closest_zbin]) * \
                (1. + redshift[i_fc])
    
    deltaz_z = (redshift[fced] - photoz[fced]) / (1. + redshift[fced])
    
    if qaplot: 
        prettyplot()                         # set up plot 
        pretty_colors = prettycolors()

        plt.figure(figsize=(8,8))
        
        bovy.scatterplot(redshift[fced], deltaz_z, 
                scatter=True, levels=[0.68, 0.95, 0.997], color=pretty_colors[1], s=3,
                xrange=[0.43, 0.7], yrange=[-0.3, 0.3], 
                xlabel='\mathtt{z_{spec}}', 
                ylabel=r'\mathtt{\frac{|z_{spec} - z_{photo}|}{z_{spec}}}')
        
        fig_file = ''.join(['figure/', 'match_deltaphotoz_specz_', catalog['name'], '.png']) 
        plt.savefig(fig_file, bbox_inches='tight') 
    
    correction['name'] = 'photoz'
    photo_cat_corr = {
            'catalog': catalog, 
            'correction': correction
            }
    photoz_file = fc_data.get_galaxy_data_file('data', **photo_cat_corr)
    print photoz_file 
    
    if catalog['name'].lower() == 'nseries': 
        np.savetxt(photoz_file, 
                np.c_[ra, dec, redshift, wfc, wcomp, zupw, upw_index, photoz], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%i', '%10.5f'], 
                delimiter='\t') 

def match_photoz_specz_cmass(): 
    ''' Compare photometric redshift to spectroscopic redshift. 
    Spherematch photometric and spectroscopic catalogs to 
    determine redshift errors of photmetric redshifts wrt spectroscopic redshifts

    Parameters 
    ----------


    Notes
    -----
    * Using pyspherematch
    * What is photometric redshift error in photometric catalog? 
    * Details on CMASS redshifts in: https://trac.sdss3.org/wiki/BOSS/clustering/WGCatalogCode
    
    '''
    prettyplot()                         # set up plot 
    pretty_colors = prettycolors()

    # read in photometric catalog
    photo_dir = '/mount/riachuelo1/hahn/photoz/'
    photo_file = ''.join([photo_dir, 
        'DR8LRG_CMASS_phtz_VF.cat']) 
    photo = np.loadtxt(photo_file, 
            unpack=True, skiprows=1, usecols=[0,1,2,3]) 
    photo_ra = photo[0] 
    photo_dec = photo[1]
    photo_z = photo[2]
    photo_zerr = photo[3]       # wha?

    # read in spectroscopic catalog
    spec_dir = '/mount/riachuelo1/hahn/data/CMASS/'
    spec_file = ''.join([spec_dir, 
        'cmass-dr12v4-N-Reid-full.dat.fits']) 
    spec = mrdfits(spec_file) 
    spec_ra = spec.ra
    spec_dec = spec.dec
    spec_z = spec.z
    spec_imatch = spec.imatch

    # spherematch upweighted galaxies to downweighted galaxies 
    match_spec, match_photo, d = pysph.spherematch(
            spec_ra, spec_dec, photo_ra, photo_dec, tol=None) 
    
    zspec, zphoto = [], [] 

    for i_spec, i_m in enumerate(match_spec):      # loop through match
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

        if i_match in (1,2):
            pass 
        else:
            continue 

        if d[i_m] > 10**-10: 
            continue 
        
        zspec.append(spec_z[i_spec])
        zphoto.append(photo_z[match_photo[i_m]]) 
    
    bovy.scatterplot(np.array(zspec), np.array(zphoto), 
            scatter=True, color=pretty_colors[1], s=3,
            xrange=[0.0, 1.0], yrange=[0.0, 1.0], 
            xlabel='\mathtt{z_{spec}}', ylabel='\mathtt{z_{photo}}')
    plt.plot(np.arange(0.0, 2.0, 0.1), np.arange(0.0, 2.0, 0.1), c='k', lw=4)
    
    fig_file = ''.join(['figure/', 'match_photoz_specz_cmass.png']) 
    plt.savefig( 
            fig_file, bbox_inches='tight' 
            ) 
    plt.clf() 

    delta_z = np.array(
            [ (zspec[i] - zphoto[i]) / (1. + zspec[i]) for i in range(len(zspec)) ]
            )
    
    bovy.scatterplot(np.array(zspec), delta_z, 
            scatter=True, levels=[0.68, 0.95, 0.997], color=pretty_colors[1], s=3,
            xrange=[0.43, 0.7], yrange=[-0.3, 0.3], 
            xlabel='\mathtt{z_{spec}}', 
            ylabel=r'\mathtt{\frac{|z_{spec} - z_{photo}|}{z_{spec}}}')
    
    fig_file = ''.join(['figure/', 'match_deltaphotoz_specz_cmass.png']) 
    plt.savefig(fig_file, bbox_inches='tight') 

def delta_photoz_specz_cmass():
    ''' Spherematch photometric and spectroscopic catalogs to 
    determine standard deviation redshift errors of photmetric redshifts 
    wrt spectroscopic redshifts

    Parameters 
    ----------

    Notes
    -----
    * Using pyspherematch
    * Details on CMASS redshifts in: https://trac.sdss3.org/wiki/BOSS/clustering/WGCatalogCode
    * Assuming delta_z/z distribution is gaussian
    
    '''
    prettyplot()                         # set up plot 
    pretty_colors = prettycolors()

    # read in photometric catalog
    photo_dir = '/mount/riachuelo1/hahn/photoz/'
    photo_file = ''.join([photo_dir, 
        'DR8LRG_CMASS_phtz_VF.cat']) 
    photo = np.loadtxt(photo_file, 
            unpack=True, skiprows=1, usecols=[0,1,2,3]) 
    photo_ra = photo[0] 
    photo_dec = photo[1]
    photo_z = photo[2]
    photo_zerr = photo[3]       # wha?

    # read in spectroscopic catalog
    spec_dir = '/mount/riachuelo1/hahn/data/CMASS/'
    spec_file = ''.join([spec_dir, 
        'cmass-dr12v4-N-Reid-full.dat.fits']) 
    spec = mrdfits(spec_file) 
    spec_ra = spec.ra
    spec_dec = spec.dec
    spec_z = spec.z
    spec_imatch = spec.imatch

    # spherematch upweighted galaxies to downweighted galaxies 
    match_spec, match_photo, d = pysph.spherematch(
            spec_ra, spec_dec, photo_ra, photo_dec, tol=None) 
    
    zspec, zphoto = [], [] 

    for i_spec, i_m in enumerate(match_spec):      # loop through match
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

        if i_match in (1,2):
            pass 
        else:
            continue 

        if d[i_m] > 10**-10: 
            continue 
        
        zspec.append(spec_z[i_spec])
        zphoto.append(photo_z[match_photo[i_m]]) 
    
    # calculate delta_z 
    delta_z = np.array(
            [ (zspec[i] - zphoto[i]) / (1. + zspec[i]) for i in range(len(zspec)) ]
            )
    
    # calculate sigma assuming delta_z is a Gaussian
    # loop through redshift bins 
    z_low = np.arange(0.43, 0.7, 0.01)
    z_mid = np.arange(0.43, 0.7, 0.01) + 0.005
    z_high = np.arange(0.43, 0.7, 0.01) + 0.01
    mu_deltaz, sigma_deltaz = [], []  
    for i_z in range(len(z_low)): 
        zbin = np.where((zspec > z_low[i_z]) & (zspec < z_high[i_z])) 
        
        mu_deltaz.append(np.mean(delta_z[zbin]))
        sigma_deltaz.append(np.std(delta_z[zbin]))
        #print np.mean(delta_z[zbin]), np.std(delta_z[zbin])
    
    # save to file: zmid, zlow, zhigh, mu_deltaz, sigma_deltaz 
    output_name = ''.join([photo_dir, 
        'cmass_photoz_error_gauss_table.dat']) 
    output_header = 'Columns: zmid, zlow, zhigh, mean(deltaz/(1+z)), stddev(deltaz/(1+z))'
    np.savetxt(output_name, 
            np.c_[
                z_mid, z_low, z_high, mu_deltaz, sigma_deltaz 
                ], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%.5e'], 
            header=output_header, delimiter='\t')

    '''
    fig = plt.figure(1)
    sub = fig.add_subplot(111)
    sub.errorbar( z_low + 0.005, mu_deltaz, yerr=sigma_deltaz, c='k') 
    sub.scatter( z_low + 0.005, mu_deltaz, c='k') 
    sub.set_xlim([0.4, 0.75]) 
    sub.set_xlabel(r'$\mathtt{z_{spec}}$', fontsize=20)
    sub.set_ylabel(r'$\sigma_\mathtt{\Delta z / z}$', fontsize=20)

    fig_file = ''.join(['figure/', 'delta_photoz_specz_cmass.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    fig.clear()
    '''

if __name__=='__main__': 
    #match_photoz_specz_cmass()
    #delta_photoz_specz_cmass()
    for i in range(1, 11): 
        cat_corr = { 
                'catalog': {'name': 'nseries', 'n_mock': i}, 
                'correction': {'name': 'upweight'}
                }
        build_fibcol_assign_photoz(qaplot=True, **cat_corr)
