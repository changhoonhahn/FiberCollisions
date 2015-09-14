"""

Plotting for Photo Z analysis 

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# --- Local ---
import bovy_plot as bovy
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from photoz.photoz import match_zspec_zphoto_cmass
from photoz.photoz import cmass_deltaz_zspec_zphoto 

def cmass_zspec_zphoto(): 
    """ Compare matching spectroscopic redshifts and photometric redshifts
    of the CMASS catalog. 
    """ 
    prettyplot()
    pretty_colors = prettycolors()

    # matching spectroscopic and photometric redshifts 
    z_spec, z_photo = match_zspec_zphoto_cmass()
    
    # scatter plot of z_spec versus z_photo
    bovy.scatterplot(
            z_spec, 
            z_photo, 
            scatter=True, 
            color=pretty_colors[1], 
            s=3,
            xrange=[0.0, 1.0], 
            yrange=[0.0, 1.0], 
            xlabel='\mathtt{z_{spec}}', 
            ylabel='\mathtt{z_{photo}}'
            )
    plt.plot(np.arange(0.0, 2.0, 0.1), np.arange(0.0, 2.0, 0.1), c='k', lw=4)   # y = x 
    
    fig_dir = '/home/users/hahn/powercode/FiberCollisions/figure/'
    fig_file = ''.join([
        fig_dir, 
        'cmass_zspec_zphoto.png'
        ]) 
    plt.savefig( 
            fig_file, 
            bbox_inches='tight' 
            ) 
    plt.clf() 

    # scatter plot of delta z/(1+z) vs z 

    delta_z = (z_spec - z_photo) / (1.0 + z_spec)
    
    bovy.scatterplot(
            z_spec, 
            delta_z, 
            scatter=True, 
            levels=[0.68, 0.95, 0.997],
            color=pretty_colors[1], 
            s=3,
            xrange=[0.43, 0.7], 
            yrange=[-0.3, 0.3], 
            xlabel='\mathtt{z_{spec}}', 
            ylabel=r'\mathtt{\frac{|z_{spec} - z_{photo}|}{1\;+\;z_{spec}}}'
            )
    
    fig_file = ''.join([
        fig_dir, 
        'cmass_delta_zphoto_zspec.png'
        ]) 
    plt.savefig(
            fig_file, 
            bbox_inches='tight'
            ) 
    return None 

def cmass_deltaz_zspec_zphoto_test(): 
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
    
    # matching spectroscopic and photometric redshifts 
    z_spec, z_photo = match_zspec_zphoto_cmass()
    
    delta_z = (z_spec - z_photo) / (1.0 + z_spec)
    
    bovy.scatterplot(
            z_spec, 
            delta_z, 
            scatter=True, 
            levels=[0.68, 0.95, 0.997],
            color=pretty_colors[1], 
            s=3,
            xrange=[0.43, 0.7], 
            yrange=[-0.3, 0.3], 
            xlabel='\mathtt{z_{spec}}', 
            ylabel=r'\mathtt{\frac{|z_{spec} - z_{photo}|}{1\;+\;z_{spec}}}'
            )
    z_mid, z_low, z_high, mu_deltaz, sigma_deltaz = cmass_deltaz_zspec_zphoto()

    plt.errorbar(
            z_mid, 
            mu_deltaz, 
            yerr=sigma_deltaz, 
            c=pretty_colors[3]
            ) 
    plt.scatter(
            z_mid, 
            mu_deltaz, 
            c=pretty_colors[3]
            )

    fig_dir = '/home/users/hahn/powercode/FiberCollisions/figure/'
    fig_file = ''.join([
        fig_dir, 
        'cmass_deltaz_zspec_zphoto.png'
        ]) 
    plt.savefig(
            fig_file, 
            bbox_inches='tight'
            ) 
    plt.close()

if __name__=="__main__": 
    cmass_deltaz_zspec_zphoto_test()
