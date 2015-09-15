'''
    
QAplots for fiber collision correction methods 
implemented in corrections module

Author(s): ChangHoon Hahn 

'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# --- Local ---
import bovy_plot as bovy
from spec.data import Data
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 

def test_photoz(catalog_name): 
    """ Test that the assigned photometric redshifts in photoz correction 
    method reproduce the CMASS photometric redshifts. 
    """

    catdict = {'name': catalog_name, 'n_mock': 1}
    corrdict = {'name': 'photoz'}
    cat_corr = {'catalog': catdict, 'correction': corrdict}

    dataclass = Data('data', cat_corr) 
    dataclass.read()
    
    # only consider galaxies with photometric redshifts 
    hasphotoz = np.where( 
            dataclass.photoz > -999.
            ) 

    deltaz_z = (dataclass.z[hasphotoz] - dataclass.photoz[hasphotoz]) / (1.0 + dataclass.z[hasphotoz])

    prettyplot()                        
    pretty_colors = prettycolors()

    plt.figure(figsize=(8,8))
            
    bovy.scatterplot(
            dataclass.z[hasphotoz], 
            deltaz_z, 
            scatter = True, 
            levels = [0.68, 0.95, 0.997], 
            color = pretty_colors[1], 
            s = 3,
            xrange = [0.43, 0.7], 
            yrange = [-0.3, 0.3], 
            xlabel = '\mathtt{z_{spec}}', 
            ylabel = r'\mathtt{\frac{|z_{spec} - z_{photo}|}{1 + z_{spec}}}' 
            )
            
    # save figure to file  
    fig_dir = '/home/users/hahn/powercode/FiberCollisions/figure/'
    fig_file = ''.join([
        fig_dir, 
        catalog_name, '_deltaz_zspec_zphoto.png'
        ]) 
    plt.savefig(fig_file, bbox_inches='tight') 

    plt.close()
    
    return None

if __name__=="__main__": 
    test_photoz('nseries')
