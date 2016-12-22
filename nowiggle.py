'''

No wiggles stuff

'''
import numpy as np 

from corr_spec.corr_spec import CorrSpec

# --- plotting ---
from matplotlib import gridspec
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot 
from ChangTools.plotting import prettycolors 


def Build_Nowiggle(correction='true', Ngrid=960): 
    ''' Build no wiggle power spectrums
    '''
    nowig_cat_corr = { 
            'catalog': {'name': 'bigmd', 'version': 'nowiggles'}, 
            'correction': {'name': correction}, 
            'spec':{'P0': 20000,  'Lbox': 3600, 'Ngrid': Ngrid, 'ell': 0}
            }
    nowig_spec = CorrSpec('pk', nowig_cat_corr)
    nowig_spec.build() 
    return None


def Plot_Nowiggle(correction='true'): 
    '''
    '''
    cat_corr = { 
            'catalog': {'name': 'bigmd'}, 
            'correction': {'name': correction}, 
            'spec':{'P0': 20000,  'Lbox': 3600, 'Ngrid':960, 'ell': 0}
            }
    spec = CorrSpec('pk', cat_corr)
    print 'With wiggles ... ' 
    print spec.file_name
    spec.read() 

    nowig_cat_corr = { 
            'catalog': {'name': 'bigmd', 'version': 'nowiggles'}, 
            'correction': {'name': correction}, 
            'spec':{'P0': 20000,  'Lbox': 3600, 'Ngrid':960, 'ell': 0}
            }
    nowig_spec = CorrSpec('pk', nowig_cat_corr)
    print 'WithOUT wiggles ... ' 
    print nowig_spec.file_name
    nowig_spec.read() 

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(1, figsize=(8,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    spec_sub = plt.subplot(gs[0])#fig.add_subplot(121)
    ratio_sub = plt.subplot(gs[1]) # fig.add_subplot(122)

    spec_sub.plot(spec.k, spec.p0k, 
            c='k', lw=3, label='With Wiggles')
    spec_sub.plot(nowig_spec.k, nowig_spec.p0k, 
            c=pretty_colors[3], lw=3, label='Without Wiggles')

    spec_sub.set_xlabel('k (h/Mpc)', fontsize=24)
    spec_sub.set_xlim([10**-2,10**0])
    spec_sub.set_xscale('log')
    spec_sub.set_xticklabels([])
    spec_sub.set_yticklabels([])
    spec_sub.set_yticks([])
    spec_sub.set_ylabel(r'$\mathtt{P_0(k)}$', fontsize=24)
    spec_sub.set_ylim([10**2,10**5.5])
    spec_sub.set_yscale('log')
    spec_sub.legend()
    
    if np.array_equal(spec.k, nowig_spec.k): 
        ratio_sub.plot(spec.k, nowig_spec.p0k/spec.p0k, c='k', lw=3)
    else: 
        raise ValueError
    ratio_sub.set_xlabel('k (h/Mpc)', fontsize=24)
    ratio_sub.set_xlim([10**-2,10**0])
    ratio_sub.set_xscale('log')
    ratio_sub.set_ylim([0.75, 1.25])
    ratio_sub.set_ylabel(r'$\mathtt{P^{no\;wiggles}_0/P^{wiggles}_0}$', fontsize=24)
    
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/', 
        'nowiggle_P0k_comparison', 
        '.', correction, '.png']) 

    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close() 


def Clean_Pk_BigMD(correction='true'): 
    ''' Produce cleaned up ASCII files without all the extra columns
    '''
    cat_corr = { 
            'catalog': {'name': 'bigmd'}, 
            'correction': {'name': correction}, 
            'spec':{'P0': 20000,  'Lbox': 3600, 'Ngrid':960, 'ell': 0}
            }
    spec = CorrSpec('pk', cat_corr)
    print 'With wiggles ... ' 
    print spec.file_name
    spec.read() 

    nowig_cat_corr = { 
            'catalog': {'name': 'bigmd', 'version': 'nowiggles'}, 
            'correction': {'name': correction}, 
            'spec':{'P0': 20000,  'Lbox': 3600, 'Ngrid':960, 'ell': 0}
            }
    nowig_spec = CorrSpec('pk', nowig_cat_corr)
    print 'WithOUT wiggles ... ' 
    print nowig_spec.file_name
    nowig_spec.read() 

    # write the files into ASCII format. Simple k Pk 
    
    # With Wiggles
    np.savetxt(
            '/mount/riachuelo1/hahn/nowiggles/'+spec.file_name.split('/')[-1], 
            (np.vstack(np.array([spec.k, spec.p0k]))).T, 
            fmt=['%.6e', '%.6e'],
            delimiter='\t', 
            header='k, P0k'
            ) 
    # Without wiggles 
    np.savetxt(
            '/mount/riachuelo1/hahn/nowiggles/'+nowig_spec.file_name.split('/')[-1], 
            (np.vstack(np.array([nowig_spec.k, nowig_spec.p0k]))).T, 
            fmt=['%.6e', '%.6e'],
            delimiter='\t', 
            header='k, P0k'
            ) 
    return None


if __name__=='__main__': 
    Clean_Pk_BigMD(correction='true')
    Clean_Pk_BigMD(correction='upweight')
    #Plot_Nowiggle(correction='true') 
    #Plot_Nowiggle(correction='upweight') 
