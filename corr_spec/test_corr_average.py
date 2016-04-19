'''

Test average spec module 

'''
import numpy as np 
import scipy as sp

from corr_average import CorrAvgSpec as AvgSpec
# plotting ----
from defutility.plotting import prettyplot 
from defutility.plotting import prettycolors 

def cosmic_variance_kbin(ell):
    '''
    Tests the correlation of the cosmic variance by calculating it 
    with different k_bins 
    '''
    specdict = {
            'P0': 20000,
            'Lbox': 3600, 
            'Ngrid': 960, 
            'ell': ell 
            }
    true_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'true'},
            'spec': specdict
            }
    true_spec = AvgSpec(20, 'pk', true_cat_corr)
    orig_stddev = len(true_spec.stddev())
    orig_k = true_spec.k 
    orig_avg_spec = true_spec.avg_spec
    rebin_stddev = len(true_spec._stddev_kbin('double'))
    rebin_k = true_spec.k 
    rebin_avg_spec = true_spec.avg_spec

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=[14,8])
    sub = fig.add_subplot(111)
    
    orig_spec_lower = orig_avg_spec - orig_stddev
    orig_spec_upper = orig_avg_spec + orig_stddev
    if orig_spec_lower.min() < 0.: 
        orig_spec_lower[np.where(orig_spec_lower < 0.)] = 10**-3 
        orig_spec_upper[np.where(orig_spec_upper < 0.)] = 10**-3 

    rebin_spec_lower = rebin_avg_spec - rebin_stddev
    rebin_spec_upper = rebin_avg_spec + rebin_stddev
    if rebin_spec_lower.min() < 0.: 
        rebin_spec_lower[np.where(rebin_spec_lower < 0.)] = 10**-3 
        rebin_spec_upper[np.where(rebin_spec_upper < 0.)] = 10**-3 

    sub.fill_between(orig_k, 
            orig_spec_lower,
            orig_spec_upper, 
            color = pretty_colors[0], 
            label = 'Original k binning')

    sub.fill_between(rebin_k, 
            rebin_spec_lower,  
            rebin_spec_upper, 
            color = pretty_colors[2], 
            label = 'Rebinned k binning')

    sub.set_xscale('log')
    sub.set_yscale('log')
    sub.set_xlim([10**-3, 10**0])
    sub.set_ylim([10**2, 10**5])

    sub.set_xlabel(r"$\mathtt{k}\;\;(\mathtt{Mpc}/h)$", fontsize=30)
    sub.set_ylabel(r"$\mathtt{P_{"+str(ell)+"}(k)}$", fontsize=30)
    
    fig.savefig(
            ''.join(['figure/'
                'P', str(ell), 'k.', 
                'cosmic_variance.kbinning_test.png']), 
            bbox_inches='tight') 
    plt.close()


if __name__=="__main__":
    cosmic_variance_kbin(0)
    cosmic_variance_kbin(2)
