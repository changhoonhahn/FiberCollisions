'''



'''
import numpy as np 

from corr_spec.corr_spec import CorrSpec
from corr_spec.corr_average import CorrAvgSpec as AvgSpec

# --- plotting ---
from matplotlib import gridspec
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot 
from ChangTools.plotting import prettycolors 


def Plot_QSO(version, ell=0): 
    '''
    '''
    eboss_dict = { 
            'catalog': {'name': 'qso_bigmd', 'version': 'ebossv1.5'}, 
            'correction': {'name': 'true'}, 
            'spec':{'P0': 6000,  'Lbox': 8000, 'Ngrid':960, 'ell': ell}
            }
    eboss_spec = AvgSpec(1, 'pk', eboss_dict, ell=ell)
    print 'eBOSS QSOs ... ' 
    print eboss_spec.file_name
    eboss_spec.read(rebin=12) 

    bigmd_dict = { 
            'catalog': {'name': 'qso_bigmd', 'version': version}, 
            'correction': {'name': 'true'}, 
            'spec':{'P0': 6000,  'Lbox': 8000, 'Ngrid':960, 'ell': ell}
            }
    bigmd_spec = AvgSpec(1, 'pk', bigmd_dict, ell=ell)
    print 'BigMD QSOs ... ' 
    print bigmd_spec.file_name
    bigmd_spec.read(rebin=12) 
    
    #if ell == 0: 
    #    # read in the jack knife P(k)
    #    jack_Pk = [] 
    #    for n_jack in range(1, 51): 
    #        jack_dict = { 
    #                'catalog': {'name': 'qso_bigmd', 'version': 'jackknife'+str(n_jack)}, 
    #                'correction': {'name': 'true'}, 
    #                'spec':{'P0': 20000,  'Lbox': 8000, 'Ngrid':960, 'ell': 0}
    #                }
    #        jack_spec = AvgSpec(1, 'pk', jack_dict, ell=0)
    #        jack_spec.read(rebin=12) 
    #        print jack_spec.file_name
    #        jack_Pk.append(jack_spec.p0k) 
    #    jack_Pk_avg = np.mean(np.array(jack_Pk), axis=0)
    #    
    #    jack_Pk_err = np.zeros(len(jack_Pk_avg))
    #    for i_jack in range(50): 
    #        jack_Pk_err += (jack_Pk[i_jack] - jack_Pk_avg)**2
    #    jack_Pk_err *= 49./50.
    #    jack_Pk_err = np.sqrt(jack_Pk_err)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(1, figsize=(8,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    spec_sub = plt.subplot(gs[0])#fig.add_subplot(121)
    ratio_sub = plt.subplot(gs[1]) # fig.add_subplot(122)
    
    if ell == 0: 
        #spec_sub.fill_between(eboss_spec.k, 
        #        eboss_spec.p0k - jack_Pk_err, eboss_spec.p0k + jack_Pk_err, 
        #        color='k', alpha=0.5, lw=0, label='Jackknife') 
        spec_sub.plot(eboss_spec.k, eboss_spec.p0k, 
                c='k', lw=1, label='eBOSS')
        spec_sub.scatter(bigmd_spec.k, bigmd_spec.p0k, 
                c=pretty_colors[3], lw=0, s=15, label='BigMD QSO')
    else: 
        spec_sub.plot(eboss_spec.k, eboss_spec.p2k, 
                c='k', lw=1, label='eBOSS')
        spec_sub.scatter(bigmd_spec.k, bigmd_spec.p2k, 
                c=pretty_colors[3], lw=0, s=15, label='BigMD QSO')

    spec_sub.set_xlabel('k (h/Mpc)', fontsize=24)
    spec_sub.set_xlim([10**-2,10**0])
    spec_sub.set_xscale('log')
    spec_sub.set_xticklabels([])
    spec_sub.set_yticklabels([])
    spec_sub.set_yticks([])
    spec_sub.set_ylabel(r'$\mathtt{P_0(k)}$', fontsize=24)
    if ell == 0: 
        spec_sub.set_ylim([10**3,10**6])
    else: 
        spec_sub.set_ylim([10**2,10**5])
    spec_sub.set_yscale('log')
    spec_sub.legend(scatterpoints=1)
    
    if np.array_equal(eboss_spec.k, bigmd_spec.k): 
        if ell == 0: 
            ratio_sub.plot(eboss_spec.k, bigmd_spec.p0k/eboss_spec.p0k, c=pretty_colors[3], lw=3)
        else: 
            ratio_sub.plot(eboss_spec.k, bigmd_spec.p2k/eboss_spec.p2k, c=pretty_colors[3], lw=3)
    else: 
        raise ValueError
    ratio_sub.set_xlabel('k (h/Mpc)', fontsize=24)
    ratio_sub.set_xlim([10**-2,10**0])
    ratio_sub.set_xscale('log')
    if ell == 0: 
        ratio_sub.set_ylim([0.5, 1.5])
    else: 
        ratio_sub.set_ylim([0., 2.])

    ratio_sub.set_ylabel(r'$\mathtt{P^{bigmd}_0/P^{eboss}_0}$', fontsize=24)
    
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/', 
        'qsoBigMD_P', str(ell), 'k_', version, '_comparison.png']) 
    fig.savefig(fig_file, bbox_inches="tight", dpi=150)
    plt.close() 


def rebin_Pk(version, ell=0): 
    '''
    '''
    eboss_dict = { 
            'catalog': {'name': 'qso_bigmd', 'version': version}, 
            'correction': {'name': 'true'}, 
            'spec':{'P0': 20000,  'Lbox': 8000, 'Ngrid':960, 'ell': ell}
            }
    eboss_spec = AvgSpec(1, 'pk', eboss_dict, ell=ell)
    print 'eBOSS QSOs ... ' 
    print eboss_spec.file_name
    eboss_spec.read(rebin=12) 


if __name__=='__main__':
    #rebin_Pk('v2-z', ell=0)
    #rebin_Pk('v2-z', ell=2)
    #rebin_Pk('v2-nsat', ell=0)
    #rebin_Pk('v2-nsat', ell=2)
    Plot_QSO('v2', ell=0)
    #Plot_QSO('v2-z', ell=0)
    #Plot_QSO('v2-nsat', ell=0)
    #Plot_QSO('v2', ell=0)
    #Plot_QSO('v2', ell=2)
