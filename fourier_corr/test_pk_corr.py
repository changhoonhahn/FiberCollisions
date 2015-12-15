import pk_extrap
import fourier_corr
from corr_spec.corr_spec import CorrSpec
from corr_spec.corr_average import CorrAvgSpec

from defutility.plotting import prettyplot 
from defutility.plotting import prettycolors 

def test_pk_extrap_scatter(ell, n_mocks, Ngrid=960, k_fit=0.25, k_fixed=0.6, **kwargs): 
    '''
    test the scatter in pk extrapolation
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure()
    sub = fig.add_subplot(111)

    alphas, ns = [], []
    for i_mock in range(1, n_mocks+1):
        # default cat-corr for Nseries
        cat_corr = {
                'catalog': {'name': 'nseriesbox', 'n_mock': i_mock}, 
                'correction': {'name': 'true'}
                }
        spec_i = CorrSpec('pk', cat_corr, ell=ell, Ngrid=Ngrid)
        print spec_i.file_name
        spec_i.read()

        spec_i_spec = getattr(spec_i, 'p'+str(ell)+'k')

        alpha_i, n_i = pk_extrap.pk_powerlaw_bestfit(spec_i.k, spec_i_spec, k_fit=k_fit, k_fixed=k_fixed)

        if i_mock == 1: 
            pk_label = 'Data'
            pk_extrap_label = 'Extrap.'
        else: 
            pk_label = None
            pk_extrap_label = None

        sub.plot(spec_i.k, np.abs(spec_i_spec), lw=2, c=pretty_colors[0], label=pk_label)
        sub.plot(
                np.arange(k_fit, 10, 0.01), 
                np.abs(pk_extrap.pk_powerlaw(np.arange(k_fit, 10, 0.01), [alpha_i, n_i], k_fixed=k_fixed)), 
                ls='--', lw=2, c=pretty_colors[1], label = pk_extrap_label
                )

        alphas.append(alpha_i)
        ns.append(n_i)

    alphas = np.array(alphas) 
    ns = np.array(ns)

    # x-axis
    sub.set_xlabel(r'$\mathtt{k}$', fontsize=30)
    sub.set_xscale('log')
    sub.set_xlim([0.1, 10.])
    # y-axis
    sub.set_ylabel(r'$\mathtt{P_'+str(ell)+'(k)}$', fontsize=30)
    sub.set_yscale('log')
    sub.set_ylim([10**0, 10**5])
    
    sum_stat = ''.join([
        r"$\bar{\alpha}, \sigma_{\alpha} = ",
        str(round(np.mean(alphas), 1)), ",\;", str(round(np.std(alphas),1)), "$", 
        '\n', 
        r"$\bar{n}, \sigma_{n} = ",
        str(round(np.mean(ns),2)), ",\;", str(round(np.std(ns),2)), "$"
        ])
    sub.text(1.0, 10**4.25, sum_stat, fontsize=20)
    sub.legend(loc='lower right')
    plt.show()

def test_delPk_corr_scatter(ell, n_mocks, Ngrid=960, k_fit=0.25, k_fixed=0.6, fs=1.0, rc=0.43, **kwargs): 
    '''
    test the scatter in del P(k)^corr from extrapolation of each individual mock catalog realizations
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure(figsize=(14,8))
    sub = fig.add_subplot(111)

    for i_mock in range(1, n_mocks+1):
        # default cat-corr for Nseries
        cat_corr = {
                'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
                'correction': {'name': 'true'}
                }
        spec_i = Spec('pk', cat_corr, ell=4, Ngrid=Ngrid)
        spec_i.read()

        specs_i = [getattr(spec_i, 'p0k'), getattr(spec_i, 'p2k'), getattr(spec_i, 'p4k')]
        
        # extrapolation parameters for P0, P2, P4
        extrap_pars = [] 
        for i_ell, ellp in enumerate([0,2,4]):
            extrap_pars.append(
                    pk_extrap.pk_powerlaw_bestfit(spec_i.k, specs_i[i_ell], k_fit=k_fit, k_fixed=k_fixed)
                    )
        print extrap_pars

        # calculate delP
        corrdelP = fourier_corr.delP_corr(spec_i.k, specs_i, ell, fs=fs, rc=rc, extrap_params=extrap_pars, k_fixed=k_fixed)

        if i_mock == 1: 
            delpk_label = r"$\mathtt{\Delta P^{corr}}$ (Extrap.)"
        else: 
            delpk_label = None

        sub.plot(spec_i.k, corrdelP, lw=2, c=pretty_colors[i_mock % 20], label=delpk_label)
    
    # upw P_l(k) - true P_l(k0
    true_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'true'}, 
            }
    true_avg_spec = AvgSpec(n_mocks, 'pk', true_cat_corr, ell=ell, Ngrid=Ngrid)
    true_avg_spec.read()
    upweight_cat_corr = true_cat_corr.copy()
    upweight_cat_corr['correction']['name'] = 'upweight'
    upweight_avg_spec = AvgSpec(n_mocks, 'pk', upweight_cat_corr, ell=ell, Ngrid=Ngrid)
    upweight_avg_spec.read()

    sub.plot(
            true_avg_spec.k, 
            getattr(upweight_avg_spec, 'p'+str(ell)+'k') - getattr(true_avg_spec, 'p'+str(ell)+'k'),
            c = 'k', 
            lw = 4,
            ls = '--', 
            label = r"$\mathtt{P^{upw}(k) - P^{true}(k)}$"
            )

    # x-axis
    sub.set_xlabel(r'$\mathtt{k}$', fontsize=30)
    sub.set_xscale('log')
    sub.set_xlim([10**-3, 10**0.])
    # y-axis
    sub.set_ylabel(r"$\mathtt{\Delta P_{"+str(ell)+"}(k)}$", fontsize=30)
    sub.set_ylim([-50., 300.])
    
    sub.legend(loc='upper left')
    fig_file = ''.join([
        'figure/', 
        'qaplot_delP_corr_extrap_scatter.', 
        str(n_mocks), 'mocks.', 
        'kfit', str(round(k_fit,2)), '.', 
        'kfixed', str(round(k_fixed,2)),
        '.png'
        ])
    fig.savefig(fig_file, bbox_inches="tight")
    plt.close()

if __name__=='__main__':
    #for k_fit in np.arange(3.5, 4.5, 0.25):
    test_pk_extrap_scatter(6, 7, Ngrid=960, k_fit=4.0, k_fixed=4.3)
        #test_delPk_corr_scatter(2, 20, Ngrid=960, k_fit=k_fit, k_fixed=0.837) 
