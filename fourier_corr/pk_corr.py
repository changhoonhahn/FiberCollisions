'''

fourier correction for P_l(k). 
Interfaces with spec modules in order to actually calculate P_l(k) for mocks/data.

'''
from spec.spec import Spec
import pk_extrap

def test_pk_extrap_scatter(ell, n_mocks, Ngrid=960, k_max=0.25, k_fixed=0.6, **kwargs): 
    '''
    test the scatter in pk extrapolation
    '''
    alphas, ns = [], []
    for i_mock in range(1, n_mocks+1):
        # default cat-corr for Nseries
        cat_corr = {
                'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
                'correction': {'name': 'true'}
                }
        spec_i = Spec('pk', cat_corr, ell=ell, Ngrid=Ngrid)
        spec_i.read()

        spec_i_spec = getattr(spec_i, 'p'+str(ell)+'k')

        alpha_i, n_i = pk_extrap.pk_powerlaw_bestfit(spec_i.k, spec_i_spec, k_max=k_max, k_fixed=k_fixed)

        alphas.append(alpha_i)
        ns.append(n_i)

    alphas = np.array(alphas) 
    ns = np.array(ns)

    print 'alphas'
    print 'mean ', np.mean(alphas)
    print 'min ', np.min(alphas), ' max ', np.max(alphas)
    print 'stddev ', np.std(alphas)

    print 'ns'
    print 'mean ', np.mean(ns)
    print 'min ', np.min(ns), ' max ', np.max(ns)
    print 'stddev ', np.std(ns)

if __name__=='__main__':
    test_pk_extrap_scatter(4, 20, Ngrid=960, k_max=0.7, k_fixed=0.8) 

