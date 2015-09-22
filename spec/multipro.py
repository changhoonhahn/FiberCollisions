'''

Multiprocessing code to calculate corrected data 
in parallel

'''
from data import Data 
from spec import Spec
from interruptible_pool import InterruptiblePool as Pewl

# --- Wrappers ---
def build_corrdata_wrapper(params): 
    ''' Wrapper for calculating corrected data   

    Parameters
    ----------
    params : [ cat_corr, kwargs ] 

    '''

    cat_corr = params[0]
    kwargs = {} 
    if len(params) > 1: 
        kwargs = params[1]
    kwargs['clobber'] = True
    
    dataclass = Data('data', cat_corr, **kwargs) 
    print dataclass.file_name
    dataclass.build() 

    return None

def build_spec_wrapper(params): 
    """ Wrapper for calculating power/bispectrum
    """
    cat_corr = params[0]
    kwargs = {} 
    if len(params) > 1: 
        kwargs = params[1]

    spectrum = Spec('pk', cat_corr, **kwargs)
    print spectrum.file()
    spectrum.build()

    return None 

# --- Multiprocessing --- 
def build_multipro(type, catalog_name, corr_name, n_mocks, Nthreads=8, **kwargs): 
    """ Calculate dLOS for catalogs in parallel using interruptible
    pool, which is multiprocessing pool that allows for interrputions

    Parameters
    ----------
    catalog_name : Name of catalog 
    corr_name : Name of correction
    n_mocks : Number of mock catalogs to calculate 
    Nthreads : Number of CPUs to use 

    """
    
    if isinstance(n_mocks, list): 
        n_mock_list = n_mocks
    else:
        n_mock_list = range(1, n_mocks + 1)

    corrdict = {} 
    if catalog_name == 'nseries':

        corrdict['name'] = corr_name

        if 'dlospeak' in corr_name: 
            # hardcoded values for bestfit dlos peak
            # parameters
            corrdict['fit'] = 'gauss'
            corrdict['sigma'] = 3.9
            corrdict['fpeak'] = 0.68 

        if 'env' in corr_name: 
            # hardcoded values for galaxy environment
            # parameters
            corrdict['n_NN'] = 5

        if 'photoz' in corr_name: 

            corrdict['d_photoz_tail_cut'] = 200 
    
    pool = Pewl(processes=Nthreads)
    mapfn = pool.map
    
    arglist = [ [{
                'catalog': {'name': catalog_name, 'n_mock': i_mock}, 
                'correction': corrdict
                }, kwargs]
            for i_mock in n_mock_list]
    
    if type == 'data': 
        mapfn( build_corrdata_wrapper, [arg for arg in arglist])
    elif type == 'spec': 
        mapfn( build_spec_wrapper, [arg for arg in arglist])

    pool.close()
    pool.terminate()
    pool.join() 

    return None 

if __name__=="__main__":
    build_multipro('spec', 'nseries', 'dlospeakphotoz', 20, Nthreads=10, clobber=True)
    #build_multipro('spec', 'nseries', 'true', 20, Nthreads=10, clobber=True)
    #build_multipro('spec', 'nseries', 'upweight', 20, Nthreads=10, clobber=True)
    #build_multipro('spec', 'nseries', 'dlospeakenv', 20, Nthreads=10, clobber=True)
    #build_multipro('data', 'nseries', 'photoz', 84, Nthreads=10)
    #build_multipro('spec', 'nseries', 'true', 84, Nthreads=10)
    #build_multipro('spec', 'nseries', 'upweight', 84, Nthreads=10)
