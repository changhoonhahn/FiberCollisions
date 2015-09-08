'''

Multiprocessing code to calculate corrected data 
in parallel

'''
import multiprocessing as mp

from data import Data 

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

def build_corrdata_multipro(catalog_name, corr_name, n_mocks, Nthreads=8): 
    """ Calculate dLOS for catalogs in parallel using 
    multiprocessing pool

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
    
    pool = mp.Pool(processes=Nthreads)
    mapfn = pool.map
    
    arglist = [ [{
                'catalog': {'name': catalog_name, 'n_mock': i_mock}, 
                'correction': corrdict
                }]
            for i_mock in n_mock_list]
    
    mapfn( build_corrdata_wrapper, [arg for arg in arglist])

    pool.close()
    pool.terminate()
    pool.join() 

    return None 

if __name__=="__main__":
    build_corrdata_multipro('nseries', 'upweight', 84, Nthreads=10)
    build_corrdata_multipro('nseries', 'dlospeak', 84, Nthreads=10)
