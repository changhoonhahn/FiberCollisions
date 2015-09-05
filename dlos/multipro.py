'''

Multiprocessing code to calculate dLOS values 
in parallel

'''
import multiprocessing as mp

from dlos import Dlos

def build_dlos_wrapper(params): 
    ''' Wrapper for calculating dLOS values  

    Parameters
    ----------
    params : [ cat_corr, kwargs ] 

    '''

    cat_corr = params[0]
    kwargs = {} 
    if len(params) > 1: 
        kwargs = params[1]
    kwargs['clobber'] = True
    
    dlosclass = Dlos(cat_corr, **kwargs) 
    dlosclass.read()

    return None

def build_dlos_multipro(catalog_name, n_mocks, Nthreads=8): 
    """ Calculate dLOS for catalogs in parallel using 
    multiprocessing pool
    """
    
    if isinstance(n_mocks, list): 
        n_mock_list = n_mocks
    else:
        n_mock_list = range(1, n_mocks + 1)
    
    pool = mp.Pool(processes=Nthreads)
    mapfn = pool.map
    
    arglist = [ [{
                'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
                'correction': {'name': 'upweight'}
                }]
            for i_mock in n_mock_list]
    
    mapfn( build_dlos_wrapper, [arg for arg in arglist])

    pool.close()
    pool.terminate()
    pool.join() 

    return None 

if __name__=="__main__":
    build_dlos_multipro('nseries', 100, Nthreads=4)
