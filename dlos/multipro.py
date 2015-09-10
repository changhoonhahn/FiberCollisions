'''

Multiprocessing code to calculate dLOS values 
in parallel

'''
from dlos import Dlos
from dlos_env import DlosEnv
from until.interruptible_pool import InterruptiblePool as Pewl

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

def build_dlosenv_wrapper(params): 
    ''' Wrapper for calculating dLOS dNN values  

    Parameters
    ----------
    params : [ cat_corr, kwargs ] 

    '''

    cat_corr = params[0]
    n_NN = params[1]
    kwargs = {} 
    if len(params) > 2: 
        kwargs = params[2]
    
    dlosclass = DlosEnv(cat_corr, n_NN=n_NN, **kwargs) 
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
    
    pool = Pewl(processes=Nthreads)
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

def build_dlosenv_multipro(catalog_name, n_mocks, Nthreads=8): 
    """ Calculate dLOS for catalogs in parallel using 
    multiprocessing pool
    """
    
    if isinstance(n_mocks, list): 
        n_mock_list = n_mocks
    else:
        n_mock_list = range(1, n_mocks + 1)

    n_NN_list = [1, 3, 5, 7, 10] 

    pool = Pewl(processes=Nthreads)
    mapfn = pool.map
    
    arglist = [] 
    for i_mock in n_mock_list: 
        for n_NN in n_NN_list: 
            arglist.append([
                {
                'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
                'correction': {'name': 'upweight'}
                }, 
                n_NN
                ])
    
    mapfn( build_dlosenv_wrapper, [arg for arg in arglist])

    pool.close()
    pool.terminate()
    pool.join() 

    return None 

if __name__=="__main__":
    #build_dlosenv_multipro('nseries', 1, Nthreads=1)
    build_dlosenv_multipro('nseries', 84, Nthreads=10)
