'''

Multiprocessing code to calculate corrected data or 
power/bispectrum in parallel!

'''
from corr_spec import CorrSpec
from corr_corrdata import CorrCorrData 
from Spectrum.interruptible_pool import InterruptiblePool as Pewl

def build_multipro(type, catalog_name, corr_name, n_mocks, Nthreads=8, ell=2, Ngrid=360, **kwargs): 
    ''' 
    Calculate dLOS for catalogs in parallel using interruptible
    pool, which is multiprocessing pool that allows for interrputions

    Parameters
    ----------
    catalog_name : Name of catalog 
    corr_name : Name of correction
    n_mocks : Number of mock catalogs to calculate 
    Nthreads : Number of CPUs to use 

    '''
    
    if isinstance(n_mocks, list): 
        n_mock_list = n_mocks
    else:
        n_mock_list = range(1, n_mocks + 1)

    corrdict = {} 
    if catalog_name == 'nseries':               # Nseries
        if isinstance(corr_name, dict): 
            corrdict = corr_name
        else:
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

                corrdict['d_photoz_tail_cut'] = 15 

            if corr_name == 'fourier_tophat': 
                corrdict['fs'] = 1.0 
                corrdict['rc'] = 0.43 
                corrdict['k_fit'] = 0.7 
                corrdict['k_fixed'] = 0.84

    elif catalog_name == 'qpm': 
        if isinstance(corr_name, dict): 
            corrdict = corr_name
        else:
            corrdict['name'] = corr_name

    elif catalog_name == 'cmass': 
        if isinstance(corr_name, dict): 
            corrdict = corr_name
        else:
            corrdict['name'] = corr_name

    else: 
        raise NotImplementedError
    
    if type == 'bk':
        arglist = [ 
                [{ 
                    'catalog': {'name': catalog_name, 'n_mock': i_mock}, 
                    'correction': corrdict, 
                    'spec': {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid} 
                    }, kwargs]
                for i_mock in n_mock_list]
    elif type == 'pk':
        arglist = [ [{
                    'catalog': {'name': catalog_name, 'n_mock': i_mock}, 
                    'correction': corrdict, 
                    'spec': { 'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell} 
                    }, ell, kwargs]
                for i_mock in n_mock_list
                ]
    elif type == 'data': 
        arglist = [ [{
                    'catalog': {'name': catalog_name, 'n_mock': i_mock}, 
                    'correction': corrdict, 
                    'spec': { 'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid, 'ell': ell} 
                    }, kwargs]
                for i_mock in n_mock_list]
    else: 
        raise NameError
    
    if Nthreads > 1: 
        pool = Pewl(processes=Nthreads)
        mapfn = pool.map
    
        if type == 'data': 
            mapfn( build_corrdata_wrapper, [arg for arg in arglist])
        elif type == 'pk': 
            mapfn( build_pk_wrapper, [arg for arg in arglist])
        elif type == 'bk': 
            mapfn( build_bk_wrapper, [arg for arg in arglist])

        pool.close()
        pool.terminate()
        pool.join() 
    else: 
        for arg in arglist: 
            if type == 'data': 
                build_corrdata_wrapper(arg)
            elif type == 'pk': 
                build_pk_wrapper(arg)
            elif type == 'bk': 
                build_bk_wrapper(arg)

    return None 

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
    
    dataclass = CorrCorrData(cat_corr, **kwargs) 
    print dataclass.file_name
    dataclass.build() 

    return None

def build_pk_wrapper(params): 
    ''' Wrapper for calculating power/bispectrum
    '''
    cat_corr = params[0]
    ell = params[1]
    kwargs = {} 
    if len(params) > 2: 
        kwargs = params[2]

    spectrum = CorrSpec('pk', cat_corr, **kwargs)
    print spectrum.file()
    spectrum.build()

    return None 

def build_bk_wrapper(params): 
    ''' Wrapper for calculating power/bispectrum
    '''
    cat_corr = params[0]
    kwargs = {} 
    if len(params) > 2: 
        kwargs = params[1]

    spectrum = CorrSpec('bk', cat_corr, **kwargs)
    print spectrum.file()
    spectrum.build()

    return None 

# --- Multiprocessing --- 
if __name__=="__main__":
    #build_multipro('data', 'cmass', 'upweight', [1], Nthreads=1, clobber=True)
    build_multipro('pk', 'qpm', 'true', range(1,11), ell=2, Nthreads=1, clobber=True, Ngrid=480)
    build_multipro('pk', 'qpm', 'upweight', range(1,11), ell=2, Nthreads=1, clobber=True, Ngrid=480)
    #build_multipro('data', 'qpm', 'upweight', range(1,51), Nthreads=5, clobber=True)
    #build_multipro('bk', 'nseries', 'upweight', range(2, 85), Nthreads=1, clobber=True, Ngrid=360)
    #build_multipro('bk', 'nseries', 'true', range(2, 85), Nthreads=1, clobber=True, Ngrid=360)
    #build_multipro('pk', 'nseries', 'fourier_tophat', [1], Nthreads=1, clobber=True, Ngrid=960)
