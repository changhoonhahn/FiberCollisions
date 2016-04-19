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
            
            if 'dlospeak' in corr_name: 
                # hardcoded values for bestfit dlos peak
                # parameters
                corrdict['fit'] = 'gauss'
                corrdict['sigma'] = 4.4 
                corrdict['fpeak'] = 0.62
    
    elif catalog_name == 'bigmd': 
        if isinstance(corr_name, dict): 
            corrdict = corr_name
        else:
            corrdict['name'] = corr_name
            
            if 'dlospeak' in corr_name: 
                # hardcoded values for bestfit dlos peak
                # parameters
                corrdict['fit'] = 'gauss'
                corrdict['sigma'] = 5.5 
                corrdict['fpeak'] = 0.6

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
    '''
    # Done Harmattan 9:07AM
    build_multipro('pk', 'nseries', 'dlospeakknown', range(8,11), Nthreads=1, ell=2, Ngrid=960)

    # Running March 13, 12:06pm on Harmattan
    build_multipro('pk', 'nseries', 'noweight', range(11,21), Nthreads=1, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'noweight', range(11,21), Nthreads=1, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'nseries', 'noweight', range(21,31), Nthreads=1, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'noweight', range(21,31), Nthreads=1, clobber=True, ell=2, Ngrid=960)

    # Running March 13, 12:04pm on Masa
    build_multipro('pk', 'nseries', 'true', [9], Nthreads=1, ell=2, Ngrid=960)
    build_multipro('pk', 'nseries', 'upweight', range(31,51), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'nseries', 'dlospeak', range(31,51), Nthreads=4, clobber=True, ell=2, Ngrid=960)

    build_multipro('pk', 'qpm', 'true', range(21,51), Nthreads=4, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'upweight', range(11,51), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'dlospeak', range(11,51), Nthreads=4, clobber=True, ell=2, Ngrid=960)

    # Finish March 15, 6:25pm on Masa 
    #build_multipro('pk', 'nseries', 'true', [10, 12, 20, 84], Nthreads=4, ell=2, Ngrid=960)

    # Running on March 15, 10:03pm on Masa
    build_multipro('pk', 'nseries', 'upweight', range(51,61), Nthreads=5, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'nseries', 'dlospeak', range(51,61), Nthreads=5, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'true', range(51,61), Nthreads=4, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'upweight', range(51,61), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'dlospeak', range(51,61), Nthreads=4, clobber=True, ell=2, Ngrid=960)

    build_multipro('pk', 'nseries', 'upweight', range(61,71), Nthreads=5, clobber=True, ell=2, Ngrid=960)

    # Running on March 16 10:50pm on Masa
    build_multipro('pk', 'nseries', 'dlospeak', range(66,71), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'true', range(61,73), Nthreads=4, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'upweight', range(61,73), Nthreads=4, clobber=True, ell=2, Ngrid=960)

    build_multipro('pk', 'nseries', 'upweight', range(71,85), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'nseries', 'dlospeak', range(71,85), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'true', range(73,81), Nthreads=4, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'upweight', range(73,81), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'dlospeak', range(73,81), Nthreads=4, clobber=True, ell=2, Ngrid=960)

    build_multipro('pk', 'qpm', 'true', range(81,101), Nthreads=4, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'upweight', range(81,101), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'dlospeak', range(81,101), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    '''

    #build_multipro('pk', 'qpm', 'dlospeak', range(61,73), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    #build_multipro('pk', 'nseries', 'noweight', range(51,63), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    #build_multipro('pk', 'qpm', 'noweight', range(51,63), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    #build_multipro('pk', 'nseries', 'noweight', range(63,75), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    #build_multipro('pk', 'qpm', 'noweight', range(63,75), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    #build_multipro('pk', 'nseries', 'noweight', range(75,85), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'noweight', range(75,101), Nthreads=4, clobber=True, ell=2, Ngrid=960)
    
    '''
    # Finished March 15, 10:04pm on Harmattan
    build_multipro('pk', 'nseries', 'noweight', range(31,41), Nthreads=1, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'noweight', range(31,41), Nthreads=1, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'nseries', 'noweight', range(41,51), Nthreads=1, clobber=True, ell=2, Ngrid=960)
    build_multipro('pk', 'qpm', 'noweight', range(41,51), Nthreads=1, clobber=True, ell=2, Ngrid=960)
    '''
