import multiprocessing as mp 
import numpy as np
import fibcol_data as fc_data


def f(DorR, cat_corr): 
    print DorR
    print cat_corr
    #fc_data.galaxy_data('data', **cat_corr) 
    return 

if __name__=="__main__": 
    pool = mp.Pool(processes=5) 
    mapfn = pool.map

    arglist = [
            ('data', 
                {'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
                    'correction': {'name': 'photozpeakshot', 'fit': 'gauss', 'sigma': 4.0, 'fpeak': 0.69}
                    }
                ) 
            for i_mock in np.arange(1,101)]

    result = mapfn(f, [ags for ags in arglist])

    pool.close()
    pool.terminate()
    pool.join() 
