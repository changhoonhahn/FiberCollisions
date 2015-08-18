import multiprocessing as mp 
import numpy as np
import fibcol_data as fc_data


def f(x): 
    print x**3
    return x**3

if __name__=="__main__": 

    pool = mp.Pool(processes=5) 
    mapfn = pool.map

    arglist = np.arange(1,101)
    arglists = [arglist]
    #arglists = np.split(arglist, 5) 

    result = mapfn(f, [ags for ags in arglists])

    pool.close()
    pool.terminate()
    pool.join() 
