'''

Make dLOS distribution peak corrected CMASS and LOWZ 
combined galaxy catalogs 

'''
from data import Data

if __name__=="__main__": 
    for zbin_str in ['_low', '_high']:
        for specifier_str in ['', 'e2', 'e3']: 
            cat_corr = {
                    'catalog': {'name': 'cmasslowz'+specifier_str+zbin_str}, 
                    'correction': {'name': 'upweight'} 
                    }
            corrdata = Data('data', cat_corr)
            print corrdata.file_name
            print corrdata.build()
