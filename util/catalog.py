'''



'''

class catalog: 

    def __init__(self, cat_corr): 
        """ Class describing simulated/data catalog 
        """
        self.cat_corr = cat_corr
        self.catalog_name = (cat_corr['catalog'])['name']

        if self.catalog_name not in ('nseries'): 
            raise NotImplementedError()

        pass

    def file(self): 
        """ File elements that pertain to the simulated/data catalog
        """
        cat = self.cat_corr['catalog']
        
        if self.catalog_name == 'nseries':         # Nseries

            data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
            file_beg = ''.join(['CutskyN', str(cat['n_mock'])])
            file_end = '.dat'
    
        return [data_dir, file_beg, file_end]
