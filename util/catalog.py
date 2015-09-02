'''



'''

class catalog: 

    def __init__(self, cat_corr): 
        """ Class describing simulated/data catalog 
        """
        self.cat_corr = cat_corr
        self.catalog_name = (cat_corr['catalog'])['name'].lower()

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

    def cosmo(self): 
        """ Survey cosmology. Note this is different than the fiducial cosmology
        used to compute powerspectrum 
        """

        if self.catalog_name == 'nseries': 
            omega_m =  0.286
        else: 
            raise NotImplementedError()
        
        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.676
        cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        self.cosmo = cosmo 

        return cosmo 

