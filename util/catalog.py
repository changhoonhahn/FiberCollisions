'''



'''

class Catalog(object): 

    def __init__(self, cat_corr): 
        """ Class describing simulated/data catalog 
        """
        self.cat_corr = cat_corr.copy()
        self.catalog_name = (cat_corr['catalog'])['name'].lower()

        if self.catalog_name in ('nseries'): 
            self.cat_col_dict_key = self.catalog_name
        elif 'cmass' in self.catalog_name: 
            self.cat_col_dict_key = 'cmass' 
        else: 
            raise NotImplementedError()
    
        # dictionary of catlaog columsn 
        self.cat_col_dict = {
                'nseries': {
                    'cols': ['ra', 'dec', 'z', 'wfc', 'comp'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                    'hdrs': "Column : ra, dec, z, wfc, comp"
                    }, 
                'cmass': {
                    'cols': ['ra', 'dec', 'z', 'nbar', 'wsys', 'wnoz', 'wfc', 'comp'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                    'hdrs': 'Columns : ra, dec, z, nbar, w_systot, w_noz, w_cp, comp'
                    }
                }

    def file(self): 
        """ File elements that pertain to the simulated/data catalog
        """
        cat = self.cat_corr['catalog']
        
        if self.catalog_name == 'nseries':         # Nseries

            data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
            file_beg = ''.join(['CutskyN', str(cat['n_mock'])])
            file_end = '.dat'

        elif 'cmass' in self.catalog_name: 

            if self.catalog_name == 'cmass': 
                data_dir = '/mount/riachuelo1/hahn/data/CMASS/'

            elif 'cmasslowz' in self.catalog_name:
                data_dir = '/mount/riachuelo1/hahn/data/CMASS/dr12v5/'
                file_beg = ''.join(['galaxy_DR12v5_', self.catalog_name.upper(), '_North'])
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

    def survey_zlimits(self): 
        """ Redshift limits of survey 
        """
        
        # survey redshift limits 
        if self.catalog_name in ('lasdamasgeo', 'ldgdownnz'):     
            survey_zmin, survey_zmax = 0.16, 0.44

        elif self.catalog_name in ('tilingmock', 'qpm', 'patchy', 'nseries'): 
            survey_zmin, survey_zmax = 0.43, 0.7 

        elif 'bigmd' in self.catalog_name:             
            survey_zmin, survey_zmax = 0.43, 0.7    

        elif 'cmass' in self.catalog_name:             

            if self.catalog_name == 'cmass': 
                survey_zmin, survey_zmax = 0.43, 0.7  

            elif 'cmasslowz' in self.catalog_name: 

                if '_high' in self.catalog_name: 
                    survey_zmin, survey_zmax = 0.5, 0.75    

                elif '_low' in self.catalog_name: 
                    survey_zmin, survey_zmax = 0.2, 0.5

            else: 
                raise NotImplementedError('CMASS or CMASSLOWZ combined sample')

        else: 
            raise NotImplementedError('Mock Catalog not included')

        return [survey_zmin, survey_zmax]

    def datacolumns(self): 
        """ Columns of catalog data
        """
        return (self.cat_col_dict[self.cat_col_dict_key])['cols']

    def datacols_fmt(self): 
        """ Data format of columns of catalog data
        """
        return (self.cat_col_dict[self.cat_col_dict_key])['fmts']

    def datacols_header(self): 
        """ Header string that describes data columsn
        """
        return (self.cat_col_dict[self.cat_col_dict_key])['hdrs']
