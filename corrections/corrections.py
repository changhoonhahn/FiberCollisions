'''

Corrections class 

'''
import cosmolopy as cosmos

from Spectrum.util.catalog import Catalog

class Corrections(object): 

    def __init__(self, cat_corr, **kwargs): 
        """ Class describing the correction to galaxy catalog of simulations or BOSS data
        """
        self.cat_corr = cat_corr
        self.kwargs = kwargs
        self.corr_str = None 
        self.catclass = Catalog(self.cat_corr)

    def build(self): 
        pass 

    def corrstr(self): 
        pass 
    
    def file(self):
        """ Name of corrected galaxy catalog file 
        """

        cat = Catalog(self.cat_corr)
        file_list = cat.file()      # list of file name elements
   
        # cosmology string 
        try: 
            if self.kwargs['cosmology'] == 'survey': 
                cosmos_str = '.sureycosmo'
            else: 
                cosmos_str = '.fidcosmo'
        except KeyError: 
            cosmos_str = '.fidcosmo'
        file_list.insert( -1, cosmos_str )

        # correction string 
        if self.corr_str != None: 
            corr_str = self.corr_str

            file_list.insert( -1, corr_str )

            return ''.join(file_list)
        else: 
            return None

    def cosmo(self): 
        """ Default cosmology used for the correction 
        """
    
        try: 
            if self.kwargs['cosmology'] == 'survey': 
                # survey cosmology
                cat = Catalog(self.cat_corr)
                self.cosmo = cat.cosmo()

                return self.cosmo
            else: 
                # default fiducial cosmology (hardcoded)
                omega_m = 0.31 

        except KeyError: 
            omega_m = 0.31  # default 

        # survey cosmology 
        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.676
        cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        self.cosmo = cosmo 

        return self.cosmo 

    def survey_zlimits(self): 
        """ Catalog survey limits
        """
        return (self.catclass).survey_zlimits()

    def datacolumns(self): 
        """ Data columns for given catalog and correction
        """
        return (self.catclass).datacolumns()

    def datacols_fmt(self): 
        """ Data format of columns of catalog data
        """
        return (self.catclass).datacols_fmt()

    def datacols_header(self): 
        """ Header string that describes data columsn
        """
        return (self.catclass).datacols_header()
