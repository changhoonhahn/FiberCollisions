'''

Corrections class 

'''
from util import catalog as cata 

class Corrections(object): 

    def __init__(self, cat_corr, **kwargs): 
        """ Class describing the correction to galaxy catalog of simulations or BOSS data
        """
        self.cat_corr = cat_corr
        self.kwargs = kwargs
        self.corr_str = None 

    def build(self): 
        pass 

    def corrstr(self): 
        pass 
    
    def file(self):
        """ Name of corrected galaxy catalog file 
        """
        cat = cata.catalog(self.cat_corr)
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
