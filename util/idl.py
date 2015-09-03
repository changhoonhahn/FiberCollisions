'''

Interface with IDL code 


'''

import os 

class Idl(object): 

    def __init__(self, purpose, **kwargs): 
        """ Class that describes the IDL code 
        """
        
        self.kwargs = kwargs
        self.idlrun = 'idl -e "'
        self.idl_dir = '/home/users/hahn/powercode/FiberCollisions/pro/'
        if purpose == 'dlos': 
            self.cmd = self.build_dlos()
        else: 
            raise NotImplementedError()

    def build_dlos(self):
        """ Build line-of-sight displacement 
        """

        if not all([x in (self.kwargs).keys() for x in ['catalog_name', 'galaxy_file', 'dlos_file']]): 
            err_msg = "Keywords 'catalog_name', 'galaxy_file', 'dlos_file' must be specified"
            raise ValueError(err_msg) 

        idl_cmd = ''.join([
            self.idlrun, 
            'build_fibcoll_dlos_cosmo, ', 
            "'", (self.kwargs)['catalog_name'], "', ", 
            "'", (self.kwargs)['galaxy_file'], "', ", 
            "'", (self.kwargs)['dlos_file'], "'", 
            '"'
            ])

        return idl_cmd 

    def run(self): 
        """ Run IDL command 
        """
        print self.cmd 
        os.system(self.cmd)

        return None
