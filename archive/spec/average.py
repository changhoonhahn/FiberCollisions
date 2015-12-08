'''

Average spectrum

'''
import numpy as np
import os.path
import subprocess
import cosmolopy as cosmos

from spec import Spec

class AvgSpec(Spec):

    def __init__(self, n_mocks, spectype, cat_corr, ell=None, **kwargs): 
        """ Child class of Spec object class in spec.py that describes 
        the average power/bispectrum 

        """
        if isinstance(n_mocks, list) or isinstance(n_mocks, np.ndarray): 
            self.n_mocks = len(n_mocks)
            self.n_mocks_list = n_mocks
        else: 
            self.n_mocks = n_mocks
            self.n_mocks_list = range(1, n_mocks+1)
        
        super(AvgSpec, self).__init__(spectype, cat_corr, ell=ell, **kwargs)
        
        self.ell = ell 
        self.k = None
        self.p0k = None
        self.p2k = None
        self.p4k = None
        self.avg_spec = None
        self.spec_dev = None

    def variance(self): 
        '''
        Calculate variance of power/bispectrum
        '''
        # calculate average if not already calculated
        if self.k is None: 
            if os.path.isfile(self.file_name): 
                self.read()
            else: 
                self.build()
    
        for i_mock in self.n_mocks_list: 
            k_i, spec_i_spec = self.spec_i(i_mock)
            
            try: 
                var += (self.avg_spec - spec_i_spec)**2
            except UnboundLocalError: 
                var = (self.avg_spec - spec_i_spec)**2
        var /= np.float(self.n_mocks) 

        self.spec_var = var
        return var

    def stddev(self): 
        '''
        Calculate standard deviation 
        '''
        if self.spec_dev is None: 
            self.variance()
        
        return np.sqrt(self.spec_var)
    
    def spec_i(self, i_mock):
        '''
        k, power/bispectrum (p0k, p2k, p4k, bk) of ith mock
        '''
        specdict = self.cat_corr['spec']

        # copy cat_corr dictionary 
        cat_corr_i = self.cat_corr.copy() 
        cat_corr_i['catalog']['n_mock'] = i_mock

        spec_i = Spec(self.type, cat_corr_i, **self.kwargs)
        spec_i.read() 

        if self.type == 'pk':
            spec_ell = ''.join(['p', str(specdict['ell']), 'k'])
        else: 
            raise NotImplementedError
        
        spec_i_spec = getattr(spec_i, spec_ell)

        return [spec_i.k, spec_i_spec]

    def file(self): 
        '''
        File name of average power/bispectrum
        '''
        specdict = self.cat_corr['spec']
        self.cat_corr['catalog']['n_mock']  = 1

        spec_file = super(AvgSpec, self).file()
        spec_file_name = spec_file.split('/')[-1]
        spec_file_core = spec_file_name.split('.')[0]
        spec_file_ending = '.'.join(spec_file_name.split('.')[1:])

        if self.cat_corr['catalog']['name'] == 'nseries': 
            avg_file = ''.join([
                '/'.join(spec_file.split('/')[:-1]), '/', 
                'AVG_P', str(specdict['ell']), 'K_', 
                spec_file_core.split('1')[0], '.', 
                str(self.n_mocks), 'mocks.', 
                spec_file_ending
                ])
        return avg_file
    
    def build(self):
        '''
        Calculate average spectrum
        '''
        specdict = self.cat_corr['spec']
    
        for i_mock in self.n_mocks_list: 
            
            spec_i_k, spec_i_spec = self.spec_i(i_mock)
            
            try: 
                spec_sum += spec_i_spec
            except UnboundLocalError:
                k = spec_i_k
                spec_sum = spec_i_spec
        
        self.k = k
        self.avg_spec = spec_sum/np.float(self.n_mocks)
        
        if self.type == 'pk':
            spec_ell = ''.join(['p', str(specdict['ell']), 'k'])
            setattr(self, spec_ell, self.avg_spec)
        else: 
            raise NotImplementedError

        return [self.k, self.avg_spec] 

    def writeout(self): 
        '''
        Write average power/bispectrum to file 
        '''
        specdict = self.cat_corr['spec']

        if self.k is None: 
            self.build()
    
        spec_ell_str = ''.join(['p', str(specdict['ell']), 'k'])
        data_list = [self.k, getattr(self, spec_ell_str)]  
        data_hdrs = ''.join(['# Columns k, ', spec_ell_str])

        output_file = self.file_name
        np.savetxt(
                output_file, 
                (np.vstack(np.array(data_list))).T, 
                fmt=['%10.5f', '%10.5f'], 
                delimiter='\t', 
                header=data_hdrs
                ) 

        return None

    def read(self): 
        '''
        Read average power/bispectrum from file 
        '''
        specdict = self.cat_corr['spec']
        
        if self.k is None: 
            self.build()
            self.writeout()

        k, avg_spec = np.loadtxt(self.file_name, skiprows=1, unpack=True, usecols=[0,1])

        self.k = k 
        spec_ell_str = ''.join(['p', str(specdict['ell']), 'k'])
        setattr(self, spec_ell_str, avg_spec)

        return None
        
if __name__=="__main__":
    cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 10}, 
            'correction': {'name': 'true'}
            }
    blah = AvgSpec(10, 'pk', cat_corr, ell=4, Ngrid=960)
    print blah.file()
    #blah.read()
    blah.build()
    #print blah.variance()
    #print blah.p4k[-10:-1]
