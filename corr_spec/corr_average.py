'''

Corrected Average spectrum

'''
import numpy as np
import os.path
import subprocess
import cosmolopy as cosmos

from corr_spec import CorrSpec

class CorrAvgSpec(CorrSpec):

    def __init__(self, n_mocks, spectype, cat_corr, ell=None, **kwargs): 
        ''' 
        Child class of Spec object class in spec.py that describes 
        the average power/bispectrum 

        ''' 
        if isinstance(n_mocks, list) or isinstance(n_mocks, np.ndarray): 
            self.n_mocks = len(n_mocks)
            self.n_mocks_list = n_mocks
        else: 
            self.n_mocks = n_mocks
            self.n_mocks_list = range(1, n_mocks+1)
        
        super(CorrAvgSpec, self).__init__(spectype, cat_corr, ell=ell, **kwargs)
        
        self.ell = ell 
        self.rebin = None
        if self.type == 'pk': 
            self.k = None
            self.p0k = None
            self.p2k = None
            self.p4k = None

            self.avg_spec = None
            self.spec_var = None

        elif self.type == 'bk': 
            self.k1 = None
            self.k2 = None
            self.k3 = None
            self.p0k1 = None
            self.p0k2 = None
            self.p0k3 = None
            self.bk = None
            self.qk = None

            self.avg_bk = None
            self.bk_var = None
            
            self.avg_qk = None
            self.qk_var = None

    def covariance(self, krange=None, rebin=None): 
        ''' Calculate the covariance matrix of the power/bispectrum 
        '''
        if self.type == 'pk': 
            # calculate average if not already calculated
            specs = [] 
            for i_mock in self.n_mocks_list: 
                if krange is not None: 
                    k_i, spec_i_spec, count_i = self.spec_i(i_mock, rebin=rebin)
                        
                    klim = np.where(
                            (k_i > krange[0]) & 
                            (k_i <= krange[1])
                            )
                    specs.append(spec_i_spec[klim])
                else: 
                    k_i, spec_i_spec, count_i = self.spec_i(i_mock, rebin=rebin)
                    specs.append(spec_i_spec)

            #print np.vstack(specs)
            covmat = np.cov(np.array(np.vstack(specs)).T)
        else: 
            raise ValueError
        return covmat

    def variance(self, rebin=None): 
        ''' Calculate variance of power/bispectrum
        '''
        if self.type == 'pk': 
            # calculate average if not already calculated
            if self.k is None: 
                if os.path.isfile(self.file_name): 
                    self.read()
                else: 
                    self.build()
        
            for i_mock in self.n_mocks_list: 
                k_i, spec_i_spec, count = self.spec_i(i_mock, rebin=rebin)
                
                try: 
                    var += (self.avg_spec - spec_i_spec)**2
                except UnboundLocalError: 
                    var = (self.avg_spec - spec_i_spec)**2
            var /= np.float(self.n_mocks-1) 

            self.spec_var = var
        
            return var

        elif self.type == 'bk': 

            if self.k1 is None: 
                if os.path.isfile(self.file_name): 
                    self.read()
                else: 
                    self.build()

            for i_mock in self.n_mocks_list: 
                k1_i, k2_i, k3_i, spec_i_bk, spec_i_qk = self.spec_i(i_mock)
                
                try: 
                    bk_var += (self.avg_bk - spec_i_bk)**2
                    qk_var += (self.avg_qk - spec_i_qk)**2
                except UnboundLocalError: 
                    bk_var = (self.avg_bk - spec_i_bk)**2
                    qk_var = (self.avg_qk - spec_i_qk)**2

            bk_var /= np.float(self.n_mocks) 
            qk_var /= np.float(self.n_mocks) 

            self.bk_var = bk_var
            self.qk_var = qk_var

            return bk_var, qk_var

    def stddev(self, rebin=None): 
        '''
        Calculate standard deviation 
        '''
        if self.type == 'pk':
            if self.spec_var is None:
                self.variance(rebin=rebin)

            return np.sqrt(self.spec_var)

        elif self.type == 'bk': 

            if self.bk_var is None: 
                self.variance()
            
            return np.sqrt(self.bk_var), npsqrt(self.qk_var)
    
    def spec_i(self, i_mock, rebin=None):
        '''
        k, power/bispectrum (p0k, p2k, p4k, bk) of ith mock
        '''
        specdict = self.cat_corr['spec']

        # copy cat_corr dictionary 
        cat_corr_i = self.cat_corr.copy() 
        cat_corr_i['catalog']['n_mock'] = i_mock
        spec_i = CorrSpec(self.type, cat_corr_i, **self.kwargs)
        spec_i.read() 

        if self.type == 'pk':
            spec_ell = ''.join(['p', str(specdict['ell']), 'k'])

            spec_i_spec = getattr(spec_i, spec_ell)

            return self._SpecRebin([spec_i.k, spec_i_spec, spec_i.count], rebin=rebin)

        elif self.type == 'bk': 
            bk_i = getattr(spec_i, 'bk')
            qk_i = getattr(spec_i, 'qk')

            return self._SpecRebin([spec_i.k1, spec_i.k2, spec_i.k3, bk_i, qk_i, spec_i.count], rebin=rebin)
        
        else: 
            raise NotImplementedError

    def _SpecRebin(self, input_list, rebin=None): 
        ''' Assumes that input_list is structured so that input_list[0] is k while input_list[-1] is the counts
        '''
        if rebin is None: 
            return input_list 

        k = input_list[0]
        counts = input_list[-1]
        specs = input_list[1:-1]
        
        tot_counts = [] 
        avg_k = [] 
        avg_specs = [[] for spec in specs]
        if isinstance(rebin, int): 

            i_threshold = (np.abs(k[::rebin] - 0.05)).argmin()
            for ii in range(i_threshold*rebin)[::2]: 
                indices = range(len(k))[ii:ii+2]

                # all mode counts 
                allcounts = np.sum(counts[indices])
                tot_counts.append(allcounts)
                avg_k.append(np.sum(k[indices] * counts[indices])/allcounts)

                for i_spec, spec in enumerate(specs): 
                    avg_specs[i_spec].append(np.sum(spec[indices] * counts[indices])/allcounts)
                #tot_counts.append(counts[ii])
                #avg_k.append(k[ii])
                #for i_spec, spec in enumerate(specs): 
                #    avg_specs[i_spec].append(spec[ii])

            for istart in range(len(k))[::rebin]: 

                if istart < i_threshold*rebin: 
                    continue 
                indices = range(len(k))[istart:istart+rebin]

                # all mode counts 
                allcounts = np.sum(counts[indices])
                tot_counts.append(allcounts)
                avg_k.append(np.sum(k[indices] * counts[indices])/allcounts)

                for i_spec, spec in enumerate(specs): 
                    avg_specs[i_spec].append(np.sum(spec[indices] * counts[indices])/allcounts)

            avg_specs = [np.array(spec) for spec in avg_specs]
            avg_specs.append(np.array(tot_counts))
            avg_specs.reverse()
            avg_specs.append(np.array(avg_k))
            avg_specs.reverse()

            return avg_specs

        if isinstance(rebin, str): 
            if rebin == 'log': 
                bins = np.logspace(-3, 0, 50)
                klow = bins[:-1]
                khigh = bins[1:]

                for ik in range(len(klow)): 
                    kbin = np.where((k > klow[ik]) & (k <= khigh[ik]))
                    if len(kbin[0]) == 0: 
                        continue 

                    # all mode counts 
                    allcounts = np.sum(counts[kbin])
                    tot_counts.append(allcounts)
                    avg_k.append(np.sum(k[kbin] * counts[kbin])/allcounts)

                    for i_spec, spec in enumerate(specs): 
                        avg_specs[i_spec].append(np.sum(spec[kbin] * counts[kbin])/allcounts)

                avg_specs = [np.array(spec) for spec in avg_specs]
                avg_specs.append(np.array(tot_counts))
                avg_specs.reverse()
                avg_specs.append(np.array(avg_k))
                avg_specs.reverse()

                return avg_specs

    def file(self, rebin=None): 
        '''
        File name of average power/bispectrum
        '''
        specdict = self.cat_corr['spec']
        self.cat_corr['catalog']['n_mock']  = 1

        spec_file = super(CorrAvgSpec, self).file()
        spec_file_name = spec_file.split('/')[-1]
        spec_file_core = spec_file_name.split('.')[0]
        spec_file_ending = '.'.join(spec_file_name.split('.')[1:])
    
        rebin_str = ''
        if rebin is not None: 
            if isinstance(rebin, int): 
                rebin_str = '.rebin'+str(rebin)+'x'
            elif isinstance(rebin, str): 
                rebin_str = '.rebin'+rebin
            else: 
                raise ValueError

        #if self.cat_corr['catalog']['name'] == 'nseries': 
        if self.type == 'pk': 
            if self.n_mocks != 1: 
                avg_file = ''.join([
                    '/'.join(spec_file.split('/')[:-1]), '/', 
                    'AVG_P', str(specdict['ell']), 'K_', 
                    spec_file_core.split('1')[0], '.', 
                    str(self.n_mocks), 'mocks', 
                    rebin_str, '.', 
                    spec_file_ending
                    ])
            else: 
                avg_file = ''.join([
                    '/'.join(spec_file.split('/')[:-1]), '/', 
                    'AVG_P', str(specdict['ell']), 'K_', 
                    spec_file_core, '.', 
                    str(self.n_mocks), 'mocks', 
                    rebin_str, '.', 
                    spec_file_ending
                    ])

        elif self.type == 'bk': 
            avg_file = ''.join([
                '/'.join(spec_file.split('/')[:-1]), '/', 
                'AVG_BK_', 
                spec_file_core.split('1')[0], '.', 
                str(self.n_mocks), 'mocks.', 
                rebin_str, 
                spec_file_ending
                ])

        return avg_file
    
    def build(self, rebin=None):
        '''
        Calculate average spectrum
        '''
        specdict = self.cat_corr['spec']
    
        if self.type == 'pk':       # powerspectrum

            for i_mock in self.n_mocks_list: 
                spec_i_k, spec_i_spec, count_i = self.spec_i(i_mock, rebin=rebin)
                
                try: 
                    count_sum += count_i
                    k_sum += spec_i_k * count_i
                    spec_sum += spec_i_spec * count_i
                except UnboundLocalError:
                    count_sum = count_i
                    k_sum = spec_i_k * count_i
                    spec_sum = spec_i_spec * count_i
            
            self.k = k_sum / count_sum
            self.avg_spec = spec_sum/count_sum
            self.count = count_sum 
        
            spec_ell = ''.join(['p', str(specdict['ell']), 'k'])
            setattr(self, spec_ell, self.avg_spec)

            return [self.k, self.avg_spec] 

        elif self.type == 'bk':     # bispectrum
            for i_mock in self.n_mocks_list: 
                spec_i_k1, spec_i_k2, spec_i_k3, spec_i_bk, spec_i_qk = self.spec_i(i_mock, rebin=rebin)

                try: 
                    spec_bk_sum += spec_i_bk
                    spec_qk_sum += spec_i_qk
                except UnboundLocalError:
                    k1 = spec_i_k1
                    k2 = spec_i_k2
                    k3 = spec_i_k3
                    spec_bk_sum = spec_i_bk
                    spec_qk_sum = spec_i_qk
            
            self.k1 = k1
            self.k2 = k2
            self.k3 = k3
            self.avg_bk = spec_bk_sum/np.float(self.n_mocks)
            self.avg_qk = spec_qk_sum/np.float(self.n_mocks)
        
            return [self.k1, self.k2, self.k3, self.avg_bk, self.avg_qk] 

        else: 
            raise NotImplementedError

    def _build_kbin(self, kbin):
        '''
        Calculate average spectrum
        '''
        specdict = self.cat_corr['spec']
    
        if self.type == 'pk':       # powerspectrum

            for i_mock in self.n_mocks_list: 
                spec_i_k, spec_i_spec = self.spec_i(i_mock)
                
                try: 
                    spec_sum += spec_i_spec
                except UnboundLocalError:
                    k = spec_i_k
                    spec_sum = spec_i_spec
             
            avg_spec = spec_sum/np.float(self.n_mocks)

            # rebin 
            if kbin == 'double': 
                if len(k[0::2]) == len(k[1::2]): 
                    rebin_k = 0.5 * (k[0::2] + k[1::2])
                    rebin_avg_spec = 0.5 * (avg_spec[0::2] + avg_spec[1::2])
                else: 
                    rebin_k = 0.5 * (k[2::2] + k[1::2])
                    rebin_avg_spec = 0.5 * (avg_spec[2::2] + avg_spec[1::2])
            else: 
                raise NotImplementedError
            
            self.k = rebin_k
            self.avg_spec = rebin_avg_spec
        
            spec_ell = ''.join(['p', str(specdict['ell']), 'k'])
            setattr(self, spec_ell, self.avg_spec)

            return [self.k, self.avg_spec] 

        else: 
            raise NotImplementedError('Only implemented for power spectrum')

    def _variance_kbin(self, kbin):
        '''
        Calculate variance of the power spectrum with different kbin
        '''
        self._build_kbin(kbin)

        if self.type == 'pk': 
            # calculate average if not already calculated
            if self.k is None: 
                if os.path.isfile(self.file_name): 
                    self.read()
                else: 
                    self.build()
        
            for i_mock in self.n_mocks_list: 
                k_i, spec_i_spec = self.spec_i(i_mock)

                if kbin == 'double': 
                    if len(k_i[0::2]) == len(k_i[1::2]): 
                        rebin_k = 0.5 * (k_i[0::2] + k_i[1::2])
                        rebin_spec_i = 0.5 * (spec_i_spec[0::2] + spec_i_spec[1::2])
                    else: 
                        rebin_k = 0.5 * (k_i[2::2] + k_i[1::2])
                        rebin_spec_i = 0.5 * (spec_i_spec[2::2] + spec_i_spec[1::2])
                else: 
                    raise NotImplementedError
                    
                try: 
                    var += (self.avg_spec - rebin_spec_i)**2
                except UnboundLocalError: 
                    var = (self.avg_spec - rebin_spec_i)**2

            var /= np.float(self.n_mocks) 

            self.spec_var = var
        
            return var
        else: 
            raise NotImplementedError

    def _stddev_kbin(self, kbin): 
        '''
        Calculate standard deviation 
        '''
        if self.type == 'pk':
            spec_var = self._variance_kbin(kbin)

            return np.sqrt(spec_var)
        else: 
            raise NotImplementedError

    def writeout(self): 
        '''
        Write average power/bispectrum to file 
        '''
        specdict = self.cat_corr['spec']
    
        if self.type == 'pk': 
            if self.k is None: 
                self.build()
    
            spec_ell_str = ''.join(['p', str(specdict['ell']), 'k'])
            data_list = [self.k, getattr(self, spec_ell_str), self.count]  
            data_hdrs = ''.join(['# Columns k, ', spec_ell_str, ', count'])
            data_fmts = ['%10.5f', '%10.5f', '%10.5f']

        elif self.type == 'bk': 
            if self.k1 is None: 
                self.build()
            
            data_list = [
                    self.k1, 
                    self.k2, 
                    self.k3, 
                    getattr(self, 'avg_bk'),
                    getattr(self, 'avg_qk')
                    ]  
            data_hdrs = ''.join(['# Columns k1, k2, k3, Bk, Qk'])
            data_fmts = ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f']

        output_file = self.file_name
        np.savetxt(
                output_file, 
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmts, 
                delimiter='\t', 
                header=data_hdrs
                ) 

        return None

    def read(self, rebin=None): 
        '''
        Read average power/bispectrum from file 
        '''
        self.file_name = self.file(rebin=rebin)

        specdict = self.cat_corr['spec']
        
        if self.type == 'pk': 
            if self.k is None: 
                self.build(rebin=rebin)
                self.writeout()

            k, avg_spec, count = np.loadtxt(self.file_name, skiprows=1, unpack=True, usecols=[0,1,2])

            self.k = k 
            spec_ell_str = ''.join(['p', str(specdict['ell']), 'k'])
            setattr(self, spec_ell_str, avg_spec)
            self.count = count

        elif self.type == 'bk': 
            if self.k1 is None: 
                self.build()
                self.writeout()

            k1, k2, k3, avg_bk, avg_qk = np.loadtxt(self.file_name, skiprows=1, unpack=True, usecols=[0,1,2,3,4])

            self.k1 = k1 
            self.k2 = k2 
            self.k3 = k3 
            self.avg_bk = avg_bk
            self.avg_qk = avg_qk

        return None
