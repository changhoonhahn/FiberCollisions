
# --- Universal --- 
import os
import pickle
import numpy as np
from scipy.interpolate import interp1d

# Spectrum 
from Spectrum.spec import Spec
from Spectrum.util.direc import direc

from corr_fft import CorrFft
from corr_corrdata import CorrCorrData

from fourier_corr import fourier_corr as fourcorr
from fourier_corr.pk_corr import fourier_tophat_Pk
#from pk_corr import fourier_tophat_Pk


class CorrSpec(Spec): 
    def __init__(self, spectype, cat_corr, ell=None, **kwargs): 
        '''
        '''
        try: 
            self.gal_data = CorrCorrData(cat_corr, **kwargs)
        except NotImplementedError:
            pass
        try: 
            self.datafft = CorrFft('data', cat_corr, **kwargs)
        except NotImplementedError:
            pass
        try: 
            self.randfft = CorrFft('random', cat_corr, **kwargs)
        except NotImplementedError:
            pass
        super(CorrSpec, self).__init__(spectype, cat_corr, ell=ell, **kwargs)

    def file(self): 
        '''
        Return power/bispectrum file name strings. Spectrum.spec plus exceptions.
        '''
        if self.cat_corr['correction']['name'].lower() == 'tophat_conv': 
            return self._TophatConv_file()
        elif self.cat_corr['correction']['name'].lower() == 'floriansn': 
            fc_file = super(CorrSpec, self).file()
            return '.floriansn.dat'.join(fc_file.split('.dat'))
        elif self.cat_corr['correction']['name'].lower() == 'hectorsn': 
            fc_file = super(CorrSpec, self).file()
            return '.hectorsn.dat'.join(fc_file.split('.dat'))
        else: 
            return super(CorrSpec, self).file()

    def read(self): 
        '''
        Read bispectrum/powerspectrum file
        '''
        if self.cat_corr['catalog']['name'].lower() == 'nseriesbox': 
            col_index = [0, -1, 2, 3, 4, 5, 6]
            data_cols = ['k', 'p0k', 'p2k', 'p4k', 'p6k', 'p8k', 'p10k']

            spec_data = np.loadtxt(
                        self.file_name, 
                        unpack = True, 
                        usecols = col_index
                        )
            for i_col, col in enumerate(data_cols): 
                if col != 'k': 
                    setattr(self, col, (2.*np.pi)**3 * spec_data[i_col])
                else: 
                    setattr(self, col, spec_data[i_col])
            return None
        elif self.cat_corr['correction']['name'].lower() == 'tophat_conv': 
            return self._TophatConv_read()
        else: 
            return super(CorrSpec, self).read()
    
    def build(self):
        '''
        '''
        return super(CorrSpec, self).build()

    def _TophatConv_file(self): 
        ''' Powerspectrum file for Tophat Convolution corrected power spectrum
        '''
        cat_dict = self.cat_corr['catalog']
        corr_dict = self.cat_corr['correction']
        spec_dict = self.cat_corr['spec']

        # only supported for Nseries
        if cat_dict['name'].lower() != 'nseries': 
            raise ValueError('Tophat Convolution is only supported for Nseries')
        if not all(x in corr_dict.keys() for x in ['fs', 'rc', 'fold', 'rebin']): 
            err_msg = ''.join(['Specify ', 
                'fs (collided fraction), ', 
                'rc (fibcollision comoving radius), \n ', 
                'fold (folding number), ', 
                'rebin (k rebinning) in the correction dictionary'])
            raise KeyError(err_msg)

        # spectrum directory
        spec_dir = direc('spec', self.cat_corr)
        # powerspectrum or bispectrum specifier 
        if self.type in ('pk'): 
            spec_str = 'POWER_'
        else: 
            raise NotImplementedError
        if (self.type == 'pk') and (self.ell != 0): 
            spec_str += 'Q_'
        # galaxy file that specifies correction 
        gal_file = ''.join([
            'CutskyN', str(cat_dict['n_mock']), '.fidcosmo.dat'
            '.', corr_dict['name'].lower(), 
            '.fs', str(round(corr_dict['fs'], 1)), 
            '.rc', str(round(corr_dict['rc'], 2)), 
            '.fold', str(corr_dict['fold']), 
            '.rebin', str(corr_dict['rebin'])
            ])
        # spectrum specifiers
        if self.type == 'pk': 
            specparam_str = ''.join([
                '.grid', str(spec_dict['Ngrid']), 
                '.P0', str(spec_dict['P0']), 
                '.box', str(spec_dict['Lbox'])
                ])
        else: 
            raise NotImplementedError
    
        file_name = ''.join([
            spec_dir, 
            spec_str,
            gal_file, 
            specparam_str
            ])
        return file_name

    def _TophatConv_read(self): 
        ''' Read tophat convolution corrected power spectrum monopole or quadrupole
        '''
        if os.path.isfile(self.file()): 
            k, corr_pk = np.loadtxt(self.file(), unpack=True, usecols=[0,1])
        else: 
            print 'Building ', self.file()
            self.build()
            k, corr_pk = np.loadtxt(self.file(), unpack=True, usecols=[0,1])
        
        self.k = k 
        setattr(self, 'p'+str(self.ell)+'k', corr_pk)
        return None
    
    def _TophatConv_build(self): 
        ''' Build Tophat Convolution corrected power spectrum monopole and quadrupole.
        '''
        cat_dict = self.cat_corr['catalog']
        corr_dict = self.cat_corr['correction']
        spec_dict = self.cat_corr['spec']
        # read the true power spectrum
        true_cat_corr = {
                'catalog': cat_dict, 
                'correction': {'name': 'true'}, 
                'spec': spec_dict
                }
        true_spec = CorrSpec('pk', true_cat_corr)
        true_spec.read()
        self.k = true_spec.k
        if self.ell == 0:   
            true_pk = true_spec.p0k
        elif self.ell == 2: 
            true_pk = true_spec.p2k
            
        # Del P^uncorr
        k_uncorr, delP_uncorr = self._TophatConv_delP_uncorr()

        # Del P^corr from P(k, mu)
        k_corr, delP_corr = self._TophatConv_delP_corr_pkmu(Ngrid=3600)

        # interpolate DelP^uncorr and DelP^corr to correct true_pk 
        delP_uncorr_interp = interp1d(k_uncorr, delP_uncorr)
        delP_corr_interp = interp1d(k_corr, delP_corr)

        uncorr_range = np.where((self.k >= k_uncorr[0]) & (self.k <= k_uncorr[-1]))
        uncorr_below = np.where(self.k < k_uncorr[0])
        uncorr_above = np.where(self.k > k_uncorr[-1])

        corr_range = np.where((self.k >= k_corr[0]) & (self.k <= k_corr[-1]))
        corr_below = np.where(self.k < k_corr[0])
        corr_above = np.where(self.k > k_corr[-1])
        
        uncorr_delPk = np.zeros(len(self.k))
        uncorr_delPk[uncorr_range] = delP_uncorr_interp(self.k[uncorr_range])
        uncorr_delPk[uncorr_below] = (uncorr_delPk[uncorr_range])[0]
        uncorr_delPk[uncorr_above] = (uncorr_delPk[uncorr_range])[-1]
        
        corr_delPk = np.zeros(len(self.k))
        corr_delPk[corr_range] = delP_corr_interp(self.k[corr_range]) 
        corr_delPk[corr_below] = corr_delPk[corr_range][0]
        corr_delPk[corr_above] = corr_delPk[corr_range][-1]

        # Corrected P_l(k)
        corr_pk = true_pk + uncorr_delPk + corr_delPk

        data_list = [self.k, corr_pk]
        data_fmts = ['%10.5f', '%10.5f']
        np.savetxt(
                self.file(),  
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmts, 
                delimiter='\t'
                ) 
        return None

    def _TophatConv_delP_uncorr(self): 
        '''Uncorrelated term of the tophat convolved Del P. The term is 
        computed as below:

        Del P_uncorr(k)_l = - (2 pi f_s) (pi rc^2) (2 l + 1)/2 Leg_l(0) W_2d(k*rc)/k
        '''
        cat_dict = self.cat_corr['catalog']
        corr_dict = self.cat_corr['correction']
        if cat_dict['name'].lower() == 'nseries':
            uncorr_file = ''.join([ 
                '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                'uncorrdelP', str(self.ell), 'k', 
                '.fs', str(round(corr_dict['fs'], 2)), 
                '.rc', str(round(corr_dict['rc'], 2)),
                '.p'])
        
        # if DelP_uncorr has not been precomputed 
        if not os.path.isfile(uncorr_file):
            # k values from generic Nseries P(k) 
            true_pk_file = ''.join(['/mount/riachuelo1/hahn/power/Nseries/', 
                'POWER_Q_CutskyN1.fidcosmo.dat.grid960.P020000.box3600'])
            k = np.loadtxt(true_pk_file, unpack = True, usecols =[0])

            delp = fourcorr.DelPuncorr(
                    k,                  # k 
                    self.ell,           # ell 
                    fs=corr_dict['fs'], # fs
                    rc=corr_dict['rc']  # rc
                    )
            pickle.dump([k, delp], open(uncorr_file, 'wb'))
        else:     
            k, delp = pickle.load(open(uncorr_file, 'rb')) 
        
        return [k, delp]
    
    def _TophatConv_delP_corr_pkmu(self, Ngrid=3600): 
        ''' Correlated term of the top hat convolution correction. Calculated
        from Nseires folded P(k, mu).  
        '''
        cat_dict = self.cat_corr['catalog']
        corr_dict = self.cat_corr['correction']
        if cat_dict['name'].lower() == 'nseries': 
            corr_file = ''.join([
                '/mount/riachuelo1/hahn/power/Nseries/Box/', 
                'delPcorr_pkmu', 
                '.l', str(self.ell), 
                '.fs', str(round(corr_dict['fs'],2)), 
                '.rc', str(round(corr_dict['rc'],2)), 
                '.Ngrid', str(Ngrid), 
                '.', str(corr_dict['fold']), 'fold', 
                '.', str(corr_dict['rebin']), 'rebin',
                '.p'])

        if not os.path.isfile(corr_file): 
            raise ValueError('Donde esta el file?') 
        else:  
            k_pkmu, delPcorr_pkmu = pickle.load(open(corr_file, 'rb'))
    
        return [k_pkmu, delPcorr_pkmu]

if __name__=='__main__': 
    #cat_corr = { 
    #        'catalog': {'name': 'qso_bigmd', 'version': 'evo'}, 
    #        'correction': {'name': 'true'}, 
    #        'spec':{
    #            'P0': 20000, #P0 
    #            'Lbox': 8000, 
    #            'Ngrid':960, 
    #            'ell': 2 
    #            }
    #        }
    #peek = CorrSpec('pk', cat_corr)
    #peek.build()
    #cat_corr = {'catalog': {'name': 'qso_bigmd', 'version': 'ebossv1.5'}, 
    cat_corr = {'catalog': {'name': 'qso_bigmd', 'version': 'v2'}, 
            'correction': {'name': 'true'}, 
            'spec':{'P0': 6000, 'Lbox': 8000, 'Ngrid':960, 'ell': 0}}
    peek = CorrSpec('pk', cat_corr)
    peek.build()
    '''
    for vers in ['v2-z', 'v2-nsat']: 
        for ell in [0]:#, 2]: 
            cat_corr = { 
                    'catalog': {'name': 'qso_bigmd', 'version': vers}, 
                    'correction': {'name': 'true'}, 
                    'spec':{
                        'P0': 6000, #P0 
                        'Lbox': 8000, 
                        'Ngrid':360, 
                        'ell': ell}}
            peek = CorrSpec('pk', cat_corr)
            peek.build()
    '''
    #for n_jack in range(1, 51): 
    #    cat_corr = { 
    #            'catalog': {'name': 'qso_bigmd', 'version': 'jackknife'+str(n_jack)}, 
    #            'correction': {'name': 'true'}, 
    #            'spec':{
    #                'P0': 20000, #P0 
    #                'Lbox': 8000, 
    #                'Ngrid':960, 
    #                'ell': 0 
    #                }
    #            }
    #    mono = CorrSpec('pk', cat_corr)
    #    mono.build()

    #for i_mock in range(11, 85):
    #    cat_corr = { 
    #            'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
    #            'correction': {'name': 'tophat_conv', 'fs': 0.6, 'rc': 0.43, 'fold': 10, 'rebin': 20}
    #            }
    #    mono = CorrSpec('pk', cat_corr, ell=0, Ngrid=960)
    #    mono.build()

    #    cat_corr = { 
    #            'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
    #            'correction': {'name': 'tophat_conv', 'fs': 0.6, 'rc': 0.43, 'fold': 10, 'rebin': 20}
    #            }
    #    quad = CorrSpec('pk', cat_corr, ell=2, Ngrid=480)
    #    quad.build()
