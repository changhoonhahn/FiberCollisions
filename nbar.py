'''


Code to analyze nbar(z) of mock catalogs


Author(s) : ChangHoon Hahn 


'''

import os.path
import numpy as np
import cosmolopy as cosmos

# --- Local ---
from spec.data import Data
from util.direc import direc

class Nbar: 

    def __init__(self, cat_corr, **kwargs):
        '''
        Class that decribes the nbar(z) of galaxy catalogs

        Parameters
        ----------
        cat_corr : 
            catalog and correction dictionary 

        Notes
        -----
        * nbar(z) is built by first using a nbar-ngal file, then scaling the nbar(z) file 
        * Upweight, Peak+shot correciton methods change nbar(z) by a negligible amount so a correction is unnecssary
        * only coded for lasdamasgeo so far

        '''
        self.cat_corr = cat_corr.copy()
        self.kwargs = kwargs

        # set zlow, zmid, zhigh (consistent with CMASS nbar(z) files) 
        self.zlow = np.arange(0.0, 1.005, 0.005)
        self.zhigh = self.zlow+0.005 
        self.zmid = self.zlow+0.0025 

        self.file_name = self.file()

    def file(self): 
        ''' Get file name of nbar(z) file 
        '''

        gal_data = Data('data', self.cat_corr, **self.kwargs)
        self.data_file = gal_data.file_name
        
        data_dir = direc('data', self.cat_corr)     # data directory 
        gal_file = (gal_data.file_name).split('/')[-1]
    
        nbar_file = ''.join([
            data_dir, 
            'nbar_', 
            gal_file
            ])

        return nbar_file 
    
    def build(self): 
        ''' 
        Calculate nbar(z) for given galaxy catalog 
        '''
        pass

class Ngal(object): 

    def __init__(self, cat_corr, DorR='data', **kwargs): 
        ''' 
        Class that describes Ngal(z_bin)

        Parameters
        ----------
        cat_corr : 
            catalog + correction dictionary 

        Notes
        -----
        * Simple to handle nbar(z) ratios in a cosmology independent way 

        '''
        self.cat_corr = cat_corr.copy()     # save catalog correction dictionary 
        self.kwargs = kwargs
        self.DorR = DorR

        # set zlow, zmid, zhigh (consistent with CMASS nbar(z) files) 
        self.zlow = np.arange(0.0, 1.005, 0.005)
        self.zhigh = self.zlow+0.005 
        self.zmid = self.zlow+0.0025 

        self.ngal = np.zeros(len(self.zlow))

        self.file_name = self.file()

    def file(self): 
        ''' 
        Get file name of Ngal(z) file 
        '''

        gal_data = Data(self.DorR, self.cat_corr, **self.kwargs)
        self.data_file = gal_data.file_name
        
        data_dir = direc(self.DorR, self.cat_corr)     # data directory 

        gal_file = (gal_data.file_name).split('/')[-1]
    
        ngal_file = ''.join([
            data_dir, 
            'ngal_', 
            gal_file
            ])

        return ngal_file 
    
    def build(self): 
        ''' Calculate Ngal(z) from the data file specified by the dictionary   
        '''
        cat_dict = self.cat_corr['catalog']

        gal_data = Data(self.DorR, self.cat_corr, **self.kwargs) 
        gal_data.read()

        redshifts = gal_data.z  # redshifts 
       
        if 'noweights' in self.kwargs.keys():   # do not use the weights
            
            if self.kwargs['noweights']: 
                weights = np.zeros(len(redshifts))
                weights = weights + 1.0 

        else:
            if cat_dict['name'] == 'nseries': 
                weights = gal_data.wfc * (1. / gal_data.comp)
            else: 
                raise NotImplementedError
            """
            if 'cmass' in cat['name'].lower(): 
                # weights and completeness accounted for CMASS like catalogs 
                if DorR == 'data': 
                    weights = data.wsys * (data.wnoz + data.wfc - 1.0) * (1./data.comp) 
                else: 
                    weights = 1./data.comp
            else: 
                raise NotImplementedError("Only CMASS implemented so far") 
            """
    
        ngal = [] 
        for zlow, zhigh in zip(self.zlow, self.zhigh): 

            # for each redshift bin calculate the weighted Ngal values
            zbin = np.where(
                    (redshifts >= zlow) & 
                    (redshifts < zhigh)
                    )
            
            ngal.append( 
                    np.sum(weights[zbin]) 
                    ) 

            #print len(weights[zbin]), np.sum(weights[zbin])

        self.ngal = np.array(ngal)
        data_hdrs = "Columns : zmid, zlow, zhigh, ngal"
        data_list = [self.zmid, self.zlow, self.zhigh, self.ngal]  
        data_fmts = ['%10.5f', '%10.5f', '%10.5f', '%10.5f']

        output_file = self.file()
        np.savetxt(
                output_file, 
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmts, 
                delimiter='\t', 
                header=data_hdrs
                ) 

        return self.ngal 

# Plotting -----
def plot_Ngal(cat_corrs, type='regular', DorR_list=None, **kwargs): 
    ''' Plot Ngal(z) given catalog and correction dictionary

    Parameters
    ----------
    cat_corr : (list of) dictionary specifying catalog and correction  
    type : 'regular', 'ratio' type of plot 

    '''
    # list of catalog and correction methods 
    if not isinstance(cat_corrs, list): 
        cat_corrs = [cat_corrs]
    cat_corr_str = ''
    if DorR_list != None: 
        if not isinstance(DorR_list, list): 
            DorR_list = [DorR_list]

        if len(DorR_list) != len(cat_corrs): 
            raise ValueError('DorR list must have the same number of elements as cat_corr')

    prettyplot()        # set up plot 
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(8,8)) 
    sub = fig.add_subplot(111) 

    for i_cc, cc in enumerate(cat_corrs):    # loop through cat corrs 
        cat = cc['catalog']
        corr = cc['correction']

        if DorR_list[i_cc] == 'random': 
            kwargs['DorR'] = 'random'
        # read Ngal 
        ng = Ngal(cc, **kwargs)
        ng.Read(**kwargs)  

        # catalog correction specifier 
        cat_corr_label = ';'.join([''.join(cat['name'].split('_')).upper(), corr['name'].upper()]) 
        cat_corr_str += ''.join(['_', cat['name'], '_', corr['name']]) 
        
        if DorR_list[i_cc] == 'random': 
            kwargs['DorR'] = 'data'
            # read Ngal 
            D_ng = Ngal(cc, **kwargs)
            D_ng.Read(**kwargs)  
            alpha = np.sum(D_ng.ngal)/np.sum(ng.ngal)
        else: 
            alpha = 1.0

        ng.ngal = alpha * ng.ngal 
        
        if type == 'regular': 
            # plot Ngal(z) 
            sub.plot( ng.zmid, ng.ngal, 
                    c = pretty_colors[i_cc+1], label = cat_corr_label, lw=2) 
        elif type == 'ratio':
            try: 
                sub.plot( ng.zmid, ng.ngal/ngal_denom, 
                        c = pretty_colors[i_cc+1], label = cat_corr_label, lw=2) 
            except NameError: 
                ngal_denom = ng.ngal
                denom_catcorr =  ';'.join([cat['name'].upper(), corr['name'].upper()])    
        
        # redshift limits 
        if cat['name'].lower() == 'lasdamasgeo': 
            cat_zmin = 0.16
            cat_zmax = 0.44
        elif 'cmass' in cat['name'].lower(): 
            if cat['name'].lower() == 'cmass': 
                cat_zmin = 0.43
                cat_zmax = 0.7
            elif cat['name'].lower() == 'cmasslowz_high': 
                cat_zmin = 0.5
                cat_zmax = 0.75
            elif cat['name'].lower() == 'cmasslowz_low': 
                cat_zmin = 0.2
                cat_zmax = 0.5
        else: 
            raise NotImplementedError("Not Implemented Yet")
        
        try: 
            zmin = np.min(zmin, cat_zmin)
            zmax = np.max(zmax, cat_zmax)
        except NameError: 
            zmin = cat_zmin
            zmax = cat_zmax

    if type == 'regular':
        y_label = r"$\mathtt{N_\mathrm{gal}}$ (galaxies)" 
        type_str = ''
    elif type == 'ratio':  
        y_label = ''.join([r"$\mathtt{N_\mathrm{gal}} / N_\mathtt{gal}^{", denom_catcorr, "}$"])
        type_str = '_ratio'

    if 'ylimit' in kwargs.keys(): 
        sub.set_ylim(kwargs['ylimit'])

    sub.set_xlim([zmin, zmax]) 
    sub.set_xlabel('Redshift (z)')
    sub.set_ylabel(y_label)

    sub.legend(loc = 'upper left') 
    
    fig_name = ''.join(['figure/', 
        'NGAL', cat_corr_str, type_str, '.png']) 
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear()
    return None 

if __name__=='__main__':
    
    for i_mock in xrange(1,11): 
        for corr in ['true', 'upweight']:
            cat_corr = {
                    'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
                    'correction': {'name': corr} 
                    }

            nseries_ngal = Ngal(cat_corr)
            nseries_ngal.build()

    """
    cat_corrs = [
            {'catalog': {'name': 'cmasslowz_high'}, 
                'correction': {'name': 'upweight'}}, 
            {'catalog': {'name': 'cmasslowz_high'}, 
                'correction': {'name': 'upweight'}}, 
            {'catalog': {'name': 'cmasslowz_high'}, 
                'correction': {'name': 'peakshot', 'fit': 'gauss', 'sigma': 6.6, 'fpeak': 0.72}}
            ]
    plot_Ngal(cat_corrs, DorR_list=['data', 'random', 'data'])
    plot_Ngal(cat_corrs, type='ratio', DorR_list=['data', 'random', 'data'], ylimit=[0.994, 1.006])
    
    cat_corrs = [
            {'catalog': {'name': 'cmasslowz_low'}, 
                'correction': {'name': 'upweight'}}, 
            {'catalog': {'name': 'cmasslowz_high'}, 
                'correction': {'name': 'upweight'}}, 
            {'catalog': {'name': 'cmasslowz_low'}, 
                'correction': {'name': 'peakshot', 'fit': 'gauss', 'sigma': 6.9, 'fpeak': 0.72}}
            ]
    plot_Ngal(cat_corrs, DorR_list=['data', 'random', 'data'])
    plot_Ngal(cat_corrs, type='ratio', DorR_list=['data', 'random', 'data'], ylimit=[0.994, 1.006])
    """
