'''


Code to analyze nbar(z) of mock catalogs


Author(s) : ChangHoon Hahn 


'''

import numpy as np
import os.path
import subprocess
import cosmolopy as cosmos

# --- Local --- 
import fibcol_data as fc_data
import fibcol_utility as fc_util

# Classes ------------------------------------------------------------
class Nbar: 
    def __init__(self, cat_corr, **kwargs):
        '''Object class to handle nbar(z) data 

        Parameters
        ----------
        cat_corr : catalog and correction dictionary 

        Notes
        -----
        * nbar(z) is built by first using a nbar-ngal file, then scaling the nbar(z) file 
        * Upweight, Peak+shot correciton methods change nbar(z) by a negligible amount so a correction is unnecssary
        * only coded for lasdamasgeo so far

        '''
        self.cat_corr = cat_corr    # save catalog correction dictionary 
        catalog = cat_corr['catalog'] 
        correction  = cat_corr['correction'] 

        # set zlow, zmid, zhigh (consistent with CMASS nbar(z) files) 
        self.zlow = np.arange(0.0, 1.005, 0.005)
        self.zhigh = self.zlow+0.005 
        self.zmid = self.zlow+0.0025 
        n_zbin = len(self.zlow) 

    def File(self, **kwargs): 
        ''' Get file name of nbar(z) file 
        '''
        cat = self.cat_corr['catalog']
        corr = self.cat_corr['correction']

        # get data file to construct nbar(z) file name 
        data_file = fc_data.get_galaxy_data_file('data', **cat_corr)
        data_dir = '/'.join(data_file.split('/')[:-1])+'/'      # directory
        data_file_name = data_file.split('/')[-1]           # file name 

        # correction specifier
        if 'cmass' in catalog['name'].lower(): 
            corr_str = '.'.join(data_file_name.split('.')[1:])
        else: 
            raise NameError('not yet coded!')

        # combine to form filename  
        file_name = ''.join([data_dir, 'NBAR-', catalog['name'].lower(), '.', corr_str])

        return file_name
    
    def calculate(self, **kwargs): 
        ''' Calculate nbar(z)
        '''

        if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo -----------------------------------------------

            self.file_name = get_nbar_file(**cat_corr)  # nbar file name 
            
            ldg_nbar = 0.0000944233     # constant nbar(z) value for true

            if correction['name'] .lower() in ('true', 'upweight', 'peakshot'): 
                # no correction 
                self.nbar = np.array([ldg_nbar for i in range(n_zbin)]) 

            else:
                # nbar(z) needs to be corrected 

                if (os.path.isfile(self.file_name) == True) and (clobber == False): 

                    self.nbar = np.loadtxt(self.file_name, unpack=True, usecols=[3])

                else: 

                    print 'Constructing ', self.file_name 

                    # check that the corrected nbar-ngal files exist
                    corr_ngal_file = get_nbar_ngal_file('allmocks', **cat_corr)

                    if (os.path.isfile(corr_ngal_file) == False) or clobber: 
                        print 'Constructing ', corr_ngal_file
                        write_nbar_ngal('allmocks', **cat_corr) 
                    else: 
                        print 'Using ', corr_ngal_file

                    # read corrected nbar_ngal
                    corr_ngal = np.loadtxt(corr_ngal_file, unpack=True, usecols=[3]) 
       
                    # check that the true nbar_ngal fiel exists
                    true_cat_corr = {'catalog':catalog, 'correction':{'name':'true'}}
                    true_ngal_file = get_nbar_ngal_file('allmocks', **true_cat_corr) 

                    if os.path.isfile(true_ngal_file) == False: 
                        print 'Constructing ', true_ngal_file
                        write_nbar_ngal('allmocks', **true_cat_corr) 
                    
                    # read true nbar_ngal
                    true_ngal = np.loadtxt(true_ngal_file, unpack=True, usecols=[3])

                    # determine corrected nbar(z) using corrected nbar_ngal file
                    self.nbar = np.zeros(len(self.zlow))
                    for i in range(len(true_ngal)): 
                        if true_ngal[i] != 0: 
                            self.nbar[i] = ldg_nbar*(corr_ngal[i]/true_ngal[i])
                            #print (corr_ngal[i]/true_ngal[i])
                        else: 
                            self.nbar[i] = 0.0

                    self.writenbar()

        else: 
            raise NameError("Not yet Coded!") 

    def writenbar(self): 
        ''' Write class object nbar to ASCII file 
        '''
        np.savetxt(self.file_name, 
                np.c_[self.zmid, self.zlow, self.zhigh, self.nbar], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e'], delimiter='\t')

class Ngal: 
    def __init__(self, cat_corr, **kwargs): 
        ''' Object to handle Ngal in redshift bins 

        Notes
        -----
        * Simple to handle nbar(z) ratios in a cosmology independent way 

        '''
        self.cat_corr = cat_corr    # save catalog correction dictionary 
        catalog = cat_corr['catalog'] 
        correction  = cat_corr['correction'] 

        # set zlow, zmid, zhigh (consistent with CMASS nbar(z) files) 
        self.zlow = np.arange(0.0, 1.005, 0.005)
        self.zhigh = self.zlow+0.005 
        self.zmid = self.zlow+0.0025 
        n_zbin = len(self.zlow) 

        self.ngal = None
        self.file_name = self.File(**kwargs)

    def File(self, **kwargs): 
        ''' Get file name of nbar(z) file 
        '''
        cat = self.cat_corr['catalog']
        corr = self.cat_corr['correction']

        # get data file to construct nbar(z) file name 
        data_file = fc_data.get_galaxy_data_file('data', **self.cat_corr)
        data_dir = '/'.join(data_file.split('/')[:-1])+'/'      # directory
        data_file_name = data_file.split('/')[-1]           # file name 

        # correction specifier
        if 'cmass' in cat['name'].lower(): 
            corr_str = '.'.join(data_file_name.split('.')[1:])
        else: 
            raise NameError('not yet coded!')

        # combine to form filename  
        file_name = ''.join([data_dir, 'NGAL-', cat['name'].lower(), '.', corr_str])

        return file_name
    
    def calculate(self, **kwargs): 
        ''' Calculate Ngal(z) from the data file specified by the dictionary   
        '''
        cat = self.cat_corr['catalog'] 

        # read data from cat_corr dictionary 
        data = fc_data.galaxy_data('data', **self.cat_corr)

        redshifts = data.z  # redshifts 
        
        if 'noweights' in kwargs.keys(): 
            if kwargs['noweights']: 
                weights = np.array([1.0 for i in range(len(redshifts))])
        else:
            if 'cmass' in cat['name'].lower(): 
                # weights and completeness accounted for CMASS like catalogs 
                weights = data.wsys * (data.wnoz + data.wfc - 1.0) * (1./data.comp) 
            else: 
                raise NotImplementedError("Only CMASS implemented so far") 
    
        ngal = [] 
        # for each redshift bin calculate the weighted Ngal values
        for zlow, zhigh in zip(self.zlow, self.zhigh): 
            zbin = np.where((redshifts >= zlow) & (redshifts < zhigh))
            
            ngal.append( np.sum(weights[zbin]) ) 

            print len(weights[zbin]), np.sum(weights[zbin])

        self.ngal = np.array(ngal)

        return self.ngal 

    def Write(self, **kwargs): 
        ''' Write Ngal to file  
        '''
        # check that there are ngal files 
        if not isinstance(self.ngal, np.ndarray): 
            raise ValueError("Ngal object does not have ngal values") 
        
        # header is hardcoded
        head_str = "# columns : zmid, zlow, zhigh, Ngal "
        # write file 
        np.savetxt(self.file_name,            
                np.c_[
                    self.zmid, self.zlow, self.zhigh, self.ngal
                    ],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t', header=head_str)

        return self.file_name

    def Read(self, **kwargs): 
        ''' Read Ngal file (make sure columns are synchronized with Write function) 
        '''
        if not os.path.isfile(self.file_name): 
            self.calculate(**kwargs)
            self.Write(**kwargs)
        if 'clobber' in kwargs.keys():
            if kwargs['clobber']: 
                self.calculate(**kwargs)
                self.Write(**kwargs)

        self.zmid, self.zlow, self.zhigh, self.ngal = np.loadtxt(self.file_name, 
                skiprows=1, unpack=True, usecols=[0,1,2,3])
        return self.file_name 

# Plotting -----
def plot_Ngal(cat_corr, type='regular', **kwargs): 
    ''' Plot Ngal(z) given catalog and correction dictionary

    Parameters
    ----------
    cat_corr : (list of) dictionary specifying catalog and correction  
    type : 'regular', 'ratio' type of plot 

    '''
    # list of catalog and correction methods 
    if isinstance(cat_corr, list): 
        cat_corrs = cat_corr 
    else: 
        cat_corrs = [cat_corr]
    cat_corr_str = ''

    prettyplot()        # set up plot 
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(8,8)) 
    sub = fig.add_subplot(111) 

    for i_cc, cc in enumerate(cat_corrs):    # loop through cat corrs 
        cat = cc['catalog']
        corr = cc['correction']
        # read Ngal 
        ng = Ngal(cc, **kwargs)
        ng.Read(**kwargs)  

        # catalog correction specifier 
        cat_corr_label = ' '.join([cat['name'], corr['name']]) 
        cat_corr_str += ''.join(['_', catalog['name'], '_', correction['name']]) 
        
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

    sub.set_xlim([zmin, zmax]) 
    sub.set_xlabel('Redshift (z)')
    sub.set_ylabel(y_label)

    sub.legend(loc = 'upper left') 
    
    fig_name = ''.join(['figure/', 
        'NGAL_', cat_corr_str, type_str, '.png']) 
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear()
    return None 

if __name__=='__main__':
    cat_corrs = [
            {'catalog': {'name': 'cmasslowz_high'}, 
                'correction': {'name': 'upweight'}}, 
            {'catalog': {'name': 'cmasslowz_high'}, 
                'correction': {'name': 'peakshot', 'fit': 'gauss', 'sigma': 6.6, 'fpeak': 0.72}}
            ]
    plot_Ngal(cat_corrs)
