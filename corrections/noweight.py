'''

Fiber collided Galaxy catalogs with NN weights. 
Aggressively not accounting for issues


'''
import numpy as np 
import matplotlib.pyplot as plt

from corrections import Corrections
from fibcollided import UpweightCorr 

class NoweightCorr(Corrections): 
    
    def __init__(self, cat_corr, **kwargs): 
        ''' Child class of the Corrections class in corrections.py
        Noweight correction. In other words, galaxies with w_fc = 0 are *not* included in 
        the sample. 

        '''
        super(NoweightCorr, self).__init__(cat_corr, **kwargs)
        self.corr_str = self.corrstr() 

    def corrstr(self): 
        ''' Specify correction string
        '''
        corr_str = '.noweight'

        return corr_str 

    def build(self): 
        '''
        '''
        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']

        # upweight corrected galaxy catalog 
        fc_cat_corr = {
                'catalog': (self.cat_corr)['catalog'], 
                'correction': {'name': 'upweight'}
                }
        fc_mock = UpweightCorr(fc_cat_corr, **self.kwargs) 
        fc_mock_file = fc_mock.file()
        fc_mock_cols = fc_mock.datacolumns()

        fc_data = np.loadtxt(
                fc_mock_file, 
                skiprows=1, 
                unpack=True, 
                usecols=range(len(fc_mock_cols))
                )

        for i_col, fc_col in enumerate(fc_mock_cols): 
            setattr(fc_mock, fc_col, fc_data[i_col])
        
        collided = np.where(fc_mock.wfc == 0)[0]   # collided
        notcoll = np.where(fc_mock.wfc > 0)     # not collided
        
        fc_mock.wfc[notcoll] = 1.       # all not collided galaxies have weight of 1 


        data_list = [] 
        for data_col in self.datacolumns(): 
            new_col = getattr(fc_mock, data_col)[notcoll]
            data_list.append(new_col)

        data_fmts = self.datacols_fmt()
        data_hdrs = self.datacols_header()

        # write to corrected data to file 
        output_file = self.file()
        np.savetxt(
                output_file, 
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmts, 
                delimiter='\t', 
                header=data_hdrs
                ) 

        return None      


def NbarTable(cat_dict): 
    ''' Given cat_corr dictionary, calculate the 'corrected' nbar table
    for noweight P(k) calculations. The nbar table has to be corrected
    because the total statistical weights of the galaxies is not conserved
    with the noweight method. Default nbar(z) values are multiplied by 
 
    f_noweight(z) = sum w_g,nw / sum w_g, w

    '''
    # load default nbar table (these are hardcoded)
    if cat_dict['name'].lower() in ['qpm', 'bigmd']:
        nb_file = ''.join([
            '/mount/riachuelo2/rs123/BOSS/QPM/cmass/', 
            'nbar-cmass-dr12v4-N-Reid-om0p31.dat'
            ])
        if cat_dict['name'].lower() == 'qpm': 
            n_mocks = 100
        else: 
            n_mocks = 1
        zmid, zlow, zhigh, nbarz = np.loadtxt(nb_file, unpack=True, skiprows=2, usecols=[0,1,2,3]) 

    elif cat_dict['name'].lower() == 'nseries': 
        nb_file = ''.join([
            '/mount/riachuelo1/hahn/data/Nseries/', 
            'nbar-nseries-fibcoll.dat'
            ])
        n_mocks = 84
        zmid, zlow, zhigh, nbarz = np.loadtxt(nb_file, unpack=True, usecols=[0,1,2,3]) 
    else: 
        raise NotImplementedError


    fc_zs, fc_ws = [], []
    now_zs, now_ws = [], []
    # first combine all the redshift and galaxy weights
    for imock in range(1,n_mocks+1): 
        
        i_cat_dict = cat_dict.copy()
        i_cat_dict['n_mock'] = imock

        # fibercollided galaxies 
        fc_cat_corr = {'catalog': i_cat_dict, 'correction': {'name': 'upweight'}}
        fc_mock = UpweightCorr(fc_cat_corr) 
        fc_mock_file = fc_mock.file()
        if imock == 1: 
            fc_mock_cols = fc_mock.datacolumns()
            fc_z_col = fc_mock_cols.index('z')
            fc_w_col = fc_mock_cols.index('wfc')
        print fc_mock_file
        fc_z, fc_w = np.loadtxt(
                fc_mock_file, 
                skiprows=1, 
                unpack=True, 
                usecols=[fc_z_col, fc_w_col]
                )
        fc_zs += list(fc_z)
        fc_ws += list(fc_w)

        # no weight galaxies
        now_cat_corr = {'catalog': i_cat_dict, 'correction': {'name': 'noweight'}}
        now_mock = NoweightCorr(fc_cat_corr) 
        now_mock_file = now_mock.file()
        if imock == 1: 
            now_mock_cols = now_mock.datacolumns()
            now_z_col = now_mock_cols.index('z')
            now_w_col = now_mock_cols.index('wfc')
        print now_mock_file
        now_z, now_w = np.loadtxt(
                now_mock_file, 
                skiprows=1, 
                unpack=True, 
                usecols=[now_z_col, now_w_col]
                )
        now_zs += list(now_z)
        now_ws += list(now_w)
    
    fc_zs = np.array(fc_zs)
    fc_ws = np.array(fc_ws)
    now_zs = np.array(now_zs)
    now_ws = np.array(now_ws)
    
    f_noweight = [] 
    for iz in range(len(zmid)): 
        fc_zbin = np.where(
                (fc_zs > zlow[iz]) &
                (fc_zs <= zhigh[iz])
                )
        now_zbin = np.where(
                (now_zs > zlow[iz]) &
                (now_zs <= zhigh[iz])
                )
        if len(fc_zbin[0]) == 0: 
            f_noweight.append(0.)
        else: 
            fc_w_zbin = np.sum(fc_ws[fc_zbin])
            now_w_zbin = np.sum(now_ws[now_zbin])
             
            f_noweight.append(now_w_zbin/fc_w_zbin)

    # scale factor
    f_noweight = np.array(f_noweight)
    # scaled nbar(z)
    now_nbarz = nbarz * f_noweight
        
    output_file = ''.join([
        '/'.join(now_mock_file.rsplit('/')[:-1]), '/', 
        'nbarz.', cat_dict['name'], '.noweight.dat'
        ])
    print output_file 

    data_list = [zmid, zlow, zhigh, now_nbarz]
    data_fmts = ['%10.5f', '%10.5f', '%10.5f', '%.5e']
    data_hdrs = ''.join(['Noweight rescaled nbar(z)', '\n' 'zmid, zlow, zhigh, nbar'])

    np.savetxt(
            output_file, 
            (np.vstack(np.array(data_list))).T, 
            fmt=data_fmts, 
            delimiter='\t', 
            header=data_hdrs
            ) 

    # also print out f_noweight 
    fnoweight_file = ''.join([
        '/'.join(now_mock_file.rsplit('/')[:-1]), '/', 
        'fnoweight_nbarz.', cat_dict['name'], '.noweight.dat'
        ])
    print fnoweight_file
    f_data_list = [zmid, zlow, zhigh, f_noweight]
    f_data_fmts = ['%10.5f', '%10.5f', '%10.5f', '%10.5f']
    f_data_hdrs = ''.join(['Noweight rescaled factor f_noweight(z)', '\n' 
        'zmid, zlow, zhigh, f_noweight'])

    np.savetxt(
            fnoweight_file, 
            (np.vstack(np.array(f_data_list))).T, 
            fmt=f_data_fmts, 
            delimiter='\t', 
            header=f_data_hdrs
            ) 
    return None

def PlotFnoweight(): 
    ''' Plot the nbar scale factor f_noweight as a function of redshift 
    '''
    fig = plt.figure()
    sub = fig.add_subplot(111)

    cat_dirs = [
            '/mount/riachuelo1/hahn/data/BigMD/', 
            '/mount/riachuelo1/hahn/data/Nseries/', 
            '/mount/riachuelo1/hahn/data/QPM/dr12d/'
            ]

    for icat, cat in enumerate(['bigmd', 'nseries', 'qpm']):  
        fnoweight_file = ''.join([
            cat_dirs[icat], 
            'fnoweight_nbarz.', cat, '.noweight', 
            '.dat'])

        z, f_noweight = np.loadtxt(fnoweight_file, unpack=True, usecols=[0, 3]) 

        sub.plot(z, f_noweight, lw=4, label=cat.upper()) 
    # x-axis 
    sub.set_xlim([0.43, 0.7]) 
    sub.set_xlabel('Redshift ($z$)', fontsize=25) 
    # y-axis
    sub.set_ylim([0.5, 1.0]) 
    sub.set_ylabel(r'$f_{noweight}$', fontsize=25)
    # legend
    sub.legend(loc='lower left') 
    plt.show() 

if __name__=='__main__': 
    PlotFnoweight()
    #for mock in ['bigmd', 'nseries', 'qpm']:
    #    cat_dict = {'name': mock} 
    #    NbarTable(cat_dict) 
