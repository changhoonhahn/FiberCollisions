'''


Line of Sight Displacement Codes

Author(s): ChangHoon Hahn 


'''

import numpy as np
from scipy.optimize import curve_fit
import sys
import os.path
import subprocess
import cosmolopy as cosmos
import multiprocessing as mp

# -- Local -- 
import fibcol_data as fc_data
import fibcol_utility as fc_util
import mpfit as mpfit
import galaxy_environment as genv
import pyspherematch as pysph

def photoz_dlos_nseries(n_mock, **cat_corr): 
    ''' Calculate dLOS using photometric redshifts; hardcoded for Nseries mocks 

    Parameters
    ----------
    * cat_corr : catalog and correction dictionary 

    Notes
    -----

    '''
    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction']

    if catalog['name'].lower() != 'nseries':
        raise NotImplementedError('Only for Nseries') 
    correction['name'] = 'photoz'
   
    for i_mock in range(1, n_mock+1): 
        catalog['n_mock'] = i_mock 
        cat_corr = {
                'catalog': catalog, 
                'correction': correction 
                } 
        data = fc_data.galaxy_data('data', **cat_corr)
        print data.file_name
        cosmo = data.cosmo      # survey cosmology 

        fcoll = np.where(data.wfc == 0) # fiber collided
        
        try: 
            tot_zupw = np.append(tot_zupw, data.zupw[fcoll]) 
            tot_z = np.append(tot_z, data.z[fcoll])
            tot_zphoto = np.append(tot_zphoto, data.z_photo[fcoll])
        except UnboundLocalError: 
            tot_zupw = data.zupw[fcoll]
            tot_z = data.z[fcoll]
            tot_zphoto = data.z_photo[fcoll]

    # Comoving distance of upweighted galaxy
    Dc_upw = cosmos.distance.comoving_distance(
            tot_zupw, **cosmo) * cosmo['h']
    # Comoving distance of fibcollided galaxy
    Dc_zspec = cosmos.distance.comoving_distance(
            np.array(tot_z), **cosmo) * cosmo['h']
    # Comoving distance of fibcollided galaxy photoz
    Dc_zphoto = cosmos.distance.comoving_distance(
            np.array(tot_zphoto), **cosmo) * cosmo['h']
    
    LOS_d = Dc_zspec - Dc_upw
    LOS_d_photo = Dc_zphoto - Dc_upw
    
    # Determine appropriate binsize for dist (not too important)
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min, x_max, binsize = -1000.0, 1000.0, 0.5
    n_bins = int((x_max-x_min)/binsize) 
    
    # calculate dLOS distributions 
    dlos_spec_hist, mpc_binedges = np.histogram(
            LOS_d, bins=n_bins, range=[x_min, x_max])
    dlos_photo_hist, mpc_binedges = np.histogram(
            LOS_d_photo, bins=n_bins, range=[x_min, x_max])
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([
        0.5*(xlow[i] + xhigh[i]) for i in range(len(xlow))
        ])

    # For galaxies actually in the peak  
    in_peak = np.where((LOS_d < 20.) & (LOS_d > -20.))
    dlos_spec_hist_peak, mpc_binedges = np.histogram(
            LOS_d[in_peak], bins=n_bins, range=[x_min, x_max])
    dlos_photo_hist_peak, mpc_binedges = np.histogram(
            LOS_d_photo[in_peak], bins=n_bins, range=[x_min, x_max])
    
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=[14,5]) 
    sub = fig.add_subplot(111) 

    sub.plot(xmid, dlos_spec_hist, 
            lw=4, color=pretty_colors[0], label=r"$d_{LOS}$; $z_\mathtt{spec}$")
    sub.plot(xmid, dlos_photo_hist, 
            lw=4, color=pretty_colors[2], label=r"$d_{LOS}$; $z_\mathtt{photo}$")

    sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
    sub.set_xlim([-100.0, 100.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_spec_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 
    
    fig_dir = 'figure/'
    fig_file = ''.join([fig_dir, 
        'photoz_dlos_specz_dlos_nseries.png'])
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear() 

    fig = plt.figure(1, figsize=[14,5])
    sub = fig.add_subplot(111) 

    sub.plot(xmid, dlos_spec_hist_peak, 
            lw=4, color=pretty_colors[0], label=r"Peak $d_{LOS}$; $z_\mathtt{spec}$")
    sub.plot(xmid, dlos_photo_hist_peak, 
            lw=4, color=pretty_colors[2], label=r"Peak $d_{LOS}$; $z_\mathtt{photo}$")

    sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
    sub.set_xlim([-200.0, 200.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_spec_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 
    
    fig_dir = 'figure/'
    fig_file = ''.join([fig_dir, 
        'photoz_dlos_specz_dlos_peak_nseries.png'])
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear() 
    
    fig = plt.figure(1, figsize=[14,5])
    sub = fig.add_subplot(111) 

    sub.plot(xmid, dlos_photo_hist, 
            lw=4, color=pretty_colors[0], label=r"$d_{LOS}$; $z_\mathtt{photo}$")
    sub.plot(xmid, dlos_photo_hist_peak, 
            lw=4, color=pretty_colors[2], label=r"Peak $d_{LOS}$; $z_\mathtt{photo}$")

    sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
    sub.set_xlim([-1000.0, 1000.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_photo_hist_peak)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 
    
    fig_dir = 'figure/'
    fig_file = ''.join([fig_dir, 
        'photoz_dlos_peak_comparison_nseries.png'])
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear() 
    
    fig = plt.figure(1, figsize=[14,5])
    sub = fig.add_subplot(111) 
    
    dlos_photo_ratio = []
    for i in range(len(dlos_photo_hist_peak)): 
        if dlos_photo_hist[i] == 0.: 
            dlos_photo_ratio.append(0.0)
        else:  
            dlos_photo_ratio.append(
                    np.float(dlos_photo_hist_peak[i])/np.float(dlos_photo_hist[i]))

    sub.plot(xmid, dlos_photo_ratio, 
            lw=4, color=pretty_colors[0], label=r"Peak $d_{LOS}$ Ratio; $z_\mathtt{photo}$")

    sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
    sub.set_xlim([-500.0, 500.0])
    sub.set_ylim([0.0, 1.0])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 
    
    fig_dir = 'figure/'
    fig_file = ''.join([fig_dir, 
        'photoz_dlos_peak_ratio_comparison_nseries.png'])
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear() 
# ------------------------------------------------------------------------------------
# Set up LOS displacement class
# ------------------------------------------------------------------------------------
class Dlos: 
    def __init__(self, cat_corr, **kwargs):
        ''' Given catalog correction dictionary, get line-of-sight displacement 
        '''
        self.cat_corr = cat_corr    # catalog correction dictionary 
        self.dlos = None
        self.targ_z = None
        self.neigh_z = None 

        self.file_name = self.File(**kwargs)    # File Name

        if 'env' in kwargs.keys(): 
            self.env = None 

    def File(self, **kwargs): 
        ''' Get file name for line-of-sight displacement file 
        '''
        cat = self.cat_corr['catalog']
        corr = self.cat_corr['correction']

        data_dir = '/mount/riachuelo1/hahn/data/'

        if cat['name'].lower() == 'nseries': 
            # Nseries 
            data_dir += 'Nseries/'

            # file_name
            if 'combined' in kwargs.keys():
                # combined dLOS
                file = ''.join(['DLOS_CutskyN_', str(kwargs['combined']), 'mocks']) 
            else: 
                file = ''.join(['DLOS_CutskyN', str(cat['n_mock'])])
            
            # extra specifiers organized this way for easy implementation in the future 
            extra_str = ''  
            if 'env' in kwargs.keys(): 
                if isinstance(kwargs['env'], str): 
                    extra_str += ''.join(['_env', kwargs['env']])
                else: 
                    raise TypeError("env argument has to be a string") 
            
            file_name = ''.join([data_dir, file, extra_str, '.dat'])
        else: 
            raise NotImplementedError("Only Nseries implemented at the moment")

        return file_name 

    def Columns(self, **kwargs):
        ''' Provide file column order for reading file 
        '''
        cat = self.cat_corr['catalog'] 

        if cat['name'].lower() in ('nseries'): 
            # Nseries
            columns = ['dlos', 'targ_ra', 'targ_dec', 'targ_z', 
                    'neigh_ra', 'neigh_dec', 'neigh_z']
            col_index = [0,1,2,3,4,5,6]
        else:
            raise NotImplementedError("Only Nseries implemented") 

        if 'env' in kwargs.keys():
            columns.append('env') 
            col_index.append(col_index[-1]+1)

        return [columns, indices]

    def Write(self, **kwargs): 
        ''' Write out values of Dlos object to file 
        '''
        pass

class dlos:
    def __init__(self, readdata=True, clobber=False, **cat_corr): 
        '''
        Given catalog_correction info, read dLOS values. If dLOS values do not exist, make them
        '''
        catalog = cat_corr['catalog']
        correction = cat_corr['correction']

        if catalog['name'].lower() == 'lasdamasgeo':        # LasDamas Geo ------------------
            file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
            #File name 
            if correction['name'].lower() == 'bigfc': 
                file_name = ''.join([file_dir, 
                    'DLOS_sdssmock_gamma_lrgFull_zm_oriana', 
                    str(catalog['n_mock']+100)[1:3], catalog['letter'], 
                    '_no.rdcz.big_fibcoll.dat']) 
            else: 
                file_name = ''.join([file_dir, 
                    'DLOS_sdssmock_gamma_lrgFull_zm_oriana', 
                    str(catalog['n_mock']+100)[1:3], catalog['letter'], '_no.rdcz.dat']) 
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False:      # if file does not exist then
                    build_dlos(**cat_corr) 
                
                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])

                self.dlos = readin_data[0]
                self.targ_z = readin_data[1]
                self.neigh_z = readin_data[2]
        elif catalog['name'].lower() == 'ldgdownnz':        # LDG downsampled --------------- 
            file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
            if correction['name'].lower() == 'bigfc': 
                file_name = ''.join([file_dir, 
                    'DLOS_sdssmock_gamma_lrgFull_zm_oriana', 
                    str(catalog['n_mock']+100)[1:3], catalog['letter'], 
                    '_no.rdcz.big_fibcoll.down_nz.dat']) 
            else: 
                file_name = ''.join([file_dir, 
                    'DLOS_sdssmock_gamma_lrgFull_zm_oriana', 
                    str(catalog['n_mock']+100)[1:3], catalog['letter'], 
                    '_no.rdcz.down_nz.dat']) 
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False:      # if file does not exist then
                    build_dlos(**cat_corr) 
                
                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])

                self.dlos = readin_data[0]
                self.targ_z = readin_data[1]
                self.neigh_z = readin_data[2]
        elif catalog['name'].lower() == 'tilingmock':       # Tiling Mock -------------------
            file_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'        # directory
            # File name 
            file_name = file_dir+'DLOS_cmass-boss5003sector-icoll012.dat'
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False: 
                    # if file does not exist then
                    build_dlos(**cat_corr) 
                
                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])

                self.dlos = readin_data[0]
                self.targ_z = readin_data[1]
                self.neigh_z = readin_data[2]
        elif catalog['name'].lower() == 'qpm':              # QPM ---------------------------
            file_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'        # directory
            # File name 
            file_name = ''.join([file_dir, 'DLOS_a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.fibcoll.dat']) 
            self.file_name = file_name 
            
            if readdata == True: 
                if (os.path.isfile(file_name) == False) or (clobber == True): 
                    print 'Constructing ', file_name 
                    # if file does not exist then
                    build_dlos(**cat_corr) 

                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5,6]) 
                self.dlos = readin_data[0]
                self.targ_ra = readin_data[1]
                self.targ_dec = readin_data[2]
                self.targ_z = readin_data[3]
                self.neigh_ra = readin_data[4] 
                self.neigh_dec = readin_data[5] 
                self.neigh_z = readin_data[6]
        elif catalog['name'].lower() == 'nseries':          # Nseries -----------------------
            file_dir = '/mount/riachuelo1/hahn/data/Nseries/' # directory
            file_name = ''.join([file_dir, 'DLOS_CutskyN', str(catalog['n_mock']), '.dat']) 
            
            self.file_name = file_name 

            if readdata == True: 
                if not os.path.isfile(file_name) or clobber: 
                    print 'Constructing ', file_name 
                    
                    build_dlos(**cat_corr)  
                    #build_dlos_nseries_test(**cat_corr) 

                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5,6]) 
                self.dlos = readin_data[0]
                self.targ_ra = readin_data[1]
                self.targ_dec = readin_data[2]
                self.targ_z = readin_data[3]
                self.neigh_ra = readin_data[4] 
                self.neigh_dec = readin_data[5] 
                self.neigh_z = readin_data[6]
        elif catalog['name'].lower() =='patchy':            # PATCHY mocks ------------------
            file_dir = '/mount/riachuelo1/hahn/data/PATCHY/dr12/v6c/'
            file_name = ''.join([file_dir, 
                'DLOS_Patchy-Mocks-DR12CMASS-N-V6C-Portsmouth-mass_', 
                str("%04d" % catalog['n_mock']), '.vetoed.fibcoll.dat']) 
            self.file_name = file_name 

            if readdata == True: 
                if (os.path.isfile(file_name) == False) or (clobber == True): 
                    print 'Constructing ', file_name 
                    build_dlos(**cat_corr)      # build dLOS file 
                
                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5,6]) 
                self.dlos = readin_data[0]
                self.targ_ra = readin_data[1]
                self.targ_dec = readin_data[2]
                self.targ_z = readin_data[3]
                self.neigh_ra = readin_data[4] 
                self.neigh_dec = readin_data[5] 
                self.neigh_z = readin_data[6]
        elif 'bigmd' in catalog['name'].lower():            # Big MD mock -------------------
            file_dir = '/mount/riachuelo1/hahn/data/BigMD/' # directory
            if catalog['name'].lower() == 'bigmd': 
                file_name = ''.join([file_dir, 'DLOS_bigMD-cmass-dr12v4_vetoed.fibcoll.dat']) 
            elif catalog['name'].lower() == 'bigmd1':
                file_name = ''.join([file_dir, 
                    'DLOS_bigMD-cmass-dr12v4-RST-standardHAM_vetoed.fibcoll.dat']) 
            elif catalog['name'].lower() == 'bigmd2': 
                file_name = ''.join([file_dir, 
                    'DLOS_bigMD-cmass-dr12v4-RST-quadru_vetoed.fibcoll.dat']) 
            elif catalog['name'].lower() == 'bigmd3': 
                file_name = ''.join([file_dir, 
                    'DLOS_BigMD-cmass-dr12v4-RST-standHAM-Vpeak_vetoed.fibcoll.dat'])
            else: 
                NotImplementedError('lkasdfkjadsf')
            
            self.file_name = file_name 

            if readdata == True: 
                if not os.path.isfile(file_name) or clobber: 
                    print 'Constructing ', file_name 
                    
                    build_dlos(**cat_corr)  
                    #build_dlos_nseries_test(**cat_corr) 

                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5,6]) 
                self.dlos = readin_data[0]
                self.targ_ra = readin_data[1]
                self.targ_dec = readin_data[2]
                self.targ_z = readin_data[3]
                self.neigh_ra = readin_data[4] 
                self.neigh_dec = readin_data[5] 
                self.neigh_z = readin_data[6]
        elif 'cmass' in catalog['name'].lower():            # CMASS -------------------------
            file_dir = '/mount/riachuelo1/hahn/data/CMASS/'        # directory
            # File name 
            if catalog['name'].lower() == 'cmass': 
                file_name = ''.join([file_dir, 
                    'DLOS_cmass-dr12v4-N-Reid-weights-zlim.dat']) 
            elif 'cmasslowz' in catalog['name'].lower():
                file_dir += 'dr12v5/'
                cmasslowz_str = ''
                if 'e2' in catalog['name'].lower(): 
                    cmasslowz_str = 'e2'
                elif 'e3' in catalog['name'].lower(): 
                    cmasslowz_str = 'e3'
                elif 'tot' in catalog['name'].lower():
                    cmasslowz_str = 'tot'
    
                if 'high' in catalog['name'].lower(): 
                    zbin_str = '-high'
                elif '_low' in catalog['name'].lower():
                    zbin_str = '-low'

                file_name = ''.join([file_dir, 
                    'DLOS_cmasslowz', cmasslowz_str, '-dr12v5-N', zbin_str, '.dat'])
            self.file_name = file_name 
            
            if readdata: 
                if not os.path.isfile(file_name): 
                    print 'Constructing ', file_name 
                    # if file does not exist then
                    build_dlos(**cat_corr) 

                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5,6]) 
                self.dlos = readin_data[0]
                self.targ_ra = readin_data[1]
                self.targ_dec = readin_data[2]
                self.targ_z = readin_data[3]
                self.neigh_ra = readin_data[4] 
                self.neigh_dec = readin_data[5] 
                self.neigh_z = readin_data[6]
                #readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2]) 
                #self.dlos = readin_data[0]
                #self.targ_z = readin_data[1]
                #self.neigh_z = readin_data[2]

        else: 
            raise NameError('not coded yet') 

    def get_dlos_hist(self, binsize=0.1, normed=False): 
        ''' 
        Compute the histogram for dlos  class
        '''
        x_min = -1000.0
        x_max = 1000.0
        n_bins = int((x_max-x_min)/binsize) 
        dlos_hist, mpc_binedges = np.histogram(self.dlos, bins=n_bins, range=[x_min, x_max], normed=normed)
        # store histogram info 
        self.dlos_hist = dlos_hist 
        self.xlow = mpc_binedges[:-1]
        self.xhigh = mpc_binedges[1:] 
        self.xmid = np.array([0.5*(self.xlow[i]+self.xhigh[i]) for i in range(len(self.xlow))])

    def expon(self, x, sig):
        return np.max(self.dlos_hist)*np.exp(-x/sig)
    def gauss(self, x, sig):
        return np.max(self.dlos_hist)*np.exp(-0.5*x**2/sig**2)

def build_dlos(**cat_corr): 
    ''' Build dLOS of fiber collided pairs using pyspherematch and cosmolopy

    Paramters
    ---------
    cat_corr : catalog and correction dictionary that specify the mock catalog

    Notes 
    -----
    * fiber collision angular scale : 62'' or 0.01722 
    * Bug in pyspherematch.py, so now back to using IDL 

    '''
    catalog = cat_corr['catalog']

    # fiber collided mock file name 
    mock_file = fc_data.get_galaxy_data_file('data', **cat_corr)

    if not os.path.isfile(mock_file): 
        mock = fc_data.galaxy_data('data', **cat_corr) 

    dLOS = dlos(readdata=False, **cat_corr) 
    dlos_file = dLOS.file_name  # dLOS file name 
    
    # idl Command 
    idl_command = ''.join(['idl -e ', 
        '"build_fibcoll_dlos_cosmo, ', 
        "'", catalog['name'], "', '", mock_file, "', '", dlos_file, "'", '"'
        ]) 
    os.system(idl_command) 
    print idl_command
    return idl_command 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    fib_angscale = 0.01722          # fiber collision angular scale

    if correction['name'].lower() != 'upweight':    # should only run on fibercollided mocks
        correction = {'name': 'upweight'} 
    
    mock = fc_data.galaxy_data('data', **cat_corr)  # import mock file
    ngal = len(mock.ra)
    cosmo = mock.cosmo      # survey cosmology 
    
    if catalog['name'].lower() in ('qpm', 'patchy'):    # fiber collision weights 
        wfc = mock.wfc  
    else: 
        wfc = mock.weight  
    
    upw_bool = wfc > 1.0
    now_bool = wfc == 0.0

    upw_index = np.arange(0, ngal)[upw_bool] 
    now_index = np.arange(0, ngal)[now_bool]
                
    # spherematch upweighted galaxies to downweighted galaxies 
    match_upw, match_now, d = pysph.spherematch(
            mock.ra[upw_index], mock.dec[upw_index], 
            mock.ra[now_index], mock.dec[now_index], 
            tol=fib_angscale) 
    
    same_match = d < 0.0001     # same match filter 
    
    # upweighted galaxies 
    i_upw = upw_index[match_upw]
    upw_ra = mock.ra[i_upw]
    upw_dec = mock.dec[i_upw]
    upw_z = mock.z[i_upw]       
    upw_Dc = cosmos.distance.comoving_distance(upw_z, **cosmo) * cosmo['h']
    # collided galaxies 
    i_coll = now_index[match_now]
    coll_ra = mock.ra[i_coll]
    coll_dec = mock.dec[i_coll]
    coll_z = mock.z[i_coll] 
    coll_Dc = cosmos.distance.comoving_distance(coll_z, **cosmo) * cosmo['h']
    
    LOS_d = coll_Dc - upw_Dc

    import os 
    
    los_disp = dlos(readdata=False, **cat_corr)  
    dlos_file = los_disp.file_name

    np.savetxt(dlos_file, 
            np.c_[LOS_d, upw_ra, upw_dec, upw_z, 
                coll_ra, coll_dec, coll_z], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
            delimiter='\t') 
    '''

def build_dlos_wrapper(params): 
    ''' Wrapper for build_dlos
    '''
    build_dlos(**params)
    return
# ------------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------------
def build_dlos_nseries_test(**cat_corr): 
    # Nseries has fiber collided galaxy redshifts
    catalog = cat_corr['catalog']
    # import original file 
    data_dir = '/mount/riachuelo1/hahn/data/Nseries/'   # directory

    # original file 
    orig_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.rdzwc']) 
    orig_z, orig_wfc, orig_zneigh = np.loadtxt(orig_file, unpack=True, usecols=[2,4,5])

    mock = fc_data.galaxy_data('data', **cat_corr)      # import mock file 
    cosmo = mock.cosmo  # survey cosmology 

    collided = orig_wfc < 1
    
    Dc_coll = cosmos.distance.comoving_distance(orig_z[collided], **cosmo) * cosmo['h']
    Dc_upw = cosmos.distance.comoving_distance(orig_zneigh[collided], **cosmo) * cosmo['h']

    LOS_d = Dc_coll - Dc_upw
    
    los_disp = dlos(readdata=False, **cat_corr)  
    dlos_file = los_disp.file_name+'_test'

    np.savetxt(dlos_file, np.c_[LOS_d], fmt=['%10.5f'], delimiter='\t') 

def build_dlos_ldg_test(**cat_corr): 
    ''' Calculate dLOS values for LasDamasGeo using z_upw and upw_index values from 
    the assign wcp routine 

    '''
    # LasDamasGeo has fiber collided galaxy redshifts
    catalog = cat_corr['catalog']
    # import original file 
    data_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'   # directory

    # original file 
    orig_file = ''.join([data_dir, 'sdssmock_gamma_lrgFull_zm_oriana', 
        str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.fibcoll.dat']) 
    orig_z, orig_wfc, orig_zneigh = np.loadtxt(orig_file, unpack=True, usecols=[2,3,4])

    mock = fc_data.galaxy_data('data', **cat_corr)      # import mock file 
    cosmo = mock.cosmo  # survey cosmology 

    collided = orig_wfc < 1
    
    Dc_coll = cosmos.distance.comoving_distance(orig_z[collided], **cosmo) * cosmo['h']
    Dc_upw = cosmos.distance.comoving_distance(orig_zneigh[collided], **cosmo) * cosmo['h']

    LOS_d = Dc_coll - Dc_upw
    
    los_disp = dlos(readdata=False, **cat_corr)  
    dlos_file = los_disp.file_name+'_test'

    np.savetxt(dlos_file, np.c_[LOS_d], fmt=['%10.5f'], delimiter='\t') 

def nseries_idl_python_dlos_test(n_mocks, clobber=False):
    ''' Compare dLOS calculated from IDL versus Python 

    Parameters
    ----------
    n_mocks : number of mock dLOS files to include 
    
    '''
    catalog = {'name': 'nseries'} 
    correction = {'name': 'upweight'} 

    for i_mock in range(1, n_mocks+1): 
        # individual catalog_correction dictonary 
        i_catalog = catalog.copy() 
        i_catalog['n_mock'] = i_mock 
        i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

        los_disp_i = dlos(clobber=clobber, **i_cat_corr)   # import DLOS of each mock 

        #nseries_dlos = np.loadtxt(los_disp_i.file_name+'_test', unpack=True, usecols=[0])
    
        #dlos_doublecheck = np.loadtxt('/mount/riachuelo1/hahn/data/Nseries/CutskyN'+str(i_mock)+'.fibcoll.gauss.peakshot.sigma4.0.fpeak0.6_fidcosmo.dat.dlosvalues', 
        #        unpack=True, usecols=[0])
        dlos_doublecheck = np.loadtxt(
                ''.join(['/mount/riachuelo1/hahn/data/Nseries/', 
                    'CutskyN', str(i_mock), '.fibcoll.gauss.peakshot.sigma3.8.fpeak0.7_fidcosmo.dat.dlosvalues']),
                unpack=True, usecols=[0])
        #'CutskyN', str(i_mock), '.fibcoll.true.peakshot.fpeak0.7_fidcosmo.dat.dlosvalues']), 

        # Combine dLOS values 
        try: 
            combined_idl_dlos
        except NameError: 
            combined_idl_dlos = los_disp_i.dlos
            #combined_py_dlos = nseries_dlos
            combined_check_dlos = dlos_doublecheck 
        else: 
            combined_idl_dlos = np.concatenate([combined_idl_dlos, los_disp_i.dlos]) 
            #combined_py_dlos = np.concatenate([combined_py_dlos, nseries_dlos])
            combined_check_dlos = np.concatenate([combined_check_dlos, dlos_doublecheck]) 

    # Determine appropriate binsize for dist (not too important)
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -1000.0
    x_max = 1000.0
    binsize = 0.1 
    n_bins = int((x_max-x_min)/binsize) 
    
    idl_dlos_hist, mpc_binedges = np.histogram(combined_idl_dlos, bins=n_bins, range=[x_min, x_max])
    #py_dlos_hist, mpc_binedges = np.histogram(combined_py_dlos, bins=n_bins, range=[x_min, x_max])
    check_dlos_hist, mpc_binedges = np.histogram(combined_check_dlos, bins=n_bins, range=[x_min, x_max])
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
    
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=[14,5]) 
    sub = fig.add_subplot(111) 

    sub.plot(xmid, idl_dlos_hist, lw=4, color=pretty_colors[1], label=r"IDL $d_{LOS}$")
    #sub.plot(xmid, py_dlos_hist, lw=4, color=pretty_colors[3], label=r"Python $d_{LOS}$")
    sub.plot(xmid, check_dlos_hist, ls='--', lw=4, color=pretty_colors[5], label=r"Generated $d_{LOS}$")

    sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
    sub.set_xlim([-50.0, 50.0])
    sub.set_ylim([0.0, 1.25*np.max(idl_dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_file = ''.join(['figure/nseries_', str(n_mocks), 'mocks_dlos_idl_python_test.png']) 
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear() 
    plt.close(fig) 

def ldg_idl_python_dlos_test(n_mocks, clobber=False):
    ''' Compare dLOS calculated from IDL versus Python 

    Parameters
    ----------
    n_mocks : number of mock dLOS files to include 
    
    '''
    catalog = {'name': 'lasdamasgeo'} 
    correction = {'name': 'upweight'} 

    for i_mock in range(1, n_mocks+1): 
        for letter in ['a', 'b', 'c', 'd']: 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_catalog['letter'] = letter
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_disp_i = dlos(clobber=clobber, **i_cat_corr)   # import DLOS of each mock 

            ldg_dlos = np.loadtxt(los_disp_i.file_name+'_test', unpack=True, usecols=[0])
        
            #dlos_doublecheck = np.loadtxt('/mount/riachuelo1/hahn/data/Nseries/CutskyN'+str(i_mock)+'.fibcoll.gauss.peakshot.sigma4.0.fpeak0.6_fidcosmo.dat.dlosvalues', 
            #        unpack=True, usecols=[0])
            #dlos_doublecheck = np.loadtxt(
            #        ''.join(['/mount/riachuelo1/hahn/data/Nseries/', 
            #            'CutskyN', str(i_mock), '.fibcoll.gauss.peakshot.sigma3.8.fpeak0.7_fidcosmo.dat.dlosvalues']),
            #        unpack=True, usecols=[0])
            #'CutskyN', str(i_mock), '.fibcoll.true.peakshot.fpeak0.7_fidcosmo.dat.dlosvalues']), 

            # Combine dLOS values 
            try: 
                combined_idl_dlos
            except NameError: 
                combined_idl_dlos = los_disp_i.dlos
                combined_py_dlos = ldg_dlos
            else: 
                combined_idl_dlos = np.concatenate([combined_idl_dlos, los_disp_i.dlos]) 
                combined_py_dlos = np.concatenate([combined_py_dlos, ldg_dlos])

    # Determine appropriate binsize for dist (not too important)
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -1000.0
    x_max = 1000.0
    binsize = 0.1 
    n_bins = int((x_max-x_min)/binsize) 
    
    idl_dlos_hist, mpc_binedges = np.histogram(combined_idl_dlos, bins=n_bins, range=[x_min, x_max])
    py_dlos_hist, mpc_binedges = np.histogram(combined_py_dlos, bins=n_bins, range=[x_min, x_max])
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
    
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=[14,5]) 
    sub = fig.add_subplot(111) 

    sub.plot(xmid, idl_dlos_hist, lw=4, color=pretty_colors[1], label=r"IDL $d_{LOS}$")
    sub.plot(xmid, py_dlos_hist, lw=4, color=pretty_colors[3], label=r"Python $d_{LOS}$")

    sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
    sub.set_xlim([-50.0, 50.0])
    sub.set_ylim([0.0, 1.25*np.max(idl_dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_file = ''.join(['figure/ldg_', str(n_mocks), 'mocks_dlos_idl_python_test.png']) 
    fig.savefig(fig_file, bbox_inches="tight")
    fig.clear() 
    plt.close(fig) 

# ------------------------------------------------------------------------------------
# Derived quantities 
# ------------------------------------------------------------------------------------
def get_dlos_curvefit_sigma(fit='expon', binsize=0.1, **cat_corr): 
    ''' calculate sigma of the bestfit exponential/gaussian to dLOS histogram 
    dlos xrange is currently hardcoded
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    los_disp = dlos(**cat_corr) 
    los_disp.get_dlos_hist(binsize=binsize, normed=True) 
    dlos_xmid = los_disp.xmid 
    dlos_hist = los_disp.dlos_hist 
    
    if fit=='expon': 
        peak_xrange =  (dlos_xmid >= 0.0) & (dlos_xmid < 15.0)         # approximate peak xrange
        popt, pcov = curve_fit(los_disp.expon, dlos_xmid[peak_xrange], dlos_hist[peak_xrange])
    elif fit=='gauss':
        peak_xrange =  (dlos_xmid > -15.0) & (dlos_xmid < 15.0)         # approximate peak xrange
        popt, pcov = curve_fit(los_disp.gauss, dlos_xmid[peak_xrange], dlos_hist[peak_xrange])
    return popt[0]

def average_dlos_sigma(fit='expon', **cat_corr): 
    '''
    calculate average(sigma) of the bestfit exponential/gaussian of fiber-collided pair dLOS for mock catalogs
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    if catalog['name'].lower() == 'lasdamasgeo': 
        # For LasDamasGeo catalog 
        n_mock = 0
        sigma_sum = 0.0 
        for i_mock in range(1,41): 
            for letter in ['a', 'b', 'c', 'd']: 
                i_catalog = catalog.copy()
                i_catalog['n_mock'] = i_mock
                i_catalog['letter'] = letter
                i_cat_corr = {'catalog': i_catalog, 'correction': correction}  
                
                sigma_i = get_dlos_curvefit_sigma(fit=fit, **i_cat_corr) 
                
                n_mock = n_mock+1
                sigma_sum = sigma_sum + sigma_i

    return sigma_sum/np.float(n_mock) 

def average_dlos_fpeak(**cat_corr): 
    '''
    Calculate the average fpeak for certain catalog
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    f_1sigma = 0.0
    n_mock = 0
    if catalog['name'].lower() == 'lasdamasgeo': 
         for i_mock in range(1,41): 
             for letter in ['a', 'b', 'c', 'd']: 
                i_catalog = catalog.copy()
                i_catalog['n_mock'] = i_mock
                i_catalog['letter'] = letter
                i_cat_corr = {'catalog': i_catalog, 'correction': correction}  
                
                dlos_i = dlos(**i_cat_corr) 
                sigma_i = get_dlos_curvefit_sigma(**i_cat_corr) 
            
                # dlos within 1 sigma 
                dlos_i_1sigma = dlos_i.dlos[(dlos_i.dlos > -2.0*sigma_i) & (dlos_i.dlos < 2.0*sigma_i)]

                f_1sigma = f_1sigma+np.float(len(dlos_i_1sigma))/np.float(len(dlos_i.dlos))
                n_mock = n_mock+1
    else: 
        raise NameError('Not yet coded') 
    return f_1sigma/np.float(n_mock) 

def cmass_dlos_fit_jackknife(**kwargs):
    ''' Fits CMASS dLOS distribution to specified function type 
    using Freedman-Diaconis binsizes and MPfit

    Parameters
    ----------
    kwargs : 
    
    Returns 
    ------- 
    '''
    cat_corr = { 'catalog': {'name': 'cmass'}, 'correction': {'name': 'upweight'}} 
   
    # Read dLOS data 
    # import DLOS values for mock 
    los_disp = dlos(**cat_corr)
    n_dlos = len(los_disp.dlos)

    for i in range(10): 
        int(n_dlos/10)
        rand_sub = np.random.randint(len(los_disp.dlos), size=2500)
        combined_dlos = (los_disp.dlos)[rand_sub]

    # Determine appropriate binsize for dist (not too important)
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -1000.0
    x_max = 1000.0
    binsize = 0.1 
    n_bins = int((x_max-x_min)/binsize) 
    
    guess_dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    guess_xlow = mpc_binedges[:-1]
    guess_xhigh = mpc_binedges[1:] 
    guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
    # determine appropriate bin size using Freedman-Diaconis Rule 
    peak_range = (guess_xmid >= 0.0) & (guess_xmid < 15.0) 
    dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
    n_sample = dlos_cumu[-1]

    iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    print 'Freedman-Diaconis binsize = ', binsize
    #-------------------------------------------------------------------------------------
    # recompute histogram using freedman-diaconis binsize 
    n_bins = int((x_max-x_min)/binsize) 

    dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])

    # fitting ------------------------------------------------------------------------------
    #if fit.lower() == 'gauss': 
    #elif fit.lower() == 'expon': 
    #    peak_xrange = (xmid > -25.0) & (xmid < 25.0) 
    #else: 
    #    raise NameError('Error')
    peak_xrange = (xmid > -25.0) & (xmid < 25.0)            # rough peak range 

    # amplitude of dLOS distribution 
    dlos_amp = np.mean(dlos_hist[(xmid >= -1.0) & (xmid < 1.0)])
    #print 'Amplitude = ', dlos_amp, np.max(dlos_hist)
  
    # MPFIT ------------------------------------------------------------------------------------
    p0 = [dlos_amp, 5.0]            # initial guess

    fa = {'x': xmid[peak_xrange], 'y': dlos_hist[peak_xrange]}
    if fit.lower() == 'expon': 
        peak_pars = mpfit.mpfit(mpfit_peak_expon, p0, functkw=fa, nprint=0)
    elif fit.lower() == 'gauss': 
        peak_pars = mpfit.mpfit(mpfit_peak_gauss, p0, functkw=fa, nprint=0)
    else: 
        raise NameError("Fit not yet coded") 
    
    bestfit_amp = peak_pars.params[0]
    sigma = peak_pars.params[1]
    print fit.lower(), ' Best FIt Sigma = ', sigma
    
    #--------------------------------------------------------------------------------------------------------------------
    # compute fpeak 
    fpeak_xmin = -3.0*sigma
    fpeak_xmax = 3.0*sigma
    xrange = (xmid >= fpeak_xmin) & (xmid < fpeak_xmax) 
    #fpeak = (np.sum(peak(xmid[xrange], popt[0])))/np.float(np.sum(dlos_hist)) 
    if fit.lower() == 'expon': 
        fpeak = (np.sum(peak_expon(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    elif fit.lower() == 'gauss': 
        fpeak = (np.sum(peak_gauss(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    else: 
        raise NameError("Fit not yet coded") 
    
    fpeak_dist = np.float(len(combined_dlos[(combined_dlos > fpeak_xmin) & (combined_dlos < fpeak_xmax)]))/np.float(len(combined_dlos))  # computed from the distribution

    print 'fpeak from fitting = ', fpeak
    print 'fpeak from distrib = ', fpeak_dist 
    
    #fpeak_dist = np.float(len(combined_dlos[(combined_dlos > -30.0) & (combined_dlos < 30.0)]))/np.float(len(combined_dlos))  # computed from the distribution

    #print 'fpeak from distrib = ', fpeak_dist 

    if sanitycheck == True:             # plot the distributions and fits for santicy check 
        prettyplot() 
        pretty_colors = prettycolors() 
        fig = plt.figure(1, figsize=[14,5]) 
        sub = fig.add_subplot(111) 

        # fitting labels 
        if fit.lower() == 'expon': 
            fit_label = "Exponential "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
        else: 
            fit_label = "Gaussian "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
    
        #guess_scale = binsize/0.1
        #sub.plot(guess_xmid, guess_scale*guess_dlos_hist, lw=2, color=pretty_colors[-1], label=r"$d_{LOS}$ binsize $=0.1$") 

        sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[0], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$")         # plot DLOS distribution

        # plot best fit 
        if fit.lower() == 'expon': 
            sub.plot(xmid, peak_expon(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
        elif fit.lower() == 'gauss': 
            sub.plot(xmid, peak_gauss(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
        else: 
            raise NameError("Fit not yet coded") 

        # indicate the fpeak range (for sanity check purposes)
        sub.vlines(fpeak_xmin, 0.0, 10.0*np.max(dlos_hist), color='k', lw=4)
        sub.vlines(fpeak_xmax, 0.0, 10.0*np.max(dlos_hist), color='k', lw=4)

        sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
        sub.text(-1.0*sigma, 0.2*np.max(dlos_hist), r"$f_{peak,dist} = "+str(fpeak_dist)+"$") 
        sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
        sub.set_xlim([-50.0, 50.0])
        #sub.set_xlim([-5.0, 5.0])
        sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
        sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 
    
        if correction['name'].lower() == 'bigfc': 
            bigfc_flag = '_bigfc'
        else: 
            bigfc_flag = ''

        fig_dir = 'figure/'
        fig_file = ''.join([fig_dir, 
            catalog['name'].lower(), '_', str(n_mocks), 
            'mocks_combined_dlos_peakfit_', fit.lower(), bigfc_flag, '.png'])
        fig.savefig(fig_file, bbox_inches="tight")
        fig.clear() 

    return [sigma, fpeak]

def combined_dlos_angsep_fit(n_mocks, fit='gauss', sanitycheck=False, **cat_corr):
    ''' Using combined dLOS values from entire catalog, we divide the dLOS values into bins of angular separation.  Then calculate appropriate binsize using Freedman-Diaconis rule, then fits the histogram
    Returns [sigma, fpeak]
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    angsep_cut = 0.01722/2.0
   
    # Read in dLOS data ------------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ------------------------------------------------
        for i_mock in range(1,41): 
            for letter in ['a', 'b', 'c', 'd']: 
                # individual catalog_correction dictonary 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock 
                i_catalog['letter'] = letter 
                i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
                
                los_disp_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
                # Combine dLOS values 
                try: 
                    combined_dlos
                except NameError: 
                    combined_dlos = los_disp_i.dlos
                else: 
                    combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos]) 

    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ------------------------------------------------
        # import DLOS values for mock 
        los_disp = dlos(**cat_corr)

        combined_dlos = los_disp.dlos
        n_mocks = 1

    elif catalog['name'].lower() == 'qpm':                      # QPM ------------------------------------------------------------
        for i_mock in range(1,n_mocks+1): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_disp_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            
            # calculate the angular separate between the fiber-collided pairs 
            fc_angsep = fc_util.ang_sep(los_disp_i.targ_ra, los_disp_i.targ_dec, los_disp_i.neigh_ra, los_disp_i.neigh_dec)

            # Combine dLOS values 
            try: 
                combined_dlos_close
                combined_dlos_far
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i.dlos
                combined_dlos_close = (los_disp_i.dlos)[fc_angsep < angsep_cut]
                combined_dlos_far = (los_disp_i.dlos)[fc_angsep >= angsep_cut]
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos])
                combined_dlos_close = np.concatenate([combined_dlos_close, (los_disp_i.dlos)[fc_angsep < angsep_cut]]) 
                combined_dlos_far = np.concatenate([combined_dlos_far, (los_disp_i.dlos)[fc_angsep >= angsep_cut]]) 

    elif catalog['name'].lower() == 'cmass': 
        los_disp = dlos(**cat_corr)
        combined_dlos = los_disp.dlos
        n_mocks = 1

    else: 
        raise NameError("not yet coded") 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Determine appropriate binsize for dist (not too important)
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -1000.0
    x_max = 1000.0
    binsize = 0.1 
    n_bins = int((x_max-x_min)/binsize) 
    
    # calculate the histograms 
    guess_dlos_hist_close, mpc_binedges = np.histogram(combined_dlos_close, bins=n_bins, range=[x_min, x_max])
    guess_dlos_hist_far, mpc_binedges = np.histogram(combined_dlos_far, bins=n_bins, range=[x_min, x_max])

    guess_xlow = mpc_binedges[:-1]
    guess_xhigh = mpc_binedges[1:] 
    guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
    # determine appropriate bin size using Freedman-Diaconis Rule 
    #peak_range = (guess_xmid >= 0.0) & (guess_xmid < 15.0) 
    #dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
    #n_sample = dlos_cumu[-1]

    #iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    #iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    #
    #binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    binsize = 0.1
    print 'Freedman-Diaconis binsize = ', binsize
    #---------------------------------------------------------------------------------------------------------------------------------------
    # recompute histogram using freedman-diaconis binsize 
    n_bins = int((x_max-x_min)/binsize) 

    dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max], normed=True)
    dlos_hist_close, mpc_binedges = np.histogram(combined_dlos_close, bins=n_bins, range=[x_min, x_max], normed=True)
    dlos_hist_far, mpc_binedges = np.histogram(combined_dlos_far, bins=n_bins, range=[x_min, x_max], normed=True)
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])

    # fitting ------------------------------------------------------------------------------------------------------------------------------------
    #if fit.lower() == 'gauss': 
    #elif fit.lower() == 'expon': 
    #    peak_xrange = (xmid > -25.0) & (xmid < 25.0) 
    #else: 
    #    raise NameError('Error')
    #peak_xrange = (xmid > -25.0) & (xmid < 25.0)            # rough peak range 

    ## amplitude of dLOS distribution 
    #dlos_amp = np.mean(dlos_hist[(xmid >= -1.0) & (xmid < 1.0)])
    ##print 'Amplitude = ', dlos_amp, np.max(dlos_hist)
  
    ## MPFIT ------------------------------------------------------------------------------------
    #p0 = [dlos_amp, 5.0]            # initial guess

    #fa = {'x': xmid[peak_xrange], 'y': dlos_hist[peak_xrange]}
    #if fit.lower() == 'expon': 
    #    peak_pars = mpfit.mpfit(mpfit_peak_expon, p0, functkw=fa, nprint=0)
    #elif fit.lower() == 'gauss': 
    #    peak_pars = mpfit.mpfit(mpfit_peak_gauss, p0, functkw=fa, nprint=0)
    #else: 
    #    raise NameError("Fit not yet coded") 
    #
    #bestfit_amp = peak_pars.params[0]
    #sigma = peak_pars.params[1]
    #print fit.lower(), ' Best FIt Sigma = ', sigma
    #
    ##--------------------------------------------------------------------------------------------------------------------
    ## compute fpeak 
    #fpeak_xmin = -5.0*sigma
    #fpeak_xmax = 5.0*sigma
    #xrange = (xmid >= fpeak_xmin) & (xmid < fpeak_xmax) 
    ##fpeak = (np.sum(peak(xmid[xrange], popt[0])))/np.float(np.sum(dlos_hist)) 
    #if fit.lower() == 'expon': 
    #    fpeak = (np.sum(peak_expon(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    #elif fit.lower() == 'gauss': 
    #    fpeak = (np.sum(peak_gauss(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    #else: 
    #    raise NameError("Fit not yet coded") 
    #
    #fpeak_dist = np.float(len(combined_dlos[(combined_dlos > fpeak_xmin) & (combined_dlos < fpeak_xmax)]))/np.float(len(combined_dlos))  # computed from the distribution

    #print 'fpeak from fitting = ', fpeak
    #print 'fpeak from distrib = ', fpeak_dist 
    #
    #fpeak_dist = np.float(len(combined_dlos[(combined_dlos > -30.0) & (combined_dlos < 30.0)]))/np.float(len(combined_dlos))  # computed from the distribution

    #print 'fpeak from distrib = ', fpeak_dist 

    if sanitycheck == True:             # plot the distributions and fits for santicy check 
        prettyplot() 
        pretty_colors = prettycolors() 
        fig = plt.figure(1, figsize=[14,5]) 
        sub = fig.add_subplot(111) 

        ## fitting labels 
        #if fit.lower() == 'expon': 
        #    fit_label = "Exponential "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
        #else: 
        #    fit_label = "Gaussian "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
    
        #guess_scale = binsize/0.1
        #sub.plot(guess_xmid, guess_scale*guess_dlos_hist, lw=2, color=pretty_colors[-1], label=r"$d_{LOS}$ binsize $=0.1$") 

        sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[0], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$ All")         # plot DLOS distribution
        sub.plot(xmid, dlos_hist_close, lw=4, color=pretty_colors[2], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$ Ang Sep < 31")         # plot DLOS distribution
        sub.plot(xmid, dlos_hist_far, lw=4, color=pretty_colors[4], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$ Ang Sep > 31")         # plot DLOS distribution

        # plot best fit 
        #if fit.lower() == 'expon': 
        #    sub.plot(xmid, peak_expon(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
        #elif fit.lower() == 'gauss': 
        #    sub.plot(xmid, peak_gauss(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
        #else: 
        #    raise NameError("Fit not yet coded") 

        # indicate the fpeak range (for sanity check purposes)
        #sub.vlines(fpeak_xmin, 0.0, 10.0*np.max(dlos_hist), color='k', lw=4)
        #sub.vlines(fpeak_xmax, 0.0, 10.0*np.max(dlos_hist), color='k', lw=4)

        #sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
        #sub.text(-1.0*sigma, 0.2*np.max(dlos_hist), r"$f_{peak,dist} = "+str(fpeak_dist)+"$") 
        sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
        sub.set_xlim([-50.0, 50.0])
        #sub.set_xlim([-5.0, 5.0])
        sub.set_ylim([0.0, 1.25*np.max(dlos_hist_close)])
        sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

        fig_dir = fc_util.get_fig_dir() 
        fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_dlos_peakfit_'+fit.lower()+'_angsep_bins.png', bbox_inches="tight")
        fig.clear() 
    #return [sigma, fpeak]

def build_combined_dlos_dist_peak(n_mock, **cat_corr):
    ''' Combine dLOS values and calculate the normalized histogram 
    for the peak of the distribution

    Parameters
    ----------
    n_mock : number of dLOS files to combine (each is a realization)
    cat_corr : catalog correction dictionary 

    Notes
    -----
    * The output normalized histogram is used for the True Peak Shot correction method 

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    n_mocks = 0 
    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ------------------------------------------------
        for i_mock in range(1,41): 
            for letter in ['a', 'b', 'c', 'd']: 
                # individual catalog_correction dictonary 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock 
                i_catalog['letter'] = letter 
                i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
                
                los_d_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
                # Combine dLOS values 
                try: 
                    combined_dlos
                except NameError: 
                    combined_dlos = los_d_i.dlos
                else: 
                    combined_dlos = np.concatenate([combined_dlos, los_d_i.dlos]) 
                n_mocks = n_mocks+1

    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ------------------------------------------------
        # import DLOS values for mock 
        los_d_i = dlos(**cat_corr)

        combined_dlos = los_d_i.dlos
        n_mocks = 1 

    elif catalog['name'].lower() in ('qpm', 'nseries'):         # QPM, N series ------------------------------------------
        for i_mock in range(1, n_mock+1): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_d_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_d_i.dlos
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_d_i.dlos]) 
            n_mocks += 1

    elif catalog['name'].lower() == 'patchy':                   # PATCHY ------------------------------------------------------------
        for i_mock in range(1,51): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_d_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_d_i.dlos
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_d_i.dlos]) 
            n_mocks = n_mocks+1

    elif 'bigmd' in catalog['name'].lower():                    # Big MD ---------------------
        los_d_i = dlos(**cat_corr) 
        
        combined_dlos = los_d_i.dlos 

        n_mocks = 1 # only one realization of BigMD

    elif catalog['name'].lower() == 'cmass': 
        los_d_i = dlos(**cat_corr)
        combined_dlos = los_d_i.dlos

        n_mocks = 1
    else: 
        raise NameError("not yet coded") 
    
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -1000.0
    x_max = 1000.0
    binsize = 0.1 
    n_bins = int((x_max-x_min)/binsize) 
    
    guess_dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    guess_xlow = mpc_binedges[:-1]
    guess_xhigh = mpc_binedges[1:] 
    guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
    # determine appropriate bin size using Freedman-Diaconis Rule 
    peak_range = (guess_xmid >= 0.0) & (guess_xmid < 15.0)            # 5000 km/s is hardcoded for now 
    dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
    n_sample = dlos_cumu[-1]

    iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    print 'Freedman-Diaconis binsize = ', binsize

    if binsize > 0.1: 
        binsize = 0.1
   
    # recompute histogram using appropriate binsize 
    # get normalized histogram ONLY near the peak
    new_x_min = -15.0
    new_x_max = 15.0
    narrow_range = (combined_dlos > new_x_min) & (combined_dlos < new_x_max)
    n_bins = int((new_x_max-new_x_min)/binsize) 

    print 'fpeak = ', np.float(len(combined_dlos[narrow_range]))/np.float(len(combined_dlos))

    dlos_hist, mpc_binedges = np.histogram(combined_dlos[narrow_range], bins=n_bins, range=[new_x_min, new_x_max], density=True)
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
    
    dlos_comb_peak_dist_file = ((los_d_i.file_name).rsplit('/',1))[0]+'/DLOS_norm_peak_dist_'+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined.dat'
    np.savetxt(dlos_comb_peak_dist_file,            
            np.c_[xmid, dlos_hist], fmt=['%10.5f', '%10.10f'], delimiter='\t')

def dlos_hist_binsize(**cat_corr): 
    '''
    Use Freedman-Diaconis rule to get the binsize for the peak distribution  
    '''
    # get line of sight dLOS histogram 
    LOSdisp = dlos(**cat_corr) 
    LOSdisp.get_dlos_hist(binsize=0.1) 
    LOSdisp_hist = LOSdisp.dlos_hist
    LOSdisp_xmid = LOSdisp.xmid
    
    # only consider what is "roughly" the peak (shouldn't matter too much because it tapers off fast) 
    peak_range = (LOSdisp_xmid >= 0.0) & (LOSdisp_xmid < 15.0) 
    LOSdisp_cumu = (LOSdisp_hist[peak_range]).cumsum() 

    # Need to determine the IQR 
    n_sample = LOSdisp_cumu[-1]
    iqr_index = fc_util.find_nearest(LOSdisp_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(LOSdisp_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)
    return binsize 

def qpm_combined_dlos_zreal():
    ''' dLOS(z_real) distribution '''
    catalog = {'name': 'qpm'} 
    correction = {'name': 'true'} 
    
    n_mocks = 0 
    for i_mock in range(1, 11): 
        # individual catalog_correction dictonary 
        i_catalog = catalog.copy() 
        i_catalog['n_mock'] = i_mock 
        i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

        los_disp_rsd_i = dlos(readdata=False, **i_cat_corr)         # import DLOS values from each mock 
        los_disp_i_file = '.zreal.'.join((los_disp_rsd_i.file_name).rsplit('.', 1)) 

        if os.path.isfile(los_disp_i_file) == False: 
            build_dlos_zreal(**i_cat_corr) 
        
        los_disp_i = np.loadtxt(los_disp_i_file, unpack=True, usecols=[0]) 

        # Combine dLOS values 
        try: 
            combined_dlos
        except NameError: 
            combined_dlos = los_disp_i
        else: 
            combined_dlos = np.concatenate([combined_dlos, los_disp_i]) 
        n_mocks = n_mocks+1

    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -1000.0
    x_max = 1000.0
    binsize = 0.1 
    n_bins = int((x_max-x_min)/binsize) 
    
    guess_dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    guess_xlow = mpc_binedges[:-1]
    guess_xhigh = mpc_binedges[1:] 
    guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
    # determine appropriate bin size using Freedman-Diaconis Rule 
    peak_range = (guess_xmid >= 0.0) & (guess_xmid < 15.0) 
    dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
    n_sample = dlos_cumu[-1]

    iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    print 'Freedman-Diaconis binsize = ', binsize
   
    # recompute histogram using appropriate binsize 
    n_bins = int((x_max-x_min)/binsize) 
    dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    print 'binsize ', binsize
    print n_bins
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
    # amplitude of dLOS distribution 
    dlos_amp = np.mean(dlos_hist[(xmid >= -1.0) & (xmid < 1.0)])
    print dlos_amp, np.max(dlos_hist)
   
    '''
    if fit.lower() == 'gauss': 
        peak_xrange = (xmid > -10.0) & (xmid < 10.0) 
    elif fit.lower() == 'expon': 
        peak_xrange = (xmid > -30.0) & (xmid < 30.0) 
    else: 
        raise NameError('Error')

    # UPDATED TO USE MPFIT
    # lambda functions for the fit  
    if fit.lower() == 'expon': 
        peak = lambda x, sig: dlos_amp*np.exp(-np.abs(x)/sig)
    elif fit.lower() == 'gauss': 
        peak = lambda x, sig: dlos_amp*np.exp(-0.5*x**2/sig**2)
    else: 
        raise NameError("Fit not yet coded") 
    
    # fit to dLOS histogram 
    popt, pcov = curve_fit(peak, xmid[peak_xrange], dlos_hist[peak_xrange]) #(updated to used MPFIT) 
    print 'SIGMA = ', popt[0]
    p0 = [dlos_amp, 5.0]            # initial guess

    fa = {'x': xmid[peak_xrange], 'y': dlos_hist[peak_xrange]}
    if fit.lower() == 'expon': 
        peak_pars = mpfit.mpfit(mpfit_peak_expon, p0, functkw=fa, nprint=0)
    elif fit.lower() == 'gauss': 
        peak_pars = mpfit.mpfit(mpfit_peak_gauss, p0, functkw=fa, nprint=0)
    else: 
        raise NameError("Fit not yet coded") 
    
    bestfit_amp = peak_pars.params[0]
    sigma = peak_pars.params[1]

    # compute fpeak 
    xrange = (xmid >= -15.0) & (xmid < 15.0) 
    #fpeak = (np.sum(peak(xmid[xrange], popt[0])))/np.float(np.sum(dlos_hist)) 
    if fit.lower() == 'expon': 
        fpeak = (np.sum(peak_expon(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    elif fit.lower() == 'gauss': 
        fpeak = (np.sum(peak_gauss(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    else: 
        raise NameError("Fit not yet coded") 
    
    # fitting labels 
    if fit.lower() == 'expon': 
        fit_label = "Exponential "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
    else: 
        fit_label = "Gaussian "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
    '''
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    guess_scale = binsize/0.1
    sub.plot(guess_xmid, guess_scale*guess_dlos_hist, lw=2, color=pretty_colors[-1], label=r"$d_{LOS}$ binsize $=0.1$") 
    sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[0], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$") 

    '''
    if fit.lower() == 'expon': 
        sub.plot(xmid, peak_expon(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
    elif fit.lower() == 'gauss': 
        sub.plot(xmid, peak_gauss(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
    else: 
        raise NameError("Fit not yet coded") 
    '''

    #sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
    sub.set_xlabel(r"$d_{LOS}$", fontsize=20) 
    sub.set_xlim([-20.0, 20.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_dlos_zreal.png', bbox_inches="tight")
    fig.clear() 

def qpm_dlos_zbins(): 
    '''
    compare QPM dlos distribution over three redshift bins in order to determine if there are redshift dependence issues 
    everything hardcoded
    '''
    # set up figure 
    prettyplot() 
    pretty_colors = prettycolors() 

    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 
    
    # for low, mid, high zbin 
    for i_zbin, zbin in enumerate(['low', 'mid', 'high']): 

        # combine 10 mock dlos files
        n_mocks = 0 
        for i_file in range(1,11): 
            file_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'        # directory
            dlos_file_name = ''.join([file_dir, 'dLOS_a0.6452_', str("%04d" % i_file), '.dr12d_cmass_ngc.vetoed.fibcoll.', zbin,'.dat']) 
            los_disp_i = np.loadtxt(dlos_file_name)

            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i]) 
            n_mocks = n_mocks+1

        # Create histogram for combined dLOS values  (binsize is just guessed)
        x_min = -1000.0
        x_max = 1000.0
        binsize = 0.1 
        n_bins = int((x_max-x_min)/binsize) 
        
        guess_dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
        guess_xlow = mpc_binedges[:-1]
        guess_xhigh = mpc_binedges[1:] 
        guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
        # determine appropriate bin size using Freedman-Diaconis Rule 
        peak_range = (guess_xmid >= 0.0) & (guess_xmid < 15.0) 
        dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
        n_sample = dlos_cumu[-1]

        iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
        iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
        
        binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
        print 'Freedman-Diaconis binsize = ', binsize
       
        # recompute histogram using appropriate binsize 
        n_bins = int((x_max-x_min)/binsize) 
        dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max], normed=True)
        print 'binsize ', binsize
        print n_bins
        xlow = mpc_binedges[:-1]
        xhigh = mpc_binedges[1:] 
        xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
        # amplitude of dLOS distribution 

        sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[i_zbin], label=r"$d_{LOS}$ "+zbin+" zbin") 

    #sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
    sub.set_xlabel(r"$d_{LOS}$", fontsize=20) 
    sub.set_xlim([-20.0, 20.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+'qpm_10mocks_combined_dlos_zbin_comparison.png', bbox_inches="tight")
    fig.clear() 

# ------------------------------------------------------------------------------------
# Plotting 
# ------------------------------------------------------------------------------------
def plot_dlos_peak(**cat_corr): 
    '''
    plot peak of dLOS distribution (testing appropriate binsize) 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    prettyplot() 
    pretty_colors = prettycolors()

    LOSdisp = dlos(**cat_corr) 
    LOSdisp.get_dlos_hist(binsize=0.25, normed=True) 
    sigma = get_dlos_curvefit_sigma(**cat_corr) 

    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    sub.plot(LOSdisp.xmid, LOSdisp.dlos_hist, color=pretty_colors[1]) 
    sub.plot(LOSdisp.xmid, LOSdisp.expon(LOSdisp.xmid, sigma), color=pretty_colors[5])
    sub.set_xlim([0.0, 5.0*sigma])  
    sub.set_ylim([0.0, 0.1])
    plt.show() 

def plot_dlos_tail(**cat_corr): 
    '''
    plot tail of dLOS distribution 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    prettyplot() 
    pretty_colors = prettycolors()

    LOSdisp = dlos(**cat_corr) 
    LOSdisp.get_dlos_hist(binsize=1.0, normed=True) 
    sigma = get_dlos_curvefit_sigma(**cat_corr) 

    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    sub.plot(LOSdisp.xmid, LOSdisp.dlos_hist, color=pretty_colors[1]) 
    sub.plot(LOSdisp.xmid, LOSdisp.expon(LOSdisp.xmid, sigma), color=pretty_colors[5])
    sub.set_xlim([0.0, 2.0*sigma])  
    sub.set_ylim([0.0, 0.2])
    plt.show() 

def plot_fcpaper_dlos(cat_corrs, **kwargs): 
    '''
    dLOS distribution of mocks for fc paper 
    Takes in multiple different files 
    '''
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=[14,5]) 
    sub = fig.add_subplot(111) 

    for icatcorr, cat_corr in enumerate(cat_corrs): 
        catalog = cat_corr['catalog']
        correction = cat_corr['correction']
        
        # number of mocks are hardcoded for the different mock catalogs 
        if catalog['name'].lower() == 'tilingmock': 
            mocks = ['']
        elif catalog['name'].lower() == 'lasdamasgeo': 
            mocks = [str(num)+letter for num in range(1,41) for letter in ['a', 'b', 'c', 'd']]
        elif catalog['name'].lower() == 'qpm': 
            mocks = range(1,10) 
        elif catalog['name'].lower() == 'nseries': 
            mocks = range(1,10)
        elif catalog['name'].lower() == 'patchy': 
            mocks = range(1,11) 
        elif catalog['name'].lower() == 'cmass': 
            mocks = ['']
        elif 'bigmd' in catalog['name'].lower(): 
            mocks = ['']
        else: 
            raise NameError('error')
    
        n_mocks = 0
        for i_mock in mocks: 
            i_catalog = catalog.copy()
            if catalog['name'].lower() != 'lasdamasgeo': 
                i_catalog['n_mock'] = i_mock 
            else: 
                i_catalog['n_mock'] = int(''.join(list(i_mock)[0:-1]))
                i_catalog['letter'] = list(i_mock)[-1]

            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
            
            los_disp_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i.dlos
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos]) 
            n_mocks = n_mocks+1

        # Create histogram for combined dLOS values  (binsize is just guessed)
        x_min = -1000.0
        x_max = 1000.0
        binsize = 0.2 
        n_bins = int((x_max-x_min)/binsize) 
    
        dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max], normed=True)
        xlow = mpc_binedges[:-1]
        xhigh = mpc_binedges[1:] 
        xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])

        lwid = 4
        if catalog['name'].lower() == 'lasdamasgeo': 
            cat_label = 'Las Damas'
            cat_color = pretty_colors[1]
        elif catalog['name'].lower() == 'qpm': 
            cat_label = 'QPM' 
            cat_color = pretty_colors[3]
        elif catalog['name'].lower() == 'nseries': 
            cat_label = 'N Series' 
            cat_color = pretty_colors[9]
            lwid=2
        elif catalog['name'].lower() == 'patchy': 
            cat_label = 'PATCHY v6c' 
            cat_color = pretty_colors[7]
        elif catalog['name'].lower() == 'tilingmock': 
            cat_label = 'Tiling Mock' 
            cat_color = pretty_colors[5]
        elif catalog['name'].lower() == 'cmass': 
            cat_label = 'BOSS CMASS'
            cat_color = pretty_colors[0]
            lwid = 2
        elif 'bigmd' in catalog['name'].lower(): 
            if catalog['name'].lower() == 'bigmd': 
                cat_label = 'Big MultiDark'
                cat_color = pretty_colors[11]
            elif catalog['name'].lower() == 'bigmd1': 
                cat_label = 'Big MultiDark RST StandardHAM'
                cat_color = pretty_colors[12]
            elif catalog['name'].lower() == 'bigmd2': 
                cat_label = 'Big MultiDark RST Quadru'
                cat_color = pretty_colors[13]
            elif catalog['name'].lower() == 'bigmd3': 
                cat_label = 'Big MultiDark RST StandardHAM Vpeak'
                cat_color = pretty_colors[14]
            else: 
                raise NameError('asdlkfjasdf') 
        else: 
            raise NameError('asdf') 

        sub.plot(xmid, dlos_hist, lw=lwid, color=cat_color, label=cat_label) 

        del combined_dlos

    sub.set_xlabel(r"d$_{\rm{LOS}}$ (Mpc)", fontsize=24) 
    sub.set_xlim([-45.0, 45.0])
    #sub.set_xlim([0.1, 100.0])
    #sub.set_ylim([0.0, 0.125])
    sub.set_ylim([0.0, 0.1])
    #sub.set_xscale('log') 
    #sub.set_yscale("log") 
    sub.legend(loc='upper left', scatterpoints=1, prop={'size':20}) 

    fig.savefig('figure/fcpaper_dlos_dist.png', bbox_inches="tight")
    fig.clear() 

def combined_catalog_dlos_fits(catalog, n_mock): 
    ''' Fit the n combined dLOS of specified catalog 

    Paramters
    ---------
    catalog : 'lasdamasgeo', 'qpm', 'patchy', 'tilingmock'
    n_mock : number of mocks 

    '''
    if 'lasdamasgeo' in catalog:  
        print 'LASDAMASGEO------------------------------------------------------'
        cat_corr = {'catalog': {'name':'lasdamasgeo'}, 'correction': {'name': 'bigfc'}}
        print 'Gauss ', combined_dlos_fit(n_mock, fit='gauss', sanitycheck=True,  **cat_corr) 
        cat_corr = {'catalog': {'name':'lasdamasgeo'}, 'correction': {'name': 'upweight'}}
        print 'Gauss ', combined_dlos_fit(n_mock, fit='gauss', sanitycheck=True,  **cat_corr) 
    
    elif 'ldgdownnz' in catalog:  
        print 'Downsampled LDG------------------------------------------------------'
        cat_corr = {'catalog': {'name':'ldgdownnz'}, 'correction': {'name': 'upweight'}}
        print 'Gauss ', combined_dlos_fit(n_mock, fit='gauss', sanitycheck=True,  **cat_corr) 
        cat_corr = {'catalog': {'name':'ldgdownnz'}, 'correction': {'name': 'bigfc'}}
        print 'Gauss ', combined_dlos_fit(n_mock, fit='gauss', sanitycheck=True,  **cat_corr) 

    elif 'qpm' in catalog: 

        print 'QPM------------------------------------------------------'
        cat_corr = {'catalog': {'name':'qpm'}, 'correction': {'name': 'upweight'}}
        #print 'Expon ', combined_dlos_fit(10, fit='expon', sanitycheck=True,  **cat_corr) 
        print 'Gauss ', combined_dlos_fit(10, fit='gauss', sanitycheck=True, 
                clobber=True, **cat_corr) 

    elif 'nseries' in catalog: 
        print 'Nseries ------------------------------------------------------'
        cat_corr = {'catalog': {'name':'nseries'}, 'correction': {'name': 'upweight'}}
        print 'Gauss ', combined_dlos_fit(n_mock, fit='gauss', sanitycheck=True, clobber=True, **cat_corr) 

    elif 'patchy' in catalog: 
        print 'PATCHY ------------------------------------------------------'
        cat_corr = {'catalog': {'name':'patchy'}, 'correction': {'name': 'upweight'}}
        #print 'Expon ', combined_dlos_fit(10, fit='expon', sanitycheck=True,  **cat_corr) 
        print 'Gauss ', combined_dlos_fit(n_mock, fit='gauss', sanitycheck=True, 
                clobber=True, **cat_corr) 

    elif 'bigmd' in catalog: 
        print 'BigMD -------------------------------------------------------'
        cat_corr = {'catalog': {'name':catalog}, 'correction': {'name': 'upweight'}}
        print 'Gauss ', combined_dlos_fit(1, fit='gauss', sanitycheck=True, 
                clobber=True, **cat_corr) 
    
    elif 'cmass' in catalog: 
        print 'CMASS -------------------------------------------------------'
        cat_corr = {'catalog': {'name':catalog}, 'correction': {'name': 'upweight'}}
        print 'Gauss ', combined_dlos_fit(1, fit='gauss', sanitycheck=True, 
                clobber=True, **cat_corr) 

    else: 
        raise NameError('asdfasdfasdf')  

if __name__=="__main__": 
    '''
    cat_corr = {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}} 
    pool = mp.Pool(processes=5)
    mapfn = pool.map
    
    arglist = [ [i, cat_corr] for i in [1,2,3,4,5,10]]
    
    mapfn( dlos_env_multiprocess, [arg for arg in arglist])

    #photoz_dlos_nseries(10, **cat_corr)

    #arglist = [ {'catalog': {'name': 'nseries', 'n_mock': i_mock}, 'correction': {'name': 'upweight'}} for i_mock in range(9,85)]
    #mapfn(build_dlos_wrapper, [arg for arg in arglist])

    pool.close()
    pool.terminate()
    pool.join() 
    '''

    #combined_dlos_dist(1, **cat_corr)
    #combined_catalog_dlos_fits('lasdamasgeo', 10)
    #combined_catalog_dlos_fits('ldgdownnz', 10)
    #combined_catalog_dlos_fits('cmass', 1)
    #combined_catalog_dlos_fits('bigmd3', 1)
    
    for cat_str in ['tot']: 
        cat_corrs = [
                {'catalog': {'name': 'cmasslowz_high'+cat_str}, 'correction': {'name': 'upweight'}}, 
                {'catalog': {'name': 'cmasslowz_low'+cat_str}, 'correction': {'name': 'upweight'}}
                ]
        for cat_corr in cat_corrs: 
            combined_catalog_dlos_fits((cat_corr['catalog'])['name'], 1)

    #{'catalog': {'name': 'bigmd'}, 'correction': {'name': 'upweight'}},
    #{'catalog': {'name': 'bigmd1'}, 'correction': {'name': 'upweight'}}, 
    #{'catalog': {'name': 'bigmd2'}, 'correction': {'name': 'upweight'}}, 
    #plot_fcpaper_dlos(cat_corrs)
